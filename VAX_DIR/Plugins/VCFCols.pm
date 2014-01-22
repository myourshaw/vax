=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 VCFCols

=head1 SYNOPSIS

 mv VCFCols.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin VCFCols

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds all input VCF columns as the first columns of the output:
 CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, [FORMAT, GT[, GT ...]].
 
 Requires that the VAX.pm module be in the Plugins directory
 
=cut

package VCFCols;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);
use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    #my $input_file = $self->{config}->{input_file};
    #my $in_file_handle = new FileHandle;
    #if($self->{config}->{input_file} =~ /\.gz$/){
    #    $in_file_handle->open($self->{config}->{compress}." ". $self->{config}->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $self->{config}->{input_file}, "\n");
    #}
    #else {
    #    $in_file_handle->open($self->{config}->{input_file} ) or die("ERROR: Could not read from input file ", $self->{config}->{input_file}, "\n");
    #}
    #open VCF, '<', $input_file or die "Can't open $input_file ($!)";
    
    #find all the columns in this VCF file
    my $VCF = get_in_file_handle($self->{config});
    while(<$VCF>){
        chomp;
        if(/^\#CHROM\s+POS\s+ID\sREF\sALT/){
            my @vcf_cols = split("\t");
            $vcf_cols[0] = substr($vcf_cols[0],1);
            $self->{_vcf_cols} = \@vcf_cols;
            @OUTPUT_COLS = (@vcf_cols, @OUTPUT_COLS);
            last;
        }
    }
    close $VCF;

    return $self;
}

sub version {
    return '74';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my $self = shift;
    
    my %vcf_cols;
    if (defined($self->{_vcf_cols})){
        foreach my $vcf_col(@{$self->{_vcf_cols}}){
            if(uc($vcf_col) eq 'CHROM'){
                $vcf_cols{CHROM} = "Chromosome (VCF)";
            }
            elsif(uc($vcf_col) eq 'POS'){
                $vcf_cols{'POS'} = "Position (VCF)";
            }
            elsif(uc($vcf_col) eq 'ID'){
                $vcf_cols{'ID'} = "ID (VCF)";
            }
            elsif(uc($vcf_col) eq 'REF'){
                $vcf_cols{'REF'} = "Reference allele (VCF)";
            }
            elsif(uc($vcf_col) eq 'ALT'){
                $vcf_cols{'ALT'} = "Alternate allele (VCF)";
            }
            elsif(uc($vcf_col) eq 'QUAL'){
                $vcf_cols{'QUAL'} = "Quality score (VCF)";
            }
            elsif(uc($vcf_col) eq 'FILTER'){
                $vcf_cols{'FILTER'} = "Filter (VCF)";
            }
            elsif(uc($vcf_col) eq 'INFO'){
                $vcf_cols{'INFO'} = "Info (VCF)";
            }
            elsif(uc($vcf_col) eq 'FORMAT'){
                $vcf_cols{'FORMAT'} = "Format of genotype columns (VCF)";
            }
            else{
                $vcf_cols{$vcf_col} = "$vcf_col genotype (VCF)";
            }
        }
        
    }
    return \%vcf_cols;
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};
    
    #add all VCF fields to the output file
    my @vcf_data = split("\t", $input_line);
    map {$line_hash->{$self->{_vcf_cols}[$_]} = $vcf_data[$_]} 0..$#{$self->{_vcf_cols}};
    return {};
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified
    else {
        #$in_file_handle = 'STDIN';
        die("The VCFCols plugin cannot function with STDIN input (or maybe you forgot to specify an input file?)...");
    }
    
    return $in_file_handle;
}

1;

