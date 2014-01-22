=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 OMIM

=head1 SYNOPSIS

 mv OMIM.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin OMIM

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new column:
 OMIM_Disorder.
 All fields are uri escaped for ';='.
 
 Requires that the vw plugin be in the Plugins directory
  
 Requires that the VAX.pm module be in the Plugins directory

 References:
    (1) Online Mendelian Inheritance in Man OMIM®
        Online Mendelian Inheritance in Man, OMIM®
         http://omim.org/
 
=cut

package OMIM;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info get_unique);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    return $self;
}

sub version {
    return '74';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature'];
}

sub get_header_info {
    my @new_output_cols = qw(
        OMIM_Disorder(locus)
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        'OMIM_Disorder(locus)' => "OMIM disorder(locus) for gene (| separated) (OMIM)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    my @DISEASES_PHENOTYPES;
    
    #use hgnc (ensg preferred when a table is created that has it)
    if (defined $bvfoa_info{hgnc}){
        my @data;
        
        my $query = "CALL $vw::vw_database.hgnc2omim('$bvfoa_info{hgnc}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
            
            push(@DISEASES_PHENOTYPES, split(/\|/,$row[0]));
        }
        $line_hash->{'OMIM_Disorder(locus)'} = @data ? join('|', @data) : '';
    }

    #add phenotypes to the DISEASES_PHENOTYPES portmanteau column
    my @diseases_phenotypes = defined($line_hash->{DISEASES_PHENOTYPES}) ? split(/\|/,$line_hash->{DISEASES_PHENOTYPES}) : ();
    my @union_diseases_phenotypes = (@diseases_phenotypes, @DISEASES_PHENOTYPES);
    my @unique_union_diseases_phenotypes = @{get_unique(\@union_diseases_phenotypes)};
    if(@unique_union_diseases_phenotypes){
        @unique_union_diseases_phenotypes = sort(@unique_union_diseases_phenotypes);
        my $diseases_phenotypes_str = join('|',@unique_union_diseases_phenotypes);
        $line_hash->{DISEASES_PHENOTYPES} = $diseases_phenotypes_str;
    }
    
    return {};
}


1;

