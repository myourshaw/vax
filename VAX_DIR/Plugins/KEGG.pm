=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 KEGG

=head1 SYNOPSIS

 mv KEGG.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin KEGG

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new column:
 KEGG_Pathway.
 
 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
 
 Requires that the VAX.pm module be in the Plugins directory
 
 References:
    (1) Nat Biotechnol 2010;28:1248-1250
        The KEGG databases at GenomeNet
        Nucleic Acids Res 2002;30:42-46

=cut

package KEGG;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info);

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
        KEGG_Pathway
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        KEGG_Pathway => "KEGG pathway(s) (| separated) for gene (KEGG)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    if (defined $bvfoa_info{gene}){
        my @data;
        
        my $query = "CALL $vw::vw_database.ensg2kegg('$bvfoa_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
        }
        $line_hash->{KEGG_Pathway} = @data ? join('|', @data) : '';
    }
    return {};
}


1;

