=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 DiseasesPhenotypes

=head1 SYNOPSIS

 mv VCFCols.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin DiseasesPhenotypes

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds column:
 DISEASES_PHENOTYPES.
 This plugin just adds the column without data. The column can be populated by other plugins
 such as HGMD, OMIM, Phenotypes, UniProt.
 
=cut

package DiseasesPhenotypes;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub version {
    return '74';
}

sub feature_types {
    return [];
}

sub get_header_info {
    my $self = shift;
    
    my @new_output_cols = qw(
        DISEASES_PHENOTYPES
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        DISEASES_PHENOTYPES => "Union of disease and phenotype columns (Ensembl,HGMD,OMIM,UniProt)",
    };

}

sub run {
    return {};
}


1;

