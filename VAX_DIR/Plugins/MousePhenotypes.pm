=head1 LICENSE

 Copyright (c) 2013 Aliz R. Rao.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Aliz R. Rao <alizrrao@gmail.com>
    
=cut

=head1 NAME

 MousePhenotypes

=head1 SYNOPSIS

 mv MousePhenotypes.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin MousePhenotypes

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new column:
 MOUSE_PHENOTYPES

 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
 
 Requires that the VAX.pm module be in the Plugins directory

=head1 PARAMETERS

=cut

package MousePhenotypes;

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
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature','Gene',];
}

sub get_header_info {
    my @new_output_cols = qw(
        MOUSE_PHENOTYPES
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);
    
return {
        MOUSE_PHENOTYPES => "Phenotypes associated with gene knockout mice (JAX MGI)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};
    
    if (defined $bvfoa_info{hgnc}){
        my @data;

        my $query = "CALL $vw::vw_database.hgnc2mouse_phenotype('$bvfoa_info{hgnc}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
        }
        $line_hash->{MOUSE_PHENOTYPES} = @data ? join('|', @data) : '';
    }
    return {};
}


1;

