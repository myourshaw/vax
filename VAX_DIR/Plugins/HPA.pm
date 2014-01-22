=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 HPA

=head1 SYNOPSIS

 mv HPA.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin HPA

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 HPA_<tissue>_<cell_type>, ..., HPA_subcellular_location.
 All fields are uri escaped for ';='.
 
 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
  
 Requires that the VAX.pm module be in the Plugins directory

 References:
    (1) Uhlen M, Oksvold P, Fagerberg L et al.
        Towards a knowledge-based Human Protein Atlas
        Nat Biotechnol 2010;28:1248-1250

=cut

package HPA;

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
    my ($self) = @_;

    my %tissue_cell_type_columns;
    $tissue_cell_type_columns{HPA_subcellular_location} = "Subcellular localisation of proteins based on immunofluorescently stained cells (Human Protein Atlas)";
    my @new_output_cols = ();
    my $query = "CALL $vw::vw_database.hpa_tissue_cell_type_list()";
    my $qh = $vw::vw_conn->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        push @new_output_cols, 'HPA_'.unspace($row[0]).'_'.unspace($row[1]);
        $tissue_cell_type_columns{'HPA_'.unspace($row[0]).'_'.unspace($row[1])} = "Expression profile for protein in human $row[0] $row[1] based on immunohistochemisty using tissue micro arrays [level of annotated protein expression: High, Medium, Low, None; level of antibody staining: Strong, Moderate, Weak, Negative] (Human Protein Atlas)";
    }
    @new_output_cols = sort(@new_output_cols);
    push @new_output_cols, 'HPA_subcellular_location';
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return \%tissue_cell_type_columns;
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    if (defined $bvfoa_info{gene}){
        my @data;
        
        my $query = "CALL $vw::vw_database.ensg2hpa_subcellular_location('$bvfoa_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
        }
        $line_hash->{HPA_subcellular_location} = @data ? join('|', @data) : '';

        $query = "CALL $vw::vw_database.ensg2hpa_tissue('$bvfoa_info{ensg}')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{'HPA_'.unspace($row[0]).'_'.unspace($row[1])} = $row[2];
        }
    }
    return {};
}

sub unspace{
    my $s = shift;
    $s =~ s/, |[\s"]/_/g;
    return $s;
}


1;

