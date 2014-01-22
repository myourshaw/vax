=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Consequences

=head1 SYNOPSIS

 mv Consequences.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Consequences

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 Consequence_severest, Consequence_rank, Consequences_all.
 
 Requires that the VAX.pm module be in the Plugins directory

=cut

package Consequences;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_consequence_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

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
    my @new_output_cols = qw(
        Consequence_severest
        Consequence_rank
        Consequences_all
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Consequence_severest => "Severest consequence to transcript of variant",
        Consequence_rank => "Rank of severest consequence (lower is more severe)",
        Consequences_all => "All consequences of variant to all transcripts",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    
    #get Consequence_severest, Consequence_rank, and Consequences_all for tis variant/transcript
    my $ci = get_consequence_info(@_);
    if(defined($ci)){
        my %consequences_info = %{$ci};
        $line_hash->{Consequence_severest} = $consequences_info{consequence_severest};
        $line_hash->{Consequence_rank} = $consequences_info{consequence_severest_rank};
        $line_hash->{Consequences_all} = join ",", @{$consequences_info{consequences_sorted}};
    }
    else{
        $line_hash->{Consequence_severest} = '';
        $line_hash->{Consequence_rank} = '';
        $line_hash->{Consequences_all} = '';
    }
    return {};
}

1;

