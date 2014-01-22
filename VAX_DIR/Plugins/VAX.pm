=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.

 For license, please contact

   myourshaw@ucla.edu

=head1 CONTACT

 Please email comments or questions to myourshaw@ucla.edu
 
=cut

=head1 NAME

VAX - Methods used by the Ensembl Variant Effect Predictor VAX plugins

=head1 SYNOPSIS

use VAX qw([get_unique] [...]);

=head1 METHODS

=cut


package VAX;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    &get_bvfoa_info
    &get_consequence_info
    &get_unique
    &replace_str
    &trim
    &ltrim
    &rtrim
    &strip
    &lstrip
    &rstrip
    &overlaps
    &get_base
    &replace_base
    &space_to_underscore
    &complement
    &reverse_complement
    &get_taxon_name
    &is_gap
);

#put diverse bits of data from a base variation feature into a useful object
sub get_bvfoa_info{
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my $config = $self->{config};
    my $foo = $line_hash;
    my %bvfoa_info;
    
    if ($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')){
        $bvfoa_info{variation} = $bvfoa->transcript_variation;
        $bvfoa_info{feature} = $bvfoa->variation_feature;
        $bvfoa_info{transcript} = $bvfoa->transcript;
        $bvfoa_info{enst} = $bvfoa_info{transcript}->stable_id;
        $bvfoa_info{gene} = $config->{ga}->fetch_by_transcript_stable_id($bvfoa_info{enst});
        $bvfoa_info{ensg} = $bvfoa_info{transcript}->{_gene_stable_id};
        $bvfoa_info{hgnc} = $bvfoa_info{transcript}->{_gene_hgnc};
        $bvfoa_info{translation} = $bvfoa_info{transcript}->translation;
        if(defined($bvfoa_info{translation})){
            $bvfoa_info{ensp} = $bvfoa_info{translation}->{stable_id};
            $bvfoa_info{protein_sequence} = $bvfoa_info{translation}->seq();
            $bvfoa_info{protein_length} = length($bvfoa_info{protein_sequence});
            $bvfoa_info{altered_aa_start} = $bvfoa_info{variation}->{translation_start};
            $bvfoa_info{altered_aa_end} = $bvfoa_info{variation}->{translation_end};
            $bvfoa_info{altered_base_start} = $bvfoa_info{variation}->{cdna_start};
            $bvfoa_info{altered_base_end} = $bvfoa_info{variation}->{cdna_end};
            $bvfoa_info{pep_allele_string} = $bvfoa->pep_allele_string;
            if (defined($bvfoa_info{pep_allele_string})) {
                #reference is first
                if ($bvfoa_info{pep_allele_string} =~ m/([^\/]+)\/([^\/]+)/) {
                    $bvfoa_info{amino_acid_reference} = $1;
                    $bvfoa_info{amino_acid_variant} = $2;
                    $bvfoa_info{changes_protein} = ($bvfoa_info{amino_acid_reference} ne $bvfoa_info{amino_acid_variant}) ? 1 : 0;
                }
            }
        }
    }
#    Hi Michael
#
#We deprecated this function because of its poor support for other features and it's sometimes non-exact definition of mid points. Our regulation team has a far better version we hope to integrate into 75. If we pull it off you'll be able to perform this query on every feature class. 
#
#Hope this helps
#
#Andy

    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele')){
        $bvfoa_info{variation} = $bvfoa->regulatory_feature_variation;
        $bvfoa_info{feature} = $bvfoa->regulatory_feature;
        #$bvfoa_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($bvfoa->base_variation_feature_overlap->{base_variation_feature});
        #$bvfoa_info{gene} = @{$bvfoa_info{nearest_genes}}[0];
    }
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::MotifFeatureVariationAllele')){
        $bvfoa_info{variation} = $bvfoa->motif_feature_variation;
        $bvfoa_info{feature} = $bvfoa->motif_feature;
        #$bvfoa_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($bvfoa->base_variation_feature_overlap->{base_variation_feature});
        #$bvfoa_info{gene} = @{$bvfoa_info{nearest_genes}}[0];
    }
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')){
        $bvfoa_info{variation} = $bvfoa->intergenic_variation;
        $bvfoa_info{feature} = $bvfoa->variation_feature;
        #$bvfoa_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($bvfoa->base_variation_feature_overlap->{base_variation_feature});
        #$bvfoa_info{gene} = @{$bvfoa_info{nearest_genes}}[0];
    }
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele')){
        $bvfoa_info{variation} = $bvfoa->base_variation_feature_overlap;
        $bvfoa_info{feature} = $bvfoa->base_variation_feature;
        #$bvfoa_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($bvfoa->base_variation_feature_overlap->{base_variation_feature});
        #$bvfoa_info{gene} = @{$bvfoa_info{nearest_genes}}[0];
    }
    else{
        warn "unrecognized tva type $bvfoa->type";
    }
    $bvfoa_info{base_variation_feature} = $bvfoa_info{variation}->{base_variation_feature};
    $bvfoa_info{_line} = $bvfoa_info{base_variation_feature}->{_line};
    if(!defined($bvfoa_info{ensg})) {
        $bvfoa_info{ensg} = $bvfoa_info{gene} ? $bvfoa_info{gene}->stable_id : undef;
    }
    if(!defined($bvfoa_info{hgnc})){
        if ( defined($bvfoa_info{gene}) && $bvfoa_info{gene}) {
            my @entries = grep {$_->database eq 'HGNC'} @{$bvfoa_info{gene}->get_all_DBEntries()};
            if(scalar @entries) {
                $bvfoa_info{hgnc} = $entries[0]->display_id;
            }
        }
    }
    $bvfoa_info{hgnc} = undef if defined($bvfoa_info{hgnc}) && ($bvfoa_info{hgnc} eq '' || $bvfoa_info{hgnc} eq '-');
    
    $bvfoa_info{chrom} = $bvfoa_info{feature}->seq_region_name;
    $bvfoa_info{chrom_start} = $bvfoa_info{feature}->start;
    $bvfoa_info{chrom_end} = $bvfoa_info{feature}->end;
    $bvfoa_info{chrom_strand} = $bvfoa_info{feature}->strand;
    if(defined($bvfoa_info{chrom})
    && defined($bvfoa_info{chrom_start})
    && defined($bvfoa_info{chrom_end})
    && defined($bvfoa_info{chrom_strand})){
        my %genomic_coords = (
            chr    => $bvfoa_info{chrom},
            start  => $bvfoa_info{chrom_start},
            end    => $bvfoa_info{chrom_end},
            strand => $bvfoa_info{chrom_strand}
        );
        $bvfoa_info{genomic_coords} = \%genomic_coords;
    }
    if(! $bvfoa->isa('Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele')){
        $bvfoa_info{reference_allele} = $bvfoa_info{variation}->{reference_allele}->variation_feature_seq;
        $bvfoa_info{non_reference_allele} = $bvfoa->variation_feature_seq;
    }
   
    return \%bvfoa_info;
} #get_bvfoa_info

#parse variant consequences
sub get_consequence_info{
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    my $config = $self->{config};
    
    my %consequences_info;
    
    my $term_method = $config->{terms}.'_term';
    my @c_so = map {$_->SO_term} @{$bvfoa->get_all_OverlapConsequences};
    if(@c_so){
        my @consequences = map {$_->$term_method || $_->SO_term} @{$bvfoa->get_all_OverlapConsequences};
        my @consequences_so = map {$_->SO_term || $_->SO_term} @{$bvfoa->get_all_OverlapConsequences};
        my @consequences_ensembl = map {$_->display_term || $_->SO_term} @{$bvfoa->get_all_OverlapConsequences};
        my @consequences_ncbi = map {$_->NCBI_term || $_->SO_term} @{$bvfoa->get_all_OverlapConsequences};
        my @ranks = map {$_->rank} @{$bvfoa->get_all_OverlapConsequences};
        my %consequences_ranks;
        my %consequences_ranks_so;
        my %consequences_ranks_ensembl;
        my %consequences_ranks_ncbi;
        map {$consequences_ranks{$consequences[$_]} = $ranks[$_]} 0..$#ranks;
        map {$consequences_ranks_so{$consequences_so[$_]} = $ranks[$_]} 0..$#ranks;
        map {$consequences_ranks_ensembl{$consequences_ensembl[$_]} = $ranks[$_]} 0..$#ranks;
        map {$consequences_ranks_ncbi{$consequences_ncbi[$_]} = $ranks[$_]} 0..$#ranks;
        my @consequences_sorted = sort {$consequences_ranks{$a} <=> $consequences_ranks{$b}} keys %consequences_ranks;
        my @consequences_so_sorted = sort {$consequences_ranks_so{$a} <=> $consequences_ranks_so{$b}} keys %consequences_ranks_so;
        my @consequences_ensembl_sorted = sort {$consequences_ranks_ensembl{$a} <=> $consequences_ranks_ensembl{$b}} keys %consequences_ranks_ensembl;
        my @consequences_ncbi_sorted = sort {$consequences_ranks_ncbi{$a} <=> $consequences_ranks_ncbi{$b}} keys %consequences_ranks_ncbi;
        
        $consequences_info{consequences_sorted} = \@consequences_sorted;
        $consequences_info{consequence_severest} = $consequences_sorted[0];
        $consequences_info{consequence_severest_rank} = $consequences_ranks{$consequences_info{consequence_severest}};
        $consequences_info{consequences_ranks} = \%consequences_ranks;
        
        $consequences_info{consequences_so_sorted} = \@consequences_so_sorted;
        $consequences_info{consequence_so_severest} = $consequences_so_sorted[0];
        $consequences_info{consequence_so_severest_rank} = $consequences_ranks_so{$consequences_info{consequence_so_severest}};
        $consequences_info{consequences_ranks_so} = \%consequences_ranks_so;
        
        $consequences_info{consequences_ensembl_sorted} = \@consequences_ensembl_sorted;
        $consequences_info{consequence_ensembl_severest} = $consequences_ensembl_sorted[0];
        $consequences_info{consequence_ensembl_severest_rank} = $consequences_ranks_ensembl{$consequences_info{consequence_ensembl_severest}};
        $consequences_info{consequences_ranks_ensembl} = \%consequences_ranks_ensembl;
        
        $consequences_info{consequences_ncbi_sorted} = \@consequences_ncbi_sorted;
        $consequences_info{consequence_ncbi_severest} = $consequences_ncbi_sorted[0];
        $consequences_info{consequence_ncbi_severest_rank} = $consequences_ranks_so{$consequences_info{consequence_ncbi_severest}};
        $consequences_info{consequences_ranks_ncbi} = \%consequences_ranks_ncbi;
    
        return \%consequences_info;
    }
    else{
        return undef;
    }
}

#get only unique array entries
sub get_unique{
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}
#replace part of a string
sub replace_str{
    #uses zero-based indexing
    my ($string, $replacement, $start, $end) = @_;
    my $a = substr($string,0,$start);
    my $z = substr($string,$end+1);
    my $replaced = $a . $replacement . $z;
    return $a . $replacement . $z;
}
# Perl trim function to remove whitespace from the start and end of the string
sub trim($){
	my $string = shift;
	$string =~ s/^\s+|\s+$//g;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($){
	my $string = shift;
	$string =~ s/^\s+//g;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($){
	my $string = shift;
	$string =~ s/\s+$//g;
	return $string;
}

# Perl strip function to remove regex from the start and end of the string
sub strip($$){
	my ($string,$strip) = @_;
	$string =~ s/^$strip+|$strip+$//g;
	return $string;
}
# Left strip function to remove leading regex
sub lstrip($$){
	my ($string,$strip) = @_;
	$string =~ s/^$strip+//g;
	return $string;
}
# Right strip function to remove trailing regex
sub rstrip($$){
	my ($string,$strip) = @_;
	$string =~ s/$strip+$//g;
	return $string;
}

sub overlaps{
    my %args = (@_);

    if (   ( $args{s1} >= $args{s2} )
        && ( $args{e1} <= $args{e2} ) )
    {
        return 1;
    }
    elsif (( $args{s1} <= $args{s2} )
        && ( $args{e1} >= $args{e2} ) )
    {
        return 1;
    }
    elsif (
        ( $args{s1} <= $args{s2} )
        && (   ( $args{e1} <= $args{e2} )
            && ( $args{e1} >= $args{s2} ) )
        )
    {
        return 1;
    }
    elsif (
        ( $args{e1} >= $args{e2} )
        && (   ( $args{s1} >= $args{s2} )
            && ( $args{s1} <= $args{e2} ) )
        )
    {
        return 1;
    }
    return 0;
}

sub get_base{

    #first base is position '1'
    my $position = shift;
    my $sequence = shift;
    return substr( $sequence, $position - 1, 1 );
}

sub replace_base{

    #first base is position '1'
    my $position    = shift;
    my $sequence    = shift;
    my $replacement = shift;

    my $upstream = substr( $sequence, 0, $position - 1 );
    my $snp      = substr( $sequence, $position - 1, 1 );
    my $downstream = substr( $sequence, $position, length($sequence) - $position );

    my $new_sequence = $upstream . $replacement . $downstream;
    return $new_sequence;
}

sub space_to_underscore{
    my $text = shift;
    $text =~ s/\s{1,}/_/g;
    return $text;
}

sub complement {
    my $seq = shift;
    $seq =~ tr/gatcryswkmbdhvnGATCRYSWKMBDHVN/ctagyrswmkvhdbnCTAGYRSWMKVHDBN/;
    return $seq;
}
sub reverse_complement {
    my $seq = shift;
    $seq = reverse complement($seq);
    return $seq;
}

sub get_taxon_name {
    #Bio::EnsEMBL::Compara::NCBITaxon
    my $taxon = shift;
    if (defined( $taxon->binomial)) {
        return $taxon->binomial;
    }
    if ( defined($taxon->name)) {
        return $taxon->name;
    }
    return '?';
}

sub is_gap {
    my $residue = shift;
    if ( $residue =~ m/[\.\-]/ ) {
        return 1;
    }
    return 0;
}


1;
