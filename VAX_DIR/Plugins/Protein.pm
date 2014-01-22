=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      
 This is a work in part modified from the NGS-SNP annotate_SNPs.pl script, 
 Grant JR, Arantes AS, Liao X, Stothard P., In-depth annotation of SNPs arising from resequencing projects using NGS-SNP, Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.                                                                      
 and is licensed as a whole at no charge to all third parties under the terms of the original GNU GENERAL PUBLIC LICENSE.

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Protein

=head1 SYNOPSIS

 mv Protein.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Protein

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 Protein_Length, Protein_Length_Decrease, Protein_Sequence_Lost, Protein_Length_Increase,
 Protein_Sequence_Gained, Reference_Splice_Site, Variant_Splice_Site, Overlapping_Protein_Domains.
 
 Requires that the VAX.pm module be in the Plugins directory
 
References:

 (1) Grant JR, Arantes AS, Liao X, Stothard P. 
     In-depth annotation of SNPs arising from resequencing projects using NGS-SNP
     Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.

=cut

package Protein;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_bvfoa_info get_consequence_info complement get_unique replace_base replace_str overlaps);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->{config};
    
    return $self;
}

sub version {
    return '74';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript'];
}

sub get_header_info {
    my @new_output_cols = qw(
        Protein_Length
        Protein_Length_Decrease(%)
        Protein_Sequence_Lost
        Protein_Length_Increase(%)
        Protein_Sequence_Gained
        Reference_Splice_Site
        Variant_Splice_Site
        Overlapping_Protein_Domains
    );
    
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Protein_Length => "Number of amino acids in translation of the transcript of the reference sequence",
        "Protein_Length_Decrease(%)" => "gives the length in amino acids of the protein segment that is lost because of an allele that introduces a stop codon (i.e. functional class is 'STOP_GAINED'). The value given in parentheses is the length of the lost protein segment expressed as a percentage of the length of the reference protein",
        Protein_Sequence_Lost => "the peptide sequence removed from the reference sequence because of an allele that introduces a stop codon (i.e. functional class is 'STOP_GAINED')",
        "Protein_Length_Increase(%)" => "gives the length in amino acids of the protein segment that is gained because of an allele that removes a stop codon (i.e. functional class is 'STOP_LOST'). The value given in parentheses is the length of the gained protein segment expressed as a percentage of the length of the reference protein",
        Protein_Sequence_Gained => "the peptide sequence added to the reference sequence because of an allele that removes a stop codon (i.e. functional class is 'STOP_LOST')",
        Reference_Splice_Site => "the sequence of the splice site that is altered by the SNP. The splice site bases (i.e. the first two and last two bases in the intron) are given as they appear on the transcribed strand of the reference sequence, i.e., typically GT~AG on either strand. This value is reported when the functional class is 'splice_donor_variant' or 'splice_acceptor_variant'",
        Variant_Splice_Site => "the sequence of the splice site that is altered by the SNP. The splice site bases (i.e. the first two and last two bases in the intron) are given as they appear on the transcribed strand of the variant sequence. This value is reported when the functional class is 'splice_donor_variant' or 'splice_acceptor_variant'",
        Overlapping_Protein_Domains => "protein domains that overlap with the position of the affected amino acid (Ensembl)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $ci = get_consequence_info(@_);
    my %consequences_info;
    if (defined($ci)){
        %consequences_info = %{get_consequence_info(@_)};
    }

    if(defined($bvfoa_info{translation})
       && defined($bvfoa_info{protein_sequence})
       && defined($bvfoa_info{transcript})
       ){
        $line_hash->{Protein_Length} = $bvfoa_info{protein_length};

        #protein domains (from Ensembl) overlapping with affected amino acid
        my $overlapping_protein_domains =
            determine_overlapping_protein_domains(
                $bvfoa_info{translation},
                $bvfoa_info{altered_aa_start},
                $bvfoa_info{altered_aa_end}
                );
        if ((defined( $overlapping_protein_domains))
            && (scalar(@{$overlapping_protein_domains}) > 0))
        {
            s/\;//g for (@{$overlapping_protein_domains});
            $line_hash->{Overlapping_Protein_Domains} =
                join( ';', @{$overlapping_protein_domains});
        }

        if(%consequences_info){
            #determine how STOP_GAINED SNP changes protein length
            #this may need to be altered for indels
            my ($protein_length_decrease, $protein_sequence_lost) =
                determine_effect_of_stop_gained_on_protein(
                    $consequences_info{consequences_ranks_so},
                    $bvfoa_info{protein_sequence},
                    $bvfoa_info{altered_aa_start}
                    );
            #Protein_Length_Decrease
            if (defined($protein_length_decrease)) {
                $line_hash->{"Protein_Length_Decrease(%)"} = $protein_length_decrease;
            }
            #Protein_Sequence_Lost
            if (defined($protein_sequence_lost)) {
                $line_hash->{Protein_Sequence_Lost} = $protein_sequence_lost;
            }
    
            #determine how STOP_LOST SNP changes protein length
            #this may need to be altered for indels
            my ($protein_sequence_gained, $protein_length_increase) =
                determine_effect_of_stop_lost_on_protein(
                    $consequences_info{consequences_ranks_so},
                    $bvfoa_info{protein_sequence},
                    $bvfoa_info{altered_base_start},
                    $bvfoa_info{transcript},
                    $bvfoa_info{amino_acid_variant} #$transcript_snp_reads
                    );
            #Protein_Length_Increase
            if (defined($protein_length_increase)) {
                $line_hash->{"Protein_Length_Increase(%)"} = $protein_length_increase;
            }
            #Protein_Sequence_Gained
            if (defined($protein_sequence_gained)) {
                $line_hash->{Protein_Sequence_Gained} = $protein_sequence_gained;
            }
        }
     } #if(defined($translation))

    if(defined($bvfoa_info{transcript})){
        if(%consequences_info){
            #determine how ESSENTIAL_SPLICE_SITE SNP changes
            #splice site
            my ($reference_splice_site, $variant_splice_site) =
                determine_effect_of_essential_splice_site_on_splice_site(
                    $consequences_info{consequences_ranks_so},
                    $bvfoa_info{transcript},
                    $bvfoa_info{genomic_coords},
                    $bvfoa_info{reference_allele}, #$Chromosome_Reference
                    $bvfoa_info{non_reference_allele}
                    );
            #Reference_Splice_Site
            if (defined($reference_splice_site)) {
                $line_hash->{Reference_Splice_Site} = $reference_splice_site;
            }
            #Variant_Splice_Site
            if (defined($variant_splice_site)) {
                $line_hash->{Variant_Splice_Site} = $variant_splice_site;
            }
        }
    }
    
    return {};
}

#Determine effect of 'ESSENTIAL_SPLICE_SITE' consequences on splice site
sub determine_effect_of_essential_splice_site_on_splice_site {
    
    my ($consequences_ranks_so, $transcript, $genomic_coords, $reference_allele, $non_reference_allele) = @_;
    
    #return values, as DD~AA
    my ($reference_splice_site, $variant_splice_site);

    #if ($consequence_ensembl eq 'ESSENTIAL_SPLICE_SITE' ) {
    if (exists $consequences_ranks_so->{splice_acceptor_variant}
        || exists $consequences_ranks_so->{splice_donor_variant}){

        my $ref = $genomic_coords->{strand}==1 ? $reference_allele : complement($reference_allele);
        my $alt = $genomic_coords->{strand}==1 ? $non_reference_allele : complement($non_reference_allele);
        $ref = $ref eq '-' ? '' : $ref;
        $alt = $alt eq '-' ? '' : $alt;

        my $introns = $transcript->get_all_Introns();

        foreach my $intron (@{$introns}) {

            $intron = $intron->transform('toplevel');

            my $intron_strand = $intron->strand();
            my $intron_start  = $intron->start();
            my $intron_end    = $intron->end();
            my $intron_length = $intron->length();
            my $intron_seq    = $intron->seq();
            
            if (($genomic_coords->{start} >= $intron_start && $genomic_coords->{start} <= $intron_end)
                || ($genomic_coords->{end} >= $intron_start && $genomic_coords->{end} <= $intron_end)){
                my $ref_seq = $intron_strand == 1 ? $intron_seq : complement($intron_seq);
                my $rel_start = $genomic_coords->{start}-$intron_start;
                my $rel_end = $genomic_coords->{end}-$intron_start;
                my $alt_seq;
                if ($ref eq '' && $genomic_coords->{start} == $genomic_coords->{end}+1){
                    #ins
                    if ($rel_start >= 0 && $rel_start <= $intron_length){
                        #insertion point adjacent to or in intron
                    if($intron_strand == -1){
                        my $rel_start_tmp = $rel_start;
                        $rel_start = $intron_length-1-$rel_end;
                        $rel_end = $intron_length-1-$rel_start_tmp;
                    }
                        $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                    }
                    else{
                        #insertion point is outside intron
                        last;
                    }
                }
                elsif ($alt eq ''){
                    #del
                    if ($rel_start < 0){
                        $rel_start = 0;
                    }
                    if ($rel_end > $intron_length-1){
                        $rel_end = $intron_length-1;
                    }
                    if($intron_strand == -1){
                        my $rel_start_tmp = $rel_start;
                        $rel_start = $intron_length-1-$rel_end;
                        $rel_end = $intron_length-1-$rel_start_tmp;
                    }
                    $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                }
                else{
                    $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                }
                if($intron_strand == -1){
                    $alt_seq = complement($alt_seq);
                }
                my $donor = substr($intron_seq,0,2);
                my $acceptor = substr($intron_seq,-2,2);
                $reference_splice_site = $donor.'~'.$acceptor;
                my $alt_donor = substr($alt_seq,0,2);
                my $alt_acceptor = substr($alt_seq,-2,2);
                $variant_splice_site = $alt_donor.'~'.$alt_acceptor;
                last;
            }
        } #foreach my $intron (@{$introns})
    } #if (exists $consequences_ranks_so->{splice_acceptor_variant} || exists $consequences_ranks_so->{splice_donor_variant})

    return ($reference_splice_site, $variant_splice_site);
}

#Determine effect of 'STOP_LOST' consequences on protein
sub determine_effect_of_stop_lost_on_protein {

    my ($consequences_ranks_so, $protein_sequence, $altered_base_start, $transcript, $transcript_snp_reads) = @_;

    #return values
    my ($protein_sequence_gained, $protein_length_increase);
    
    #NOTE: $codon_table should be set to match codon table used to
    #generate reference protein.
    my $codon_table = 1;

    #if ($consequence_ensembl eq 'STOP_LOST') {
    if (exists $consequences_ranks_so->{stop_lost}){
        if (defined($protein_sequence)) {
            my $reference_protein_length
                = length($protein_sequence);
            if ($protein_sequence =~ m/\*$/) {
                $reference_protein_length--;
            }

            my $transcript_sequence = $transcript->seq()->seq();
            if ((defined($transcript_sequence))
                && (defined($altered_base_start)))
            {

                my $translation_start = $transcript->cdna_coding_start;
#TODO: indel?
                my $variant_transcript_sequence = replace_base(
                    $altered_base_start,
                    $transcript_sequence,
                    $transcript_snp_reads
                );

                #get the variant sequence starting with the start codon
                my $variant_transcript_starting_with_start_codon
                    = substr($variant_transcript_sequence, $translation_start - 1, 1)
                    . substr($variant_transcript_sequence, $translation_start,
                    length($variant_transcript_sequence) - $translation_start);

                #create a sequence object to represent variant sequence and translate
                my $variant_transcript_starting_with_start_codon_seq_obj
                    = Bio::Seq->new(
                    -seq => $variant_transcript_starting_with_start_codon,
                    -alphabet => 'dna'
                    );
                my $variant_transcript_starting_with_start_codon_prot_obj
                    = $variant_transcript_starting_with_start_codon_seq_obj
                    ->translate(undef, undef, undef, $codon_table);
                my $variant_transcript_starting_with_start_codon_protein_sequence
                    = $variant_transcript_starting_with_start_codon_prot_obj
                    ->seq();

                #get the reference sequence starting with the start codon
                my $reference_transcript_starting_with_start_codon
                    = substr($transcript_sequence, $translation_start - 1, 1)
                    . substr( $transcript_sequence, $translation_start,
                    length($transcript_sequence) - $translation_start);

                #create a sequence object to represent reference sequence and translate
                my $reference_transcript_starting_with_start_codon_seq_obj
                    = Bio::Seq->new(
                    -seq => $reference_transcript_starting_with_start_codon,
                    -alphabet => 'dna'
                    );
                my $reference_transcript_starting_with_start_codon_prot_obj
                    = $reference_transcript_starting_with_start_codon_seq_obj
                    ->translate(undef, undef, undef, $codon_table);
                my $reference_transcript_starting_with_start_codon_protein_sequence
                    = $reference_transcript_starting_with_start_codon_prot_obj
                    ->seq();

                #need to examine translations to identify added region
                if ($reference_transcript_starting_with_start_codon_protein_sequence
                    =~ m/^([^\*]+)/)
                {
                    my $reference_cds_translation = $1;

                    if (!(length($reference_cds_translation)
                            == $reference_protein_length))
                    {
                        #something is wrong with translations--don't finish calculation
                        return;
                    }

                    if ($variant_transcript_starting_with_start_codon_protein_sequence
                        =~ m/\Q$reference_cds_translation\E([^\*]*\*?)/)
                    {
                        $protein_sequence_gained = $1;
                        $protein_length_increase = length($protein_sequence_gained);
                        if ($protein_sequence_gained =~ m/\*$/ )
                        {
                            $protein_length_increase--;
                        }

                        my $percentage_length_change = sprintf("%.0f",
                            ($protein_length_increase / $reference_protein_length) * 100);
                        $protein_length_increase =
                            $protein_length_increase . "($percentage_length_change)";
                    }
                }
            }
        }
    }
    return ($protein_sequence_gained, $protein_length_increase);
}

#Determine effect of 'STOP_GAINED' consequence on protein
sub determine_effect_of_stop_gained_on_protein {

    my ($consequences_ranks_so, $protein_sequence, $altered_aa_start) = @_;
    
    #return values
    my ($protein_length_decrease, $protein_sequence_lost);

    #if ($consequence_ensembl eq 'STOP_GAINED') {
    if (exists $consequences_ranks_so->{stop_gained}){
        if (defined($protein_sequence)) {
            my $reference_protein_length
                = length($protein_sequence);
            if ( $protein_sequence =~ m/\*$/) {
                $reference_protein_length--;
            }
            if (defined($altered_aa_start)) {
                my $length_lost = $reference_protein_length
                    - $altered_aa_start + 1;
                my $percentage_length_change = sprintf( "%.0f",
                    ($length_lost / $reference_protein_length) * 100);
                $protein_length_decrease
                    = $length_lost . "($percentage_length_change)";

                $protein_sequence_lost = substr($protein_sequence,
                    $altered_aa_start - 1, 1)
                    . substr($protein_sequence, $altered_aa_start,
                    length($protein_sequence) - $altered_aa_start);
            }
        }
    }
    return ($protein_length_decrease, $protein_sequence_lost);
}

sub determine_overlapping_protein_domains {
    my ($translation, $altered_aa_start, $altered_aa_end) = @_;
    if ((!defined($translation))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end)))
    {
        return;
    }

    my @overlapping_domains = ();

    my $pfeatures = $translation->get_all_DomainFeatures();
    while (my $pfeature = shift @{$pfeatures}) {
        my $overlaps    = 0;
        my $start       = $pfeature->start();
        my $end         = $pfeature->end();
        my $description = $pfeature->analysis()->logic_name();

        #interpro_ac() and idesc() seem to give whitespace sometimes
        if ((defined($pfeature->interpro_ac()))
            && ($pfeature->interpro_ac() =~ m/\S/))
        {
            $description = $description . ' ' . $pfeature->interpro_ac();
        }
        if ((defined($pfeature->idesc()))
            && ($pfeature->idesc() =~ m/\S/))
        {
            $description = $description . ' ' . $pfeature->idesc();
        }

        if (overlaps(
                s1 => $start,
                e1 => $end,
                s2 => $altered_aa_start,
                e2 => $altered_aa_end)
            )
        {

#   print "Protein feature '$description' starts at $start and ends at $end and overlaps with " .$altered_aa_start} . " to " .$altered_aa_end} . "\n";
            push( @overlapping_domains, $description );
        }
        else {

#   print "Protein feature '$description' starts at $start and ends at $end and does not overlap with " .$altered_aa_start} . " to " .$altered_aa_end} . "\n";
        }
    }
    

    return (scalar(@overlapping_domains) > 0) ?
        get_unique( \@overlapping_domains ) :
        undef;
}


1;

