=head1 LICENSE

 Copyright (c) 2011-2013 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 AlleleFrequencies

=head1 SYNOPSIS

 mv AlleleFrequencies.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin AlleleFrequencies,cortex.local,3306,vw,vw,mysql,vw

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 dbsnp137_ID,ain_dbsnp,dbsnp_pop,in_nhlbi,in_niehs,in_1kg,in_1kg_pop,G5,G5A,common,common_not_med,clinvar,CLNSIG,total_allele_count,ref_allele_count,ref_allele_frequency,alt_allele_count,alt_allele_frequency,total_sample_count,wt_sample_count,wt_frequency,het_sample_count,het_frequency,hom_sample_count,hom_frequency,gene_disruptive_wt_count,gene_disruptive_het_count,gene_disruptive_hom_count,gene_disruptive_het_freq,gene_disruptive_hom_freq,gene_missense_deleterious_wt_count,gene_missense_deleterious_het_count,gene_missense_deleterious_hom_count,gene_missense_deleterious_het_freq,gene_missense_deleterious_hom_freq,clinvar_count.
 
=head1 PARAMETERS

=cut

package AlleleFrequencies;

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
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my @new_output_cols = qw(
        dbsnp137_ID
        in_dbsnp
        in_dbsnp_pop
        in_nhlbi
        in_niehs
        in_1kg
        in_1kg_pop
        G5
        G5A
        common
        common_not_med
        clinvar
        CLNSIG
        total_allele_count
        ref_allele_count
        ref_allele_frequency
        alt_allele_count
        alt_allele_frequency
        total_sample_count
        wt_sample_count
        wt_genotype_frequency
        het_sample_count
        het_genotype_frequency
        hom_sample_count
        hom_genotype_frequency
        gene_disruptive_wt_count
        gene_disruptive_het_count
        gene_disruptive_hom_count
        gene_disruptive_het_freq
        gene_disruptive_hom_freq
        gene_missense_deleterious_wt_count
        gene_missense_deleterious_het_count
        gene_missense_deleterious_hom_count
        gene_missense_deleterious_het_freq
        gene_missense_deleterious_hom_freq
        clinvar_count
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        dbsnp137_ID => "rs ID of ALT allele (dbSNP137)",
        in_dbsnp => "ALT allele is in dbSNP137 (file 00-All.vcf.gz dated 2012-06-16, downloaded 2012-08-16) (dbSNP137)",
        in_dbsnp_pop => "ALT allele is in dbSNP137 genotyped populations (files ByPopulation/*-*-*.vcf.gz dated 2012-06-17) (dbSNP137)",
        in_nhlbi => "ALT allele is in Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) [ESP6500 accessed 2012-09-03], 6503 samples (13,006 chromosomes) (NHLBI)",
        in_niehs => "ALT allele is in NIEHS Environmental Genome Project, Seattle, WA (URL: http://evs.gs.washington.edu/niehsExome/) [v.0.0.8. (April 22, 2012) accessed 2012-09-03], 95 individuals (190 chromosomes) (NIEHS)",
        in_1kg => "ALT allele is in 1000 Genomes Phase 1 Integrated Variant Call Set release version 3.20120430 (file ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz dated 2012-04-26) (1000Genomes)",
        in_1kg_pop => "ALT allele is in genotyped populations of 1000 Genomes Phase 1 Integrated Variant Call Set release version 3.20120430 (files ALL.chr*.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz dated 2012-04-24/25) (1000Genomes)",
        G5 => ">5% minor ALT allele frequency in 1+ populations (dbSNP137)",
        G5A => ">5% minor ALT allele frequency in each and all populations (dbSNP137)",
        common => "ALT allele is polymorphic in at least one population (1000 Genomes, 1000GENOMES, CSHL-HAPMAP, EGP_SNPS, NHLBI-ESP, PGA-UW-FHCRC). A variation is polymorphic if the minor allele frequency is at least 0.01 and the minor allele is present in at least two samples. (file common_all.vcf.gz dated 2012-06-26) (dbSNP137)",
        common_not_med => "ALT allele variations from common_all.vcf.gz that do not meet the clinical criteria. A clinical variation is one the appears in clinvar_YYYYMMDD.vcf.gz with at least one of the following clinical significance codes: 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other (file common_no_known_medical_impact_20120616.vcf.gz dated 2012-06-26) (dbSNP137)",
        clinvar => "ALT allele variations from clinvar (presumably variant is clinical (LSDB,OMIM,TPA,Diagnostic)) (file clinvar_20120616.vcf.gz dated 2012-06-16) (dbSNP137)",
        CLNSIG => "ALT allele clinical significance, 0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other (note: multiple values are possible separated by , and/or | and not obviously corresponding to a particular ALT) (dbSNP137)",
        total_allele_count => "total number of alleles observed in all called genotypes (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        ref_allele_count => "number of REF alleles in called genotypes (dbSNP137,1000Genomes,NHLBI,NIEHS)",
        ref_allele_frequency => "REF allele frequency (dbSNP137,1000Genomes,NHLBI,NIEHS)",
        alt_allele_count => "number of ALT alleles in called genotypes (dbSNP137,1000Genomes,NHLBI,NIEHS)",
        alt_allele_frequency => "ALT allele frequency (dbSNP137,1000Genomes,NHLBI,NIEHS)",
        total_sample_count => "total number of samples in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        wt_sample_count => "number of samples with REF/REF genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        wt_genotype_frequency => "frequency of REF/REF genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        het_sample_count => "number of samples with REF/ALT genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        het_genotype_frequency => "frequency of REF/ALT genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        hom_sample_count => "number of samples with ALT/ALT genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        hom_genotype_frequency => "frequency of ALT/ALT genotype in all genotyped populations (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_disruptive_wt_count => "number of wild-type genotypes in gene at loci with at least one worse-than-missense variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_disruptive_het_count => "number of REF/ALT genotypes in gene with a worse-than-missense variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_disruptive_hom_count => "number of ALT/ALT genotypes in gene with a worse-than-missense variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_disruptive_het_freq => "frequency of REF/ALT genotypes in gene with a worse-than-missense variant (gene_disruptive_het_count/(gene_disruptive_wt_count+gene_disruptive_het_count+gene_disruptive_hom_count)) (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_disruptive_hom_freq => "frequency of ALT/ALT genotypes in gene with a worse-than-missense variant (gene_disruptive_hom_count/(gene_disruptive_wt_count+gene_disruptive_het_count+gene_disruptive_hom_count)) (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_missense_deleterious_wt_count => "number of wild-type genotypes in gene at loci with at least one missense (Condel deleterious) variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_missense_deleterious_het_count => "number of REF/ALT genotypes in gene with a missense (Condel deleterious) variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_missense_deleterious_hom_count => "number of ALT/ALT genotypes in gene with a missense (Condel deleterious) variant (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_missense_deleterious_het_freq => "frequency of REF/ALT genotypes in gene with a missense (Condel deleterious) variant (gene_missense_deleterious_het_count/(gene_missense_deleterious_wt_count+gene_missense_deleterious_het_count+gene_missense_deleterious_hom_count)) (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        gene_missense_deleterious_hom_freq => "frequency of ALT/ALT genotypes in gene with a missense (Condel deleterious) variant (gene_missense_deleterious_hom_count/(gene_missense_deleterious_wt_count+gene_missense_deleterious_het_count+gene_missense_deleterious_hom_count)) (dbSNP137,NHLBI,NIEHS,1000Genomes)",
        clinvar_count => "number of clinvar variations in gene (presumably variant is clinical (LSDB,OMIM,TPA,Diagnostic)) (file clinvar_20120616.vcf.gz dated 2012-06-16) (dbSNP137)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};

    my ($chrom,$pos,$id,$ref,$alt) = split(/\t/,$input_line);

    if (defined $chrom && defined $pos && defined $ref && defined $alt){
        my $query = "CALL $vw::vw_database.coord2variant_frequencies('$chrom',$pos,'$ref','$alt')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{in_dbsnp} = defined($row[4]) ? $row[4] : '';
            $line_hash->{in_dbsnp_pop} = defined($row[5]) ? $row[5] : '';
            $line_hash->{in_nhlbi} = defined($row[6]) ? $row[6] : '';
            $line_hash->{in_niehs} = defined($row[7]) ? $row[7] : '';
            $line_hash->{in_1kg} = defined($row[8]) ? $row[8] : '';
            $line_hash->{in_1kg_pop} = defined($row[9]) ? $row[9] : '';
            $line_hash->{G5} = defined($row[10]) ? $row[10] : '';
            $line_hash->{G5A} = defined($row[11]) ? $row[11] : '';
            $line_hash->{common} = defined($row[12]) ? $row[12] : '';
            $line_hash->{common_not_med} = defined($row[13]) ? $row[13] : '';
            $line_hash->{clinvar} = defined($row[14]) ? $row[14] : '';
            $line_hash->{CLNSIG} = defined($row[15]) ? $row[15] : '';
            $line_hash->{total_allele_count} = defined($row[16]) ? $row[16] : '';
            $line_hash->{alt_allele_count} = defined($row[17]) ? $row[17] : '';
            $line_hash->{alt_allele_frequency} = defined($row[18]) ? $row[18] : '';
            $line_hash->{dbsnp137_ID} = defined($row[19]) ? $row[19] : '';
            last;
        }
        $query = "CALL $vw::vw_database.coord2variant_frequencies('$chrom',$pos,'$ref','$ref')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{ref_allele_count} = defined($row[17]) ? $row[17] : '';
            $line_hash->{ref_allele_frequency} = defined($row[18]) ? $row[18] : '';
        }
        $query = "CALL $vw::vw_database.coord2genotype_frequencies('$chrom',$pos,'$ref','$ref','$ref')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{total_sample_count} = defined($row[10]) ? $row[10] : '';
            $line_hash->{wt_sample_count} = defined($row[5]) ? $row[5] : '';
            $line_hash->{wt_genotype_frequency} = defined($row[11]) ? $row[11] : '';
            last;
        }
        $query = "CALL $vw::vw_database.coord2genotype_frequencies('$chrom',$pos,'$ref','$ref','$alt')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{het_sample_count} = defined($row[5]) ? $row[5] : '';
            $line_hash->{het_genotype_frequency} = defined($row[11]) ? $row[11] : '';
            last;
        }
        $query = "CALL $vw::vw_database.coord2genotype_frequencies('$chrom',$pos,'$ref','$alt','$alt')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{hom_sample_count} = defined($row[5]) ? $row[5] : '';
            $line_hash->{hom_genotype_frequency} = defined($row[11]) ? $row[11] : '';
            last;
        }
    }
    if (defined $bvfoa_info{gene}){
        my $query = "CALL $vw::vw_database.ensg2fucked_up_gene('$bvfoa_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{gene_disruptive_hom_count} = defined($row[2]) ? $row[2] : '';
            $line_hash->{gene_disruptive_het_count} = defined($row[3]) ? $row[3] : '';
            $line_hash->{gene_disruptive_wt_count} = defined($row[4]) ? $row[4] : '';
            $line_hash->{gene_disruptive_hom_freq} = defined($row[5]) ? $row[5] : '';
            $line_hash->{gene_disruptive_het_freq} = defined($row[6]) ? $row[6] : '';
            $line_hash->{gene_missense_deleterious_hom_count} = defined($row[7]) ? $row[7] : '';
            $line_hash->{gene_missense_deleterious_het_count} = defined($row[8]) ? $row[8] : '';
            $line_hash->{gene_missense_deleterious_wt_count} = defined($row[9]) ? $row[9] : '';
            $line_hash->{gene_missense_deleterious_hom_freq} = defined($row[10]) ? $row[10] : '';
            $line_hash->{gene_missense_deleterious_het_freq} = defined($row[11]) ? $row[11] : '';
            $line_hash->{clinvar_count} = defined($row[12]) ? $row[12] : '';
        }
    }    
    return {};
}


1;

