=head1 LICENSE

 Copyright (c) 2012 Aliz R. Rao.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Aliz R. Rao <alizrrao@gmail.com>
    
=cut

=head1 NAME

 AlleleFrequencies

=head1 SYNOPSIS

 mv EthnicAlleleFrequencies.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin EthnicAlleleFrequencies,cortex.local,3306,vw,vw,mysql,vw

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 ASN_alt_allele_frequency,ASN_sample_count,ASN_allele_count,ASN_ref_allele_count,ASN_alt_allele_count,ASN_het_sample_count,ASN_hom_sample_count,EUR_alt_allele_frequency,EUR_sample_count,EUR_allele_count,EUR_ref_allele_count,EUR_alt_allele_count,EUR_het_sample_count,EUR_hom_sample_count,AFR_alt_allele_frequency,AFR_sample_count,AFR_allele_count,AFR_ref_allele_count,AFR_alt_allele_count,AFR_het_sample_count,AFR_hom_sample_count,AMR_alt_allele_frequency,AMR_sample_count,AMR_allele_count,AMR_ref_allele_count,AMR_alt_allele_count,AMR_het_sample_count,AMR_hom_sample_count,SAN_alt_allele_frequency,SAN_sample_count,SAN_allele_count,SAN_ref_allele_count,SAN_alt_allele_count,SAN_het_sample_count,SAN_hom_sample_count 

=head1 PARAMETERS

=cut

package EthnicAlleleFrequencies;

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
        ASW_alt_allele_frequency
        ASW_sample_count
        ASW_allele_count
        ASW_ref_allele_count
        ASW_alt_allele_count
        ASW_het_sample_count
        ASW_hom_sample_count
        CEU_alt_allele_frequency
        CEU_sample_count
        CEU_allele_count
        CEU_ref_allele_count
        CEU_alt_allele_count
        CEU_het_sample_count
        CEU_hom_sample_count
        CHB_alt_allele_frequency
        CHB_sample_count
        CHB_allele_count
        CHB_ref_allele_count
        CHB_alt_allele_count
        CHB_het_sample_count
        CHB_hom_sample_count
        CHD_alt_allele_frequency
        CHD_sample_count
        CHD_allele_count
        CHD_ref_allele_count
        CHD_alt_allele_count
        CHD_het_sample_count
        CHD_hom_sample_count
        GIH_alt_allele_frequency
        GIH_sample_count
        GIH_allele_count
        GIH_ref_allele_count
        GIH_alt_allele_count
        GIH_het_sample_count
        GIH_hom_sample_count
        HCB_alt_allele_frequency
        HCB_sample_count
        HCB_allele_count
        HCB_ref_allele_count
        HCB_alt_allele_count
        HCB_het_sample_count
        HCB_hom_sample_count
        JPT_alt_allele_frequency
        JPT_sample_count
        JPT_allele_count
        JPT_ref_allele_count
        JPT_alt_allele_count
        JPT_het_sample_count
        JPT_hom_sample_count
        LWK_alt_allele_frequency
        LWK_sample_count
        LWK_allele_count
        LWK_ref_allele_count
        LWK_alt_allele_count
        LWK_het_sample_count
        LWK_hom_sample_count
        MEX_alt_allele_frequency
        MEX_sample_count
        MEX_allele_count
        MEX_ref_allele_count
        MEX_alt_allele_count
        MEX_het_sample_count
        MEX_hom_sample_count
        MKK_alt_allele_frequency
        MKK_sample_count
        MKK_allele_count
        MKK_ref_allele_count
        MKK_alt_allele_count
        MKK_het_sample_count
        MKK_hom_sample_count
        TSI_alt_allele_frequency
        TSI_sample_count
        TSI_allele_count
        TSI_ref_allele_count
        TSI_alt_allele_count
        TSI_het_sample_count
        TSI_hom_sample_count
        YRI_alt_allele_frequency
        YRI_sample_count
        YRI_allele_count
        YRI_ref_allele_count
        YRI_alt_allele_count
        YRI_het_sample_count
        YRI_hom_sample_count
        min_alt_allele_frequency_dbsnp
        max_alt_allele_frequency_dbsnp
        AC_1kg
        AF_1kg
        AFR_AF_1kg
        AMR_AF_1kg
        AN_1kg
        ASN_AF_1kg
        EUR_AF_1kg
        AA_AC_ALT_nhlbi
        AA_AC_REF_nhlbi
        EA_AC_ALT_nhlbi
        EA_AC_REF_nhlbi
        All_AC_ALT_nhlbi
        All_AC_REF_nhlbi
        AA_GTC_hom_alt_nhlbi
        AA_GTC_het_nhlbi
        AA_GTC_hom_ref_nhlbi
        AA_GTC_di_alt_nhlbi
        EA_GTC_hom_alt_nhlbi
        EA_GTC_het_nhlbi
        EA_GTC_hom_ref_nhlbi
        EA_GTC_di_alt_nhlbi
        GTC_hom_alt_nhlbi
        GTC_het_nhlbi
        GTC_hom_ref_nhlbi
        GTC_di_alt_nhlbi
        AA_ALT_AF_nhlbi
        EA_ALT_AF_nhlbi
        All_ALT_AF_nhlbi
        AA_REF_AF_nhlbi
        EA_REF_AF_nhlbi
        All_REF_AF_nhlbi
        AA_MAF_nhlbi
        EA_MAF_nhlbi
        All_MAF_nhlbi
        AC_niehs
        AF_niehs
        AN_niehs
        ASN_alt_allele_frequency
        ASN_sample_count
        ASN_allele_count
        ASN_ref_allele_count
        ASN_alt_allele_count
        ASN_het_sample_count
        ASN_hom_sample_count
        EUR_alt_allele_frequency
        EUR_sample_count
        EUR_allele_count
        EUR_ref_allele_count
        EUR_alt_allele_count
        EUR_het_sample_count
        EUR_hom_sample_count
        AFR_alt_allele_frequency
        AFR_sample_count
        AFR_allele_count
        AFR_ref_allele_count
        AFR_alt_allele_count
        AFR_het_sample_count
        AFR_hom_sample_count
        AMR_alt_allele_frequency
        AMR_sample_count
        AMR_allele_count
        AMR_ref_allele_count
        AMR_alt_allele_count
        AMR_het_sample_count
        AMR_hom_sample_count
        SAN_alt_allele_frequency
        SAN_sample_count
        SAN_allele_count
        SAN_ref_allele_count
        SAN_alt_allele_count
        SAN_het_sample_count
        SAN_hom_sample_count
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        ASW_alt_allele_frequency => "frequency of alt allele in ASW dbsnp135 population",
        ASW_sample_count => "number of samples in ASW dbsnp135 population",
        ASW_allele_count => "number of alleles in called genotypes in ASW dbsnp135 population",
        ASW_ref_allele_count => "number of reference alleles in called genotypes in ASW dbsnp135 population",
        ASW_alt_allele_count => "number of alternate alleles in called genotypes in ASW dbsnp135 population",
        ASW_het_sample_count => "number of heterozygous samples in ASW dbsnp135 population",
        ASW_hom_sample_count => "number of homozygous alternate samples in ASW dbsnp135 population",
        CEU_alt_allele_frequency => "frequency of alt allele in CEU dbsnp135 population",
        CEU_sample_count => "number of samples in CEU dbsnp135 population",
        CEU_allele_count => "number of alleles in called genotypes in CEU dbsnp135 population",
        CEU_ref_allele_count => "number of reference alleles in called genotypes in CEU dbsnp135 population",
        CEU_alt_allele_count => "number of alternate alleles in called genotypes in CEU dbsnp135 population",
        CEU_het_sample_count => "number of heterozygous samples in CEU dbsnp135 population",
        CEU_hom_sample_count => "number of homozygous alternate samples in CEU dbsnp135 population",
        CHB_alt_allele_frequency => "frequency of alt allele in CHB dbsnp135 population",
        CHB_sample_count => "number of samples in CHB dbsnp135 population",
        CHB_allele_count => "number of alleles in called genotypes in CHB dbsnp135 population",
        CHB_ref_allele_count => "number of reference alleles in called genotypes in CHB dbsnp135 population",
        CHB_alt_allele_count => "number of alternate alleles in called genotypes in CHB dbsnp135 population",
        CHB_het_sample_count => "number of heterozygous samples in CHB dbsnp135 population",
        CHB_hom_sample_count => "number of homozygous alternate samples in CHB dbsnp135 population",
        CHD_alt_allele_frequency => "frequency of alt allele in CHD dbsnp135 population",
        CHD_sample_count => "number of samples in CHD dbsnp135 population",
        CHD_allele_count => "number of alleles in called genotypes in CHD dbsnp135 population",
        CHD_ref_allele_count => "number of reference alleles in called genotypes in CHD dbsnp135 population",
        CHD_alt_allele_count => "number of alternate alleles in called genotypes in CHD dbsnp135 population",
        CHD_het_sample_count => "number of heterozygous samples in CHD dbsnp135 population",
        CHD_hom_sample_count => "number of homozygous alternate samples in CHD dbsnp135 population",
        GIH_alt_allele_frequency => "frequency of alt allele in GIH dbsnp135 population",
        GIH_sample_count => "number of samples in GIH dbsnp135 population",
        GIH_allele_count => "number of alleles in called genotypes in GIH dbsnp135 population",
        GIH_ref_allele_count => "number of reference alleles in called genotypes in GIH dbsnp135 population",
        GIH_alt_allele_count => "number of alternate alleles in called genotypes in GIH dbsnp135 population",
        GIH_het_sample_count => "number of heterozygous samples in GIH dbsnp135 population",
        GIH_hom_sample_count => "number of homozygous alternate samples in GIH dbsnp135 population",
        HCB_alt_allele_frequency => "frequency of alt allele in HCB dbsnp135 population",
        HCB_sample_count => "number of samples in HCB dbsnp135 population",
        HCB_allele_count => "number of alleles in called genotypes in HCB dbsnp135 population",
        HCB_ref_allele_count => "number of reference alleles in called genotypes in HCB dbsnp135 population",
        HCB_alt_allele_count => "number of alternate alleles in called genotypes in HCB dbsnp135 population",
        HCB_het_sample_count => "number of heterozygous samples in HCB dbsnp135 population",
        HCB_hom_sample_count => "number of homozygous alternate samples in HCB dbsnp135 population",
        JPT_alt_allele_frequency => "frequency of alt allele in JPT dbsnp135 population",
        JPT_sample_count => "number of samples in JPT dbsnp135 population",
        JPT_allele_count => "number of alleles in called genotypes in JPT dbsnp135 population",
        JPT_ref_allele_count => "number of reference alleles in called genotypes in JPT dbsnp135 population",
        JPT_alt_allele_count => "number of alternate alleles in called genotypes in JPT dbsnp135 population",
        JPT_het_sample_count => "number of heterozygous samples in JPT dbsnp135 population",
        JPT_hom_sample_count => "number of homozygous alternate samples in JPT dbsnp135 population",
        LWK_alt_allele_frequency => "frequency of alt allele in LWK dbsnp135 population",
        LWK_sample_count => "number of samples in LWK dbsnp135 population",
        LWK_allele_count => "number of alleles in called genotypes in LWK dbsnp135 population",
        LWK_ref_allele_count => "number of reference alleles in called genotypes in LWK dbsnp135 population",
        LWK_alt_allele_count => "number of alternate alleles in called genotypes in LWK dbsnp135 population",
        LWK_het_sample_count => "number of heterozygous samples in LWK dbsnp135 population",
        LWK_hom_sample_count => "number of homozygous alternate samples in LWK dbsnp135 population",
        MEX_alt_allele_frequency => "frequency of alt allele in MEX dbsnp135 population",
        MEX_sample_count => "number of samples in MEX dbsnp135 population",
        MEX_allele_count => "number of alleles in called genotypes in MEX dbsnp135 population",
        MEX_ref_allele_count => "number of reference alleles in called genotypes in MEX dbsnp135 population",
        MEX_alt_allele_count => "number of alternate alleles in called genotypes in MEX dbsnp135 population",
        MEX_het_sample_count => "number of heterozygous samples in MEX dbsnp135 population",
        MEX_hom_sample_count => "number of homozygous alternate samples in MEX dbsnp135 population",
        MKK_alt_allele_frequency => "frequency of alt allele in MKK dbsnp135 population",
        MKK_sample_count => "number of samples in MKK dbsnp135 population",
        MKK_allele_count => "number of alleles in called genotypes in MKK dbsnp135 population",
        MKK_ref_allele_count => "number of reference alleles in called genotypes in MKK dbsnp135 population",
        MKK_alt_allele_count => "number of alternate alleles in called genotypes in MKK dbsnp135 population",
        MKK_het_sample_count => "number of heterozygous samples in MKK dbsnp135 population",
        MKK_hom_sample_count => "number of homozygous alternate samples in MKK dbsnp135 population",
        TSI_alt_allele_frequency => "frequency of alt allele in TSI dbsnp135 population",
        TSI_sample_count => "number of samples in TSI dbsnp135 population",
        TSI_allele_count => "number of alleles in called genotypes in TSI dbsnp135 population",
        TSI_ref_allele_count => "number of reference alleles in called genotypes in TSI dbsnp135 population",
        TSI_alt_allele_count => "number of alternate alleles in called genotypes in TSI dbsnp135 population",
        TSI_het_sample_count => "number of heterozygous samples in TSI dbsnp135 population",
        TSI_hom_sample_count => "number of homozygous alternate samples in TSI dbsnp135 population",
        YRI_alt_allele_frequency => "frequency of alt allele in YRI dbsnp135 population",
        YRI_sample_count => "number of samples in YRI dbsnp135 population",
        YRI_allele_count => "number of alleles in called genotypes in YRI dbsnp135 population",
        YRI_ref_allele_count => "number of reference alleles in called genotypes in YRI dbsnp135 population",
        YRI_alt_allele_count => "number of alternate alleles in called genotypes in YRI dbsnp135 population",
        YRI_het_sample_count => "number of heterozygous samples in YRI dbsnp135 population",
        YRI_hom_sample_count => "number of homozygous alternate samples in YRI dbsnp135 population",
        min_alt_allele_frequency_dbsnp => "minimum alternate allele frequency in dbSNP135 populationsmin population",
        max_alt_allele_frequency_dbsnp => "maximum alternate allele frequency in dbSNP135 populations",
        AC_1kg => "Alternate Allele Count in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        AF_1kg => "Global Allele Frequency based on AC/AN in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        AFR_AF_1kg => "Allele Frequency for samples from AFR based on AC/AN in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        AMR_AF_1kg => "Allele Frequency for samples from AMR based on AC/AN in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        AN_1kg => "Total Allele Count in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        ASN_AF_1kg => "Allele Frequency for samples from ASN based on AC/AN in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        EUR_AF_1kg => "Allele Frequency for samples from EUR based on AC/AN in 1000 Genomes ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf",
        AA_AC_ALT_nhlbi => "African American ALT Allele Count in NHLBI ESPS5400.snps",
        AA_AC_REF_nhlbi => "African American REF Allele Count in NHLBI ESPS5400.snps",
        EA_AC_ALT_nhlbi => "European American ALT Allele Count in NHLBI ESPS5400.snps",
        EA_AC_REF_nhlbi => "European American REF Allele Count in NHLBI ESPS5400.snps",
        All_AC_ALT_nhlbi => "Total ALT Allele Count in NHLBI ESPS5400.snps",
        All_AC_REF_nhlbi => "Total REF Allele Count in NHLBI ESPS5400.snps",
        AA_GTC_hom_alt_nhlbi => "African American homozygous ALT Genotype Count in NHLBI ESPS5400.snps",
        AA_GTC_het_nhlbi => "African American heterozygous REF/ALT Genotype Count in NHLBI ESPS5400.snps",
        AA_GTC_hom_ref_nhlbi => "African American homozygous REFGenotype Count in NHLBI ESPS5400.snps",
        AA_GTC_di_alt_nhlbi => "African American heterozygous ALT/ALTx Genotype Count in NHLBI ESPS5400.snps",
        EA_GTC_hom_alt_nhlbi => "European American homozygous ALT Genotype Count in NHLBI ESPS5400.snps",
        EA_GTC_het_nhlbi => "European Americanheterozygous REF/ALT Genotype Count in NHLBI ESPS5400.snps",
        EA_GTC_hom_ref_nhlbi => "European American homozygous REFGenotype Count in NHLBI ESPS5400.snps",
        EA_GTC_di_alt_nhlbi => "European American heterozygous ALT/ALTx Genotype Count in NHLBI ESPS5400.snps",
        GTC_hom_alt_nhlbi => "Total homozygous ALT Genotype Count in NHLBI ESPS5400.snps",
        GTC_het_nhlbi => "Total heterozygous REF/ALT Genotype Count in NHLBI ESPS5400.snps",
        GTC_hom_ref_nhlbi => "Total homozygous REFGenotype Count in NHLBI ESPS5400.snps",
        GTC_di_alt_nhlbi => "Total heterozygous ALT/ALTx Genotype Count in NHLBI ESPS5400.snps",
        AA_ALT_AF_nhlbi => "African American alternate allele frequency in NHLBI ESPS5400.snps",
        EA_ALT_AF_nhlbi => "European American alternate allele frequency in NHLBI ESPS5400.snps",
        All_ALT_AF_nhlbi => "Total alternate allele frequency in NHLBI ESPS5400.snps",
        AA_REF_AF_nhlbi => "African American reference allele frequency in NHLBI ESPS5400.snps",
        EA_REF_AF_nhlbi => "European American reference allele frequency in NHLBI ESPS5400.snps",
        All_REF_AF_nhlbi => "Total reference allele frequency in NHLBI ESPS5400.snps",
        AA_MAF_nhlbi => "African American minor allele frequency in NHLBI ESPS5400.snps",
        EA_MAF_nhlbi => "European American minor allele frequency in NHLBI ESPS5400.snps",
        All_MAF_nhlbi => "Total minor allele frequency in NHLBI ESPS5400.snps",
        AC_niehs => "Allele count in genotypesin NIEHS niehs.95samples.polymorphic.filtered.{snps,indels}.vcf",
        AF_niehs => "Allele Frequency for ALT allelein NIEHS niehs.95samples.polymorphic.filtered.{snps,indels}.vcf",
        AN_niehs => "Total number of alleles in called genotypesin NIEHS niehs.95samples.polymorphic.filtered.{snps,indels}.vcf",
        ASN_alt_allele_frequency => "frequency of alt allele in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_sample_count => "number of samples in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_allele_count => "number of alleles in called genotypes in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_ref_allele_count => "number of reference alleles in called genotypes in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_alt_allele_count => "number of alternate alleles in called genotypes in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_het_sample_count => "number of heterozygous samples in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        ASN_hom_sample_count => "number of homozygous alternate samples in East Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_alt_allele_frequency => "frequency of alt allele in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_sample_count => "number of samples in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_allele_count => "number of alleles in called genotypes in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_ref_allele_count => "number of reference alleles in called genotypes in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_alt_allele_count => "number of alternate alleles in called genotypes in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_het_sample_count => "number of heterozygous samples in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        EUR_hom_sample_count => "number of homozygous alternate samples in European populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_alt_allele_frequency => "frequency of alt allele in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_sample_count => "number of samples in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_allele_count => "number of alleles in called genotypes in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_ref_allele_count => "number of reference alleles in called genotypes in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_alt_allele_count => "number of alternate alleles in called genotypes in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_het_sample_count => "number of heterozygous samples in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AFR_hom_sample_count => "number of homozygous alternate samples in African populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_alt_allele_frequency => "frequency of alt allele in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_sample_count => "number of samples in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_allele_count => "number of alleles in called genotypes in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_ref_allele_count => "number of reference alleles in called genotypes in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_alt_allele_count => "number of alternate alleles in called genotypes in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_het_sample_count => "number of heterozygous samples in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        AMR_hom_sample_count => "number of homozygous alternate samples in Ad Mixed American populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_alt_allele_frequency => "frequency of alt allele in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_sample_count => "number of samples in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_allele_count => "number of alleles in called genotypes in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_ref_allele_count => "number of reference alleles in called genotypes in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_alt_allele_count => "number of alternate alleles in called genotypes in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_het_sample_count => "number of heterozygous samples in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
        SAN_hom_sample_count => "number of homozygous alternate samples in South Asian populations in 1000 Genomes, dbSNP135 and NHLBI",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};

    my ($chrom,$pos,$id,$ref,$alt) = split(/\t/,$input_line);

    if (defined $chrom && defined $pos && defined $ref && defined $alt){
        my $query = "CALL $vw::vw_database.coord2ethnicaf('$chrom',$pos,'$ref','$alt')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{ASW_alt_allele_frequency} = defined($row[0]) ? $row[0] : '';
            $line_hash->{ASW_sample_count} = defined($row[1]) ? $row[1] : '';
            $line_hash->{ASW_allele_count} = defined($row[2]) ? $row[2] : '';
            $line_hash->{ASW_ref_allele_count} = defined($row[3]) ? $row[3] : '';
            $line_hash->{ASW_alt_allele_count} = defined($row[4]) ? $row[4] : '';
            $line_hash->{ASW_het_sample_count} = defined($row[5]) ? $row[5] : '';
            $line_hash->{ASW_hom_sample_count} = defined($row[6]) ? $row[6] : '';
            $line_hash->{CEU_alt_allele_frequency} = defined($row[7]) ? $row[7] : '';
            $line_hash->{CEU_sample_count} = defined($row[8]) ? $row[8] : '';
            $line_hash->{CEU_allele_count} = defined($row[9]) ? $row[9] : '';
            $line_hash->{CEU_ref_allele_count} = defined($row[10]) ? $row[10] : '';
            $line_hash->{CEU_alt_allele_count} = defined($row[11]) ? $row[11] : '';
            $line_hash->{CEU_het_sample_count} = defined($row[12]) ? $row[12] : '';
            $line_hash->{CEU_hom_sample_count} = defined($row[13]) ? $row[13] : '';
            $line_hash->{CHB_alt_allele_frequency} = defined($row[14]) ? $row[14] : '';
            $line_hash->{CHB_sample_count} = defined($row[15]) ? $row[15] : '';
            $line_hash->{CHB_allele_count} = defined($row[16]) ? $row[16] : '';
            $line_hash->{CHB_ref_allele_count} = defined($row[17]) ? $row[17] : '';
            $line_hash->{CHB_alt_allele_count} = defined($row[18]) ? $row[18] : '';
            $line_hash->{CHB_het_sample_count} = defined($row[19]) ? $row[19] : '';
            $line_hash->{CHB_hom_sample_count} = defined($row[20]) ? $row[20] : '';
            $line_hash->{CHD_alt_allele_frequency} = defined($row[21]) ? $row[21] : '';
            $line_hash->{CHD_sample_count} = defined($row[22]) ? $row[22] : '';
            $line_hash->{CHD_allele_count} = defined($row[23]) ? $row[23] : '';
            $line_hash->{CHD_ref_allele_count} = defined($row[24]) ? $row[24] : '';
            $line_hash->{CHD_alt_allele_count} = defined($row[25]) ? $row[25] : '';
            $line_hash->{CHD_het_sample_count} = defined($row[26]) ? $row[26] : '';
            $line_hash->{CHD_hom_sample_count} = defined($row[27]) ? $row[27] : '';
            $line_hash->{GIH_alt_allele_frequency} = defined($row[28]) ? $row[28] : '';
            $line_hash->{GIH_sample_count} = defined($row[29]) ? $row[29] : '';
            $line_hash->{GIH_allele_count} = defined($row[30]) ? $row[30] : '';
            $line_hash->{GIH_ref_allele_count} = defined($row[31]) ? $row[31] : '';
            $line_hash->{GIH_alt_allele_count} = defined($row[32]) ? $row[32] : '';
            $line_hash->{GIH_het_sample_count} = defined($row[33]) ? $row[33] : '';
            $line_hash->{GIH_hom_sample_count} = defined($row[34]) ? $row[34] : '';
            $line_hash->{HCB_alt_allele_frequency} = defined($row[35]) ? $row[35] : '';
            $line_hash->{HCB_sample_count} = defined($row[36]) ? $row[36] : '';
            $line_hash->{HCB_allele_count} = defined($row[37]) ? $row[37] : '';
            $line_hash->{HCB_ref_allele_count} = defined($row[38]) ? $row[38] : '';
            $line_hash->{HCB_alt_allele_count} = defined($row[39]) ? $row[39] : '';
            $line_hash->{HCB_het_sample_count} = defined($row[40]) ? $row[40] : '';
            $line_hash->{HCB_hom_sample_count} = defined($row[41]) ? $row[41] : '';
            $line_hash->{JPT_alt_allele_frequency} = defined($row[42]) ? $row[42] : '';
            $line_hash->{JPT_sample_count} = defined($row[43]) ? $row[43] : '';
            $line_hash->{JPT_allele_count} = defined($row[44]) ? $row[44] : '';
            $line_hash->{JPT_ref_allele_count} = defined($row[45]) ? $row[45] : '';
            $line_hash->{JPT_alt_allele_count} = defined($row[46]) ? $row[46] : '';
            $line_hash->{JPT_het_sample_count} = defined($row[47]) ? $row[47] : '';
            $line_hash->{JPT_hom_sample_count} = defined($row[48]) ? $row[48] : '';
            $line_hash->{LWK_alt_allele_frequency} = defined($row[49]) ? $row[49] : '';
            $line_hash->{LWK_sample_count} = defined($row[50]) ? $row[50] : '';
            $line_hash->{LWK_allele_count} = defined($row[51]) ? $row[51] : '';
            $line_hash->{LWK_ref_allele_count} = defined($row[52]) ? $row[52] : '';
            $line_hash->{LWK_alt_allele_count} = defined($row[53]) ? $row[53] : '';
            $line_hash->{LWK_het_sample_count} = defined($row[54]) ? $row[54] : '';
            $line_hash->{LWK_hom_sample_count} = defined($row[55]) ? $row[55] : '';
            $line_hash->{MEX_alt_allele_frequency} = defined($row[56]) ? $row[56] : '';
            $line_hash->{MEX_sample_count} = defined($row[57]) ? $row[57] : '';
            $line_hash->{MEX_allele_count} = defined($row[58]) ? $row[58] : '';
            $line_hash->{MEX_ref_allele_count} = defined($row[59]) ? $row[59] : '';
            $line_hash->{MEX_alt_allele_count} = defined($row[60]) ? $row[60] : '';
            $line_hash->{MEX_het_sample_count} = defined($row[61]) ? $row[61] : '';
            $line_hash->{MEX_hom_sample_count} = defined($row[62]) ? $row[62] : '';
            $line_hash->{MKK_alt_allele_frequency} = defined($row[63]) ? $row[63] : '';
            $line_hash->{MKK_sample_count} = defined($row[64]) ? $row[64] : '';
            $line_hash->{MKK_allele_count} = defined($row[65]) ? $row[65] : '';
            $line_hash->{MKK_ref_allele_count} = defined($row[66]) ? $row[66] : '';
            $line_hash->{MKK_alt_allele_count} = defined($row[67]) ? $row[67] : '';
            $line_hash->{MKK_het_sample_count} = defined($row[68]) ? $row[68] : '';
            $line_hash->{MKK_hom_sample_count} = defined($row[69]) ? $row[69] : '';
            $line_hash->{TSI_alt_allele_frequency} = defined($row[70]) ? $row[70] : '';
            $line_hash->{TSI_sample_count} = defined($row[71]) ? $row[71] : '';
            $line_hash->{TSI_allele_count} = defined($row[72]) ? $row[72] : '';
            $line_hash->{TSI_ref_allele_count} = defined($row[73]) ? $row[73] : '';
            $line_hash->{TSI_alt_allele_count} = defined($row[74]) ? $row[74] : '';
            $line_hash->{TSI_het_sample_count} = defined($row[75]) ? $row[75] : '';
            $line_hash->{TSI_hom_sample_count} = defined($row[76]) ? $row[76] : '';
            $line_hash->{YRI_alt_allele_frequency} = defined($row[77]) ? $row[77] : '';
            $line_hash->{YRI_sample_count} = defined($row[78]) ? $row[78] : '';
            $line_hash->{YRI_allele_count} = defined($row[79]) ? $row[79] : '';
            $line_hash->{YRI_ref_allele_count} = defined($row[80]) ? $row[80] : '';
            $line_hash->{YRI_alt_allele_count} = defined($row[81]) ? $row[81] : '';
            $line_hash->{YRI_het_sample_count} = defined($row[82]) ? $row[82] : '';
            $line_hash->{YRI_hom_sample_count} = defined($row[83]) ? $row[83] : '';
            $line_hash->{min_alt_allele_frequency_dbsnp} = defined($row[84]) ? $row[84] : '';
            $line_hash->{max_alt_allele_frequency_dbsnp} = defined($row[85]) ? $row[85] : '';
            $line_hash->{AC_1kg} = defined($row[86]) ? $row[86] : '';
            $line_hash->{AF_1kg} = defined($row[87]) ? $row[87] : '';
            $line_hash->{AFR_AF_1kg} = defined($row[88]) ? $row[88] : '';
            $line_hash->{AMR_AF_1kg} = defined($row[89]) ? $row[89] : '';
            $line_hash->{AN_1kg} = defined($row[90]) ? $row[90] : '';
            $line_hash->{ASN_AF_1kg} = defined($row[91]) ? $row[91] : '';
            $line_hash->{EUR_AF_1kg} = defined($row[92]) ? $row[92] : '';
            $line_hash->{AA_AC_ALT_nhlbi} = defined($row[93]) ? $row[93] : '';
            $line_hash->{AA_AC_REF_nhlbi} = defined($row[94]) ? $row[94] : '';
            $line_hash->{EA_AC_ALT_nhlbi} = defined($row[95]) ? $row[95] : '';
            $line_hash->{EA_AC_REF_nhlbi} = defined($row[96]) ? $row[96] : '';
            $line_hash->{All_AC_ALT_nhlbi} = defined($row[97]) ? $row[97] : '';
            $line_hash->{All_AC_REF_nhlbi} = defined($row[98]) ? $row[98] : '';
            $line_hash->{AA_GTC_hom_alt_nhlbi} = defined($row[99]) ? $row[99] : '';
            $line_hash->{AA_GTC_het_nhlbi} = defined($row[100]) ? $row[100] : '';
            $line_hash->{AA_GTC_hom_ref_nhlbi} = defined($row[101]) ? $row[101] : '';
            $line_hash->{AA_GTC_di_alt_nhlbi} = defined($row[102]) ? $row[102] : '';
            $line_hash->{EA_GTC_hom_alt_nhlbi} = defined($row[103]) ? $row[103] : '';
            $line_hash->{EA_GTC_het_nhlbi} = defined($row[104]) ? $row[104] : '';
            $line_hash->{EA_GTC_hom_ref_nhlbi} = defined($row[105]) ? $row[105] : '';
            $line_hash->{EA_GTC_di_alt_nhlbi} = defined($row[106]) ? $row[106] : '';
            $line_hash->{GTC_hom_alt_nhlbi} = defined($row[107]) ? $row[107] : '';
            $line_hash->{GTC_het_nhlbi} = defined($row[108]) ? $row[108] : '';
            $line_hash->{GTC_hom_ref_nhlbi} = defined($row[109]) ? $row[109] : '';
            $line_hash->{GTC_di_alt_nhlbi} = defined($row[110]) ? $row[110] : '';
            $line_hash->{AA_ALT_AF_nhlbi} = defined($row[111]) ? $row[111] : '';
            $line_hash->{EA_ALT_AF_nhlbi} = defined($row[112]) ? $row[112] : '';
            $line_hash->{All_ALT_AF_nhlbi} = defined($row[113]) ? $row[113] : '';
            $line_hash->{AA_REF_AF_nhlbi} = defined($row[114]) ? $row[114] : '';
            $line_hash->{EA_REF_AF_nhlbi} = defined($row[115]) ? $row[115] : '';
            $line_hash->{All_REF_AF_nhlbi} = defined($row[116]) ? $row[116] : '';
            $line_hash->{AA_MAF_nhlbi} = defined($row[117]) ? $row[117] : '';
            $line_hash->{EA_MAF_nhlbi} = defined($row[118]) ? $row[118] : '';
            $line_hash->{All_MAF_nhlbi} = defined($row[119]) ? $row[119] : '';
            $line_hash->{AC_niehs} = defined($row[120]) ? $row[120] : '';
            $line_hash->{AF_niehs} = defined($row[121]) ? $row[121] : '';
            $line_hash->{AN_niehs} = defined($row[122]) ? $row[122] : '';
            $line_hash->{ASN_alt_allele_frequency} = defined($row[123]) ? $row[123] : '';
            $line_hash->{ASN_sample_count} = defined($row[124]) ? $row[124] : '';
            $line_hash->{ASN_allele_count} = defined($row[125]) ? $row[125] : '';
            $line_hash->{ASN_ref_allele_count} = defined($row[126]) ? $row[126] : '';
            $line_hash->{ASN_alt_allele_count} = defined($row[127]) ? $row[127] : '';
            $line_hash->{ASN_het_sample_count} = defined($row[128]) ? $row[128] : '';
            $line_hash->{ASN_hom_sample_count} = defined($row[129]) ? $row[129] : '';
            $line_hash->{EUR_alt_allele_frequency} = defined($row[130]) ? $row[130] : '';
            $line_hash->{EUR_sample_count} = defined($row[131]) ? $row[131] : '';
            $line_hash->{EUR_allele_count} = defined($row[132]) ? $row[132] : '';
            $line_hash->{EUR_ref_allele_count} = defined($row[133]) ? $row[133] : '';
            $line_hash->{EUR_alt_allele_count} = defined($row[134]) ? $row[134] : '';
            $line_hash->{EUR_het_sample_count} = defined($row[135]) ? $row[135] : '';
            $line_hash->{EUR_hom_sample_count} = defined($row[136]) ? $row[136] : '';
            $line_hash->{AFR_alt_allele_frequency} = defined($row[137]) ? $row[137] : '';
            $line_hash->{AFR_sample_count} = defined($row[138]) ? $row[138] : '';
            $line_hash->{AFR_allele_count} = defined($row[139]) ? $row[139] : '';
            $line_hash->{AFR_ref_allele_count} = defined($row[140]) ? $row[140] : '';
            $line_hash->{AFR_alt_allele_count} = defined($row[141]) ? $row[141] : '';
            $line_hash->{AFR_het_sample_count} = defined($row[142]) ? $row[142] : '';
            $line_hash->{AFR_hom_sample_count} = defined($row[143]) ? $row[143] : '';
            $line_hash->{AMR_alt_allele_frequency} = defined($row[144]) ? $row[144] : '';
            $line_hash->{AMR_sample_count} = defined($row[145]) ? $row[145] : '';
            $line_hash->{AMR_allele_count} = defined($row[146]) ? $row[146] : '';
            $line_hash->{AMR_ref_allele_count} = defined($row[147]) ? $row[147] : '';
            $line_hash->{AMR_alt_allele_count} = defined($row[148]) ? $row[148] : '';
            $line_hash->{AMR_het_sample_count} = defined($row[149]) ? $row[149] : '';
            $line_hash->{AMR_hom_sample_count} = defined($row[150]) ? $row[150] : '';
            $line_hash->{SAN_alt_allele_frequency} = defined($row[151]) ? $row[151] : '';
            $line_hash->{SAN_sample_count} = defined($row[152]) ? $row[152] : '';
            $line_hash->{SAN_allele_count} = defined($row[153]) ? $row[153] : '';
            $line_hash->{SAN_ref_allele_count} = defined($row[154]) ? $row[154] : '';
            $line_hash->{SAN_alt_allele_count} = defined($row[155]) ? $row[155] : '';
            $line_hash->{SAN_het_sample_count} = defined($row[156]) ? $row[156] : '';
            $line_hash->{SAN_hom_sample_count} = defined($row[157]) ? $row[157] : '';
            last;
        }
    }
    return {};
}


1;

