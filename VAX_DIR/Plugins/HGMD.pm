=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 HGMD

=head1 SYNOPSIS

 mv HGMD.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin HGMD[,host,port,user,password,mysql,database]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 HGMD_gene_diseases, HGMD_acc_num, HGMD_locus_disease, HGMD_tag, HGMD_base, HGMD_hgvs, HGMD_codon, HGMD_amino, HGMD_deletion, HGMD_insertion, HGMD_descr.
 
 The plugin requires the commercial HGMD Pro database and the following two additional stored procedures:
  
 Requires that the VAX.pm module be in the Plugins directory

USE `hgmd_pro`;
DROP procedure IF EXISTS `coord2hgmd`;

DELIMITER $$
USE `hgmd_pro`$$
CREATE DEFINER=`hgmd`@`%` PROCEDURE `coord2hgmd`(chromosome VARCHAR(2), coordSTART INT(11), coordEND INT(11))
BEGIN
select distinct
`allmut`.`acc_num`,
`allmut`.`disease`,
`allmut`.`tag`,
`allmut`.`base`,
`allmut`.`hgvs`,
`allmut`.`codon`,
`allmut`.`amino`,
`allmut`.`deletion`,
`allmut`.`insertion`,
`allmut`.`descr`
FROM `hgmd_pro`.`hg19_coords`
JOIN `hgmd_pro`.`allmut`
ON `hg19_coords`.`acc_num` = `allmut`.`acc_num`
where `hg19_coords`.`chromosome` = chromosome
and (`hg19_coords`.`coordSTART` BETWEEN  coordSTART and coordEND
or `hg19_coords`.`coordEND` BETWEEN  coordSTART and coordEND);
END$$

DELIMITER ;


USE `hgmd_pro`;
DROP procedure IF EXISTS `gene2hgmd_disease`;

DELIMITER $$
USE `hgmd_pro`$$
CREATE DEFINER=`hgmd`@`%` PROCEDURE `gene2hgmd_disease`(gene VARCHAR(10))
BEGIN
SELECT distinct disease
from `hgmd_pro`.`allgenes`
where `allgenes`.`gene` = gene;
END$$

DELIMITER ;

 
=head1 PARAMETERS

    host (default: cortex.local),port (default: 3306),user (default: hgmd),password (default: hgmd),platform (default: mysql),database (default: hgmd_pro)

=cut

package HGMD;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_bvfoa_info get_unique);

my $hgmd_version = '2013.1';

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my ($host,$port,$user,$password,$platform,$database) = @{$self->{params}};
    $host ||= 'cortex.local';
    $port ||= 3306;
    $user ||= 'hgmd';
    $password ||= 'hgmd';
    $platform ||= 'mysql';
    $database ||= 'hgmd_pro';
    $self->{params} = [$host,$port,$user,$password,$platform,$database];
    
    my $dsn = "dbi:$platform:$database:$host:$port";
    my $conn = DBI->connect($dsn, $user, $password)
        or die "Unable to connect: $DBI::errstr\n";
    $self->{conn} = $conn;

    return $self;
}

sub version {
    return '74';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic'];
}

sub get_header_info {
    my @new_output_cols = qw(
        HGMD_gene_diseases
        HGMD_acc_num
        HGMD_locus_disease
        HGMD_tag
        HGMD_base
        HGMD_hgvs
        HGMD_codon
        HGMD_amino
        HGMD_deletion
        HGMD_insertion
        HGMD_descr
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        HGMD_gene_diseases => "Composite list of all the phenotypes for mutations in a given gene. (| separated) (HGMD $hgmd_version)",
        HGMD_acc_num => "The HGMD mutation accession number for the mutation. (HGMD $hgmd_version)",
        HGMD_locus_disease => "The name for the disease or condition associated with the mutation. (HGMD $hgmd_version)",
        HGMD_tag => "This field categorizes mutations and polymorphisms. There are four possible values, DM, DP, DFP, and FP. DM are disease-causing mutations, pathological mutations reported to be disease causing in the original literature report. The other three tags are used for polymorphisms. DP are disease-associated polymorphisms. These are reported to be in significant association with disease (p<0.05) and are assumed to be functional (e.g. as a consequence of location, evolutionary conservation, replication studies etc), although there may be as yet no direct evidence (e.g. from an expression study) of function. DFP are disease-associated polymorphisms with additional supporting functional evidence. These are reported to be in significant association with disease (p<0.05) and to have evidence for being of direct functional importance (e.g. as a consequence of altered expression, mRNA studies etc). FP are in vitro/laboratory or in vivo functional polymorphisms. These are reported to affect the structure, function or expression of the gene (or gene product), but with no disease association reported as yet. (HGMD $hgmd_version)",
        HGMD_base => "A one-letter code determining which class (and table) the mutation belongs to. D deletion, E amplet, G grosdel, I insertion, M mutation, N grosins, P complex, R prom, S splice, and X indel.(HGMD $hgmd_version)",
        HGMD_hgvs => "Composite HGVS nomenclature for the mutation. (HGMD $hgmd_version)",
        HGMD_codon => "The number of the altered codon mapped to the HGMD cDNA sequence. (HGMD $hgmd_version)",
        HGMD_amino => "This field is specific to single base pair substitutions and contains the description of the nucleotide change. This is presented in terms of a triplet change with an additional flanking base included if the mutated base lies in either the first or third position in the triplet. For example, TACg-TAT represents a change of the last nucleotide C in the triplet to a T. TGT-TAT represents a change of the middle nucleotide G to an A. Note that the triplet itself is in upper-case letters, while the additional flanking base is in lower-case. (HGMD $hgmd_version)",
        HGMD_deletion => "Deletions are presented in terms of the deleted bases in lower case plus, in upper case, 10 bp DNA sequence flanking both sides of the lesion. Intron/exon boundary information may be provided where identified (e.g. _I12E13_). The codon number in the CODON field represents the last whole codon preceding the deletion, and is marked in the given sequence by the caret character (^). (HGMD $hgmd_version)",
        HGMD_insertion => "Insertions are presented in terms of the inserted bases in lower case plus, in upper case, 10 bp DNA sequence flanking both sides of the lesion. The numbered codon from the AMINO field is preceded in the given sequence by the caret character (^). (HGMD $hgmd_version)",
        HGMD_descr => "A narrative text description of the lesion. (HGMD $hgmd_version)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    
    my @DISEASES_PHENOTYPES;
    
    my $hgmd_chromosome = $bvfoa_info{chrom};
    my $hgmd_coordSTART = $bvfoa_info{chrom_start};
    my $hgmd_coordEND = $bvfoa_info{chrom_start}<=$bvfoa_info{chrom_end} ? $bvfoa_info{chrom_end} : $bvfoa_info{chrom_start}+1;
    my $hgmd_strand = $bvfoa_info{chrom_strand}==1 ? '+' : '-';

    my ($host,$port,$user,$password,$platform,$database) = @{$self->{params}};

    if (defined $bvfoa_info{hgnc}){
        my @data;
        
        my $query = "CALL $database.gene2hgmd_disease('$bvfoa_info{hgnc}')";
        my $qh = $self->{conn}->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            if(defined($row[0])){
                push @data, $row[0];
                push(@DISEASES_PHENOTYPES, split(/\|/,$row[0]));
            }
        }
        $line_hash->{HGMD_gene_diseases} = @data ? join('|', @data) : '';
    }

    my @HGMD_acc_num;
    my @HGMD_locus_disease;
    my @HGMD_tag;
    my @HGMD_base;
    my @HGMD_hgvs;
    my @HGMD_codon;
    my @HGMD_amino;
    my @HGMD_deletion;
    my @HGMD_insertion;
    my @HGMD_descr;
    
    my $query = "CALL $database.coord2hgmd('$hgmd_chromosome',$hgmd_coordSTART,$hgmd_coordEND)";
    my $qh = $self->{conn}->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        push @HGMD_acc_num, defined($row[0]) ? $row[0] : '';
        push @HGMD_locus_disease, defined($row[1]) ? $row[1] : '';
        push @HGMD_tag, defined($row[2]) ? $row[2] : '';
        push @HGMD_base, defined($row[3]) ? $row[3] : '';
        push @HGMD_hgvs, defined($row[4]) ? $row[4] : '';
        push @HGMD_codon, defined($row[5]) ? $row[5] : '';
        push @HGMD_amino, defined($row[6]) ? $row[6] : '';
        push @HGMD_deletion, defined($row[7]) ? $row[7] : '';
        push @HGMD_insertion, defined($row[8]) ? $row[8] : '';
        push @HGMD_descr, defined($row[9]) ? $row[9] : '';
        
        if(defined($row[1])){
            push(@DISEASES_PHENOTYPES, split(/\|/,$row[1]));
        }
    }
    $line_hash->{HGMD_acc_num} = @HGMD_acc_num ? join('|', @HGMD_acc_num) : '';
    $line_hash->{HGMD_locus_disease} = @HGMD_locus_disease ? join('|', @HGMD_locus_disease) : '';
    $line_hash->{HGMD_tag} = @HGMD_tag ? join('|', @HGMD_tag) : '';
    $line_hash->{HGMD_base} = @HGMD_base ? join('|', @HGMD_base) : '';
    $line_hash->{HGMD_hgvs} = @HGMD_hgvs ? join('|', @HGMD_hgvs) : '';
    $line_hash->{HGMD_codon} = @HGMD_codon ? join('|', @HGMD_codon) : '';
    $line_hash->{HGMD_amino} = @HGMD_amino ? join('|', @HGMD_amino) : '';
    $line_hash->{HGMD_deletion} = @HGMD_deletion ? join('|', @HGMD_deletion) : '';
    $line_hash->{HGMD_insertion} = @HGMD_insertion ? join('|', @HGMD_insertion) : '';
    $line_hash->{HGMD_descr} = @HGMD_descr ? join('|', @HGMD_descr) : '';
    
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

