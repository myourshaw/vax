=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      
 This is a work derived in part from the NGS-SNP annotate_SNPs.pl script, 
 Grant JR, Arantes AS, Liao X, Stothard P., In-depth annotation of SNPs arising from resequencing projects using NGS-SNP, Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.                                                                      
 and is licensed as a whole at no charge to all third parties under the terms of the original GNU GENERAL PUBLIC LICENSE.

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Phenotypes

=head1 SYNOPSIS

 mv Phenotypes.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Phenotypes

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new column:
 Phenotypes_gene, Phenotypes_locus, COSMIC_phenotypes_gene, COSMIC_phenotypes_locus, HGMD_PUBLIC_phenotypes_gene, HGMD_PUBLIC_phenotypes_locus.
 
 Requires that the VAX.pm module be in the Plugins directory
 
 References:

 (1) Grant JR, Arantes AS, Liao X, Stothard P. 
     In-depth annotation of SNPs arising from resequencing projects using NGS-SNP
     Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.

=head1 PARAMETERS

    parameters are in the form key=value, comma-separated
    'compara_db=s', #the name of the compara database to be used, default 'Multi'
    'model=s', #the model species to use when filling in the 'Model_Annotations' column in the output,default value is 'homo_sapiens'

=cut

package Phenotypes;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_bvfoa_info get_unique lstrip trim);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->{config};
    my %params;
    foreach my $param(@{$self->{params}}){
        my ($key,$value) = split(/=/,$param,1);
        $params{$key} = $value;
    }
    $params{model} ||= $config->{species};
    $params{compara_db} ||= 'Multi';
    $self->{params} = \%params;
    
    #adaptors used by NGS-SNP
    $self->{ma} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Member');
    $self->{ha} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Homology');
    #$self->{fa} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Family');
    $self->{goa} = $config->{reg}->get_adaptor($params{compara_db}, 'Ontology', 'GOTerm');
    #$self->{translation_adaptor} = $config->{reg}->get_adaptor($config->{species}, 'core', 'translation');
    #these adaptors are used for the Model_Annotations field.
    #if they cannot be created then this field is not filled in.
    if (defined($params{model})) {
        $self->{model_translation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'translation');
        $self->{model_transcript_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'transcript');
        #$self->{model_gene_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'gene');
        $self->{model_variation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variation');
        $self->{model_variationfeature_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variationfeature');
        $self->{model_slice_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'slice');
        #$self->{model_variationannotation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variationannotation');
        $self->{model_phenotypefeature_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'phenotypefeature');

        #if (   ( !defined( $config->{model_translation_adaptor} ) )
        #    || ( !defined( $config->{model_transcript_adaptor} ) )
        ##    || ( !defined( $config->{model_gene_adaptor} ) )
        #    || ( !defined( $config->{model_variation_adaptor} ) )
        #    || ( !defined( $config->{model_variationfeature_adaptor} ) )
        #    || ( !defined( $config->{model_slice_adaptor} ) )
        #    || ( !defined( $config->{model_variationannotation_adaptor} ) ) )
        #{
        #    message( $config->{verbose}, $config->{log_file},
        #        "Unable to get adaptors for the species '$config->{model}' specified using the '-model' option.\n"
        #    );
        #}
    }

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
        Phenotypes_gene
        Phenotypes_locus
        COSMIC_phenotypes_gene
        COSMIC_phenotypes_locus
        HGMD_PUBLIC_phenotypes_gene
        HGMD_PUBLIC_phenotypes_locus
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Phenotypes_gene => "phenotype associated with gene. (| separated) (Ensembl)",
        Phenotypes_locus => "phenotypic information associated with known variation in the model species at the variant locus. If the input SNPs are from the model species then their locations are used to identify known variations at the same locations, and any phenotypic information linked to these variations is reported. If the input SNPs are not from the model species and the input SNP alters a protein, then protein alignment is used to find the orthologous genomic region from the model. Model species variations that alter this region are obtained, and any phenotypic information linked to these variations is reported. (| separated) (Ensembl)",
        COSMIC_phenotypes_gene => "cancer phenotype associated with gene in COSMIC database. (| separated) (Ensembl)",
        COSMIC_phenotypes_locus => "cancer phenotypic information associated with known variation in the model species at the variant locus in COSMIC database. (| separated) (Ensembl)",
        HGMD_PUBLIC_phenotypes_gene => "ID of phenotype associated with gene in HGMD-PUBLIC database. (| separated) (Ensembl)",
        HGMD_PUBLIC_phenotypes_locus => "ID of phenotypic information associated with known variation in the model species at the variant locus in HGMD-PUBLIC database. (| separated) (Ensembl)",
        };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my $config = $self->{config};
    my %params = %{$self->{params}};
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    
    my @DISEASES_PHENOTYPES;
    
    #get phenotypes based on gene
    if(defined($bvfoa_info{hgnc})){
        my ($phenotype_gene_annotations, $cosmic_phenotype_gene_annotations, $hgmd_public_phenotype_gene_annotations) =
            get_gene_phenotypes_by_hgnc($self, $bvfoa_info{hgnc});
        if ((defined($phenotype_gene_annotations))
            && (scalar(@{$phenotype_gene_annotations}) > 0)
            )
        {
            s/[\|\;]//g for (@{$phenotype_gene_annotations});
            my $phenotypes_gene = join( '|', @{$phenotype_gene_annotations});
            
            push(@DISEASES_PHENOTYPES, @{$phenotype_gene_annotations});
            $line_hash->{Phenotypes_gene} = $phenotypes_gene;
        }
        if ((defined($cosmic_phenotype_gene_annotations))
            && (scalar(@{$cosmic_phenotype_gene_annotations}) > 0)
            )
        {
            s/[\|\;]//g for (@{$cosmic_phenotype_gene_annotations});
            my $cosmic_phenotypes_gene = join( '|', @{$cosmic_phenotype_gene_annotations});
            $line_hash->{COSMIC_phenotypes_gene} = $cosmic_phenotypes_gene;
        }
        if ((defined($hgmd_public_phenotype_gene_annotations))
            && (scalar(@{$hgmd_public_phenotype_gene_annotations}) > 0)
            )
        {
            s/[\|\;]//g for (@{$hgmd_public_phenotype_gene_annotations});
            my $hgmd_public_phenotypes_gene = join( '|', @{$hgmd_public_phenotype_gene_annotations});
            $line_hash->{HGMD_PUBLIC_phenotypes_gene} = $hgmd_public_phenotypes_gene;
        }
    }

    #get phenotypes base on locus
    #using a model organism as a source of phenotypes may be possible  BUT IS COMPLETELY UNTESTED
    if(defined($bvfoa_info{translation})
       && defined($bvfoa_info{ensg})
       && defined($bvfoa_info{ensp})
       && defined($bvfoa_info{altered_aa_start})
       && defined($bvfoa_info{altered_aa_end})
       && defined($bvfoa_info{genomic_coords})
       ){

        #get a list of aligned orthologous resides from model species
        #will be undef if species and model are same
        #results are used in get_overlapping_protein_features_from_UniProt_orthologues
        #and get_model_phenotypes
        my $aligned_model_orthologue_proteins =
            get_model_orthologous_residues(
                $self,
                $bvfoa_info{ensg},
                $bvfoa_info{ensp},
                $bvfoa_info{altered_aa_start},
                $bvfoa_info{altered_aa_end}
            );
        #Phenotypic information associated with known variation at site (or orthologous site) in the model species
        #using information from Ensembl
        if (defined($params{model})) {
            my ($model_phenotype_annotations, $cosmic_model_phenotype_annotations, $hgmd_public_model_phenotype_annotations) =
                get_model_phenotypes(
                    $self,
                    $bvfoa_info{genomic_coords},
                    $bvfoa_info{ensp},
                    $bvfoa_info{altered_aa_start},
                    $bvfoa_info{altered_aa_end},
                    $bvfoa_info{aligned_model_orthologue_proteins}
                    );
            if ((defined($model_phenotype_annotations))
                && (scalar(@{$model_phenotype_annotations}) > 0)
                )
            {
                s/[\|\;]//g for (@{$model_phenotype_annotations});
                my $phenotypes_locus = join( '|', map {@$_} @{$model_phenotype_annotations});
                push(@DISEASES_PHENOTYPES, map {@$_} @{$model_phenotype_annotations});
                $line_hash->{Phenotypes_locus} = $phenotypes_locus;
            }
            if ((defined($cosmic_model_phenotype_annotations))
                && (scalar(@{$cosmic_model_phenotype_annotations}) > 0)
                )
            {
                s/[\|\;]//g for (@{$cosmic_model_phenotype_annotations});
                my $cosmic_phenotypes_locus = join( '|', @{$cosmic_model_phenotype_annotations});
                $line_hash->{COSMIC_phenotypes_locus} = $cosmic_phenotypes_locus;
            }
            if ((defined($hgmd_public_model_phenotype_annotations))
                && (scalar(@{$hgmd_public_model_phenotype_annotations}) > 0)
                )
            {
                s/[\|\;]//g for (@{$hgmd_public_model_phenotype_annotations});
                my $hgmd_public_phenotypes_locus = join( '|', @{$hgmd_public_model_phenotype_annotations});
                $line_hash->{HGMD_PUBLIC_phenotypes_locus} = $hgmd_public_phenotypes_locus;
            }
        }
    } #if(defined($translation))
    
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

sub get_gene_phenotypes_by_hgnc {
#returns array of arrays of phenotypes associated with gene
    my ($self, $hgnc) = @_;

    my @phenotypes = ();

    #my $vas = $self->{model_variationannotation_adaptor}
    #    ->fetch_all_by_associated_gene($hgnc);
    my $pfs = $self->{model_phenotypefeature_adaptor}
        ->fetch_all_by_associated_gene($hgnc);

    #return get_phenotypes_from_variant_annotations($vas);
    return get_phenotypes_from_phenotype_features($pfs);
}

sub get_model_phenotypes {
#adds information about known model phenotypes to consequence
#if model and species are the same: get_model_phenotypes_from_model_genomic_coords
#   undef if model and species are the same ?
#else: get_model_genomic_coords_of_model_amino_acid, get_model_phenotypes_from_model_genomic_coords
    my ($self, $genomic_coords, $ensp, $altered_aa_start, $altered_aa_end, $aligned_model_orthologue_proteins) = @_;

    #my @model_phenotype_annotations;
    #my @cosmic_model_phenotype_annotations;
    #my @hgmd_public_model_phenotype_annotations;
    my @aligned_phenotypes = ();
    my @cosmic_aligned_phenotypes = ();
    my @hgmd_public_aligned_phenotypes = ();

    if ((!defined($self->{model_translation_adaptor}))
        || (!defined( $self->{model_transcript_adaptor}))
        || (!defined( $self->{model_variation_adaptor}))
        || (!defined( $self->{model_variationfeature_adaptor}))
        || (!defined( $self->{model_slice_adaptor}))
        #|| (!defined( $self->{model_variationannotation_adaptor})))
        || (!defined( $self->{model_phenotypefeature_adaptor})))
    {
        return;
    }

    #if input SNPs are from model then look for known variants
    #and phenotypes at that position
    if (lc( $self->{config}->{species}) eq lc($self->{params}->{model})) {
        #return get_model_phenotypes_from_model_genomic_coords($self, $genomic_coords);
        my ($phenotypes, $cosmic_phenotypes, $hgmd_public_phenotypes)
            = get_model_phenotypes_from_model_genomic_coords($self, $genomic_coords);
        if(scalar(@{$phenotypes}) > 0){
            push(@aligned_phenotypes,$phenotypes);
        }
        if(scalar(@{$cosmic_phenotypes}) > 0){
            push(@cosmic_aligned_phenotypes,@$cosmic_phenotypes);
        }
        if(scalar(@{$hgmd_public_phenotypes}) > 0){
            push(@hgmd_public_aligned_phenotypes,@$hgmd_public_phenotypes);
        }
    }

    #otherwise transform non-model coordinate to model
    else {
        if ((!defined($ensp))
            || (!defined($altered_aa_start))
            || (!defined($altered_aa_end)))
        {
            return;
        }

        foreach my $model_orthologous_residue (
            @{$aligned_model_orthologue_proteins})
        {
            my $genomic_coords
                = get_model_genomic_coords_of_model_amino_acid($self, $model_orthologous_residue);
            if (defined($genomic_coords)) {
                my ($phenotypes, $cosmic_phenotypes, $hgmd_public_phenotypes)
                    = get_model_phenotypes_from_model_genomic_coords($self, $genomic_coords);
                if(scalar(@{$phenotypes}) > 0){
                    @aligned_phenotypes = (@aligned_phenotypes,$phenotypes);
                }
                if(scalar(@{$cosmic_phenotypes}) > 0){
                    @cosmic_aligned_phenotypes = (@cosmic_aligned_phenotypes,$cosmic_phenotypes);
                }
                if(scalar(@{$hgmd_public_phenotypes}) > 0){
                    @hgmd_public_aligned_phenotypes = (@hgmd_public_aligned_phenotypes,$hgmd_public_phenotypes);
                }
            }
        }
    }
    
    return (\@{get_unique(\@aligned_phenotypes)}, \@{get_unique(\@cosmic_aligned_phenotypes)}, \@{get_unique(\@hgmd_public_aligned_phenotypes)});

    #
    #my $unique_phenotypes = get_unique( \@aligned_phenotypes );
    #push(
    #    @model_phenotype_annotations,
    #    @{$unique_phenotypes}
    #);
    #return \@model_phenotype_annotations;
}

sub get_model_orthologous_residues {
#returns array of hashes describing model orthologous residues
#apparently this doesn't do anything if species and model are the same
    my ($self, $gene_id, $ensp, $altered_aa_start, $altered_aa_end) = @_;
    if ((lc( $self->{config}->{species}) eq lc($self->{params}->{model}))
        || (!defined($gene_id))
        || (!defined($ensp))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end))
        || (!defined( $self->{model_translation_adaptor})))
    {
        return;
    }
 
    #return value
    my @model_orthologous_residues = ();

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    my $member = $self->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];

    #fetch the homology relationships where the given member is implicated
    #in pair with another member from the paired species. Member species and
    #paired species should be different
    push(
        @{$homologies},
        @{$self->{ha}->fetch_all_by_Member_paired_species($member,
                $self->{params}->{model}, ['ENSEMBL_ORTHOLOGUES'])}
    );

    foreach my $homology (@{$homologies}) {

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
     #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef;

        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if ( !( defined($aln) ) ) {
            next;
        }

        #confirm that reference protein is in the alignment
        if (!(scalar($aln->each_seq_with_id($ensp)) == 1))
        {
            next;
        }

#TODO:verify indel
        #get the column containing the residue encoded by the SNP-containing codon
        my $col = undef;
        my $col_end = undef;

        #This function gives the position in the alignment
        #(i.e. column number) of the given residue number in the
        #sequence with the given name.
        #column_from_residue_number can throw an exception if the altered residue
        #in the reference is the one that encodes the stop codon
        eval {
            $col = $aln->column_from_residue_number(
                $ensp,
                $altered_aa_start
            );
        };
        if ($@) {
            next;
        }
        eval {
            $col_end = $aln->column_from_residue_number(
                $ensp,
                $altered_aa_end
            );
        };
        if ($@) {
            next;
        }

        #Bio::SimpleAlign->slice
        #Creates a slice from the alignment inclusive of start and
        #end columns.  Sequences with no residues in the slice are
        #excluded from the new alignment and a warning is printed.
        #Slice beyond the length of the sequence does not do
        #padding.
        #get an alignment containing only the desired column
        #my $sub_align = $aln->slice( $col, $col );
        my $sub_align = undef;
        #prevent annoying warnings; seems to require Bioperl-live, not Bioperl-1.2.3
        if (!defined($self->{config}->{verbose})){
            $SIG{'__WARN__'} = sub {return 1};
        }
        eval { $sub_align = $aln->slice($col, $col_end); };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        if (!( scalar( $sub_align->each_seq()) == 2)) {
            next;
        }

        my $seq1 = $sub_align->get_seq_by_pos(1);
        my $seq2 = $sub_align->get_seq_by_pos(2);

        #in all the examples I've examined, $seq1 is the reference sequence
        if ($seq1->display_id eq $ensp) {

            #want to get Entrez Gene ID from gene object
            my $model_gene = $$homologues[1]->get_Gene;
            my $entrez_gene_id;
            my $refseq_protein_id;
            foreach my $xref (@{$model_gene->get_all_DBLinks()}) {
                if ($xref->dbname() eq 'EntrezGene') {
                    $entrez_gene_id = $xref->primary_id();
                }
                elsif ($xref->dbname() eq 'RefSeq_peptide') {
                    $refseq_protein_id = $xref->display_id();
                }
            }

            #obtain the the Ensembl ID of the model protein
            my $model_translation_ID = $seq2->display_id;

            #obtain the position in the model protein of the amino acid that aligns with the SNP-altered residue
            my $model_aligned_aa_start = $seq2->start;
            my $model_aligned_aa_end   = $seq2->end;

            my %orthologue = (
                id             => $model_translation_ID,
                start          => $model_aligned_aa_start,
                end            => $model_aligned_aa_end,
                sequence       => undef,
                entrez_gene_id => $entrez_gene_id,
                refseq_peptide => undef,   #will only be used to try to get NCBI Gene record
                uniprot_id => [],
                go         => undef
            );

            #attempt to fill in usefull IDs and other info
            my $model_translation_object
                = $self->{model_translation_adaptor}
                ->fetch_by_stable_id($model_translation_ID);
            $orthologue{sequence} = $model_translation_object->seq();

            my $xrefs = $model_translation_object->get_all_DBLinks();

            my @go_names = ();
            foreach my $xref (@{$xrefs}) {
                if ($xref->dbname() eq 'Uniprot/SWISSPROT') {
                    push(@{$orthologue{uniprot_id}}, $xref->display_id());
                }
                elsif ($xref->dbname() eq 'RefSeq_peptide') {
                    $orthologue{refseq_peptide} = $xref->display_id();
                }
                elsif ($xref->dbname() eq 'goslim_goa') {
                    my $go_id = $xref->display_id();
                    my $go_term;
                    my $go_name;
                    my $go_definition;
                    if (defined($self->{goa})) {
                        $go_term = $self->{goa}->fetch_by_accession($go_id);
                        $go_name       = $go_term->name();
                        $go_definition = $go_term->definition();
                    }

                    if (defined($go_name)) {
                        push(@go_names, "[$go_id]:$go_name");
                    }
                    else {
                        push(@go_names, "$go_id");
                    }
                }
            }
            $orthologue{go} = get_unique( \@go_names );

            if (!defined($orthologue{refseq_peptide})) {
                $orthologue{refseq_peptide} = $refseq_protein_id;
            }

            push(@model_orthologous_residues, \%orthologue);

        }
        else {
            die("Unexpected sequence ordering in alignment.");
        }
    }
    return \@model_orthologous_residues;
}

sub get_model_phenotypes_from_model_genomic_coords {
#returns an array of arrays of phenotypes associated with model genomic coordinates
    my ($self, $genomic_coords) = @_;

    my @phenotypes = ();
    my @cosmic_phenotypes = ();
    my @hgmd_public_phenotypes = ();

    my $slice = $self->{model_slice_adaptor}->fetch_by_region(
        undef,
        $genomic_coords->{chr},
        $genomic_coords->{start},
        $genomic_coords->{end},
        $genomic_coords->{strand}
    );

    my $vfs = $self->{model_variationfeature_adaptor}
        ->fetch_all_by_Slice($slice);
    my $somatic_vfs = $self->{model_variationfeature_adaptor}->fetch_all_somatic_by_Slice($slice);
    my @all_vfs = (@{$vfs}, @{$somatic_vfs});

    foreach my $vf ( @all_vfs ) {

        my $variation_name = $vf->variation_name();
        my $variation = $self->{model_variation_adaptor}
            ->fetch_by_name($variation_name);

        if(ref($variation) && $variation->isa('Bio::EnsEMBL::Variation::Variation')) {
            my $pfs = $self->{model_phenotypefeature_adaptor}
                ->fetch_all_by_Variation($variation);
            my ($phenotypes_pfs, $cosmic_phenotypes_pfs, $hgmd_public_phenotypes_pfs)
                = get_phenotypes_from_phenotype_features($pfs);
            if (scalar(@{$phenotypes_pfs}) > 0) {
                    @phenotypes = (@phenotypes, @{$phenotypes_pfs});
            }
            if (scalar(@{$cosmic_phenotypes_pfs}) > 0) {
                    @cosmic_phenotypes = (@cosmic_phenotypes, @{$cosmic_phenotypes_pfs});
            }
            if (scalar(@{$hgmd_public_phenotypes_pfs}) > 0) {
                    @hgmd_public_phenotypes = (@hgmd_public_phenotypes, @{$hgmd_public_phenotypes_pfs});
            }
        }
    }

    return (\@{get_unique(\@phenotypes)}, \@{get_unique(\@cosmic_phenotypes)}, \@{get_unique(\@hgmd_public_phenotypes)});

}

sub get_model_genomic_coords_of_model_amino_acid {
#returns hash containing chromosome, start, end, and strand of genomic sequence encoding
#amino acid described by $model_translation_id, $amino_acid_start, $amino_acid_end
    my ($self, $model_orthologous_residue) = @_;

    my $model_translation_id = $model_orthologous_residue->{id};
    my $amino_acid_start     = $model_orthologous_residue->{start};
    my $amino_acid_end       = $model_orthologous_residue->{end};

    #obtain the model protein translation object
    my $model_translation = $self->{model_translation_adaptor}
        ->fetch_by_stable_id($model_translation_id);
    my $model_transcript = $self->{model_transcript_adaptor}
        ->fetch_by_translation_stable_id($model_translation_id);

    #    print "\$amino_acid_start is $amino_acid_start\n";
    #    print "\$amino_acid_end is $amino_acid_end\n";

    #determine the position on the transcript
    my $trmapper = $model_transcript->get_TranscriptMapper();

    my $model_slice = $self->{model_slice_adaptor}
        ->fetch_by_transcript_stable_id($model_transcript->stable_id);

    my $model_chromosome_name = $model_slice->seq_region_name();
    my $model_slice_start     = $model_slice->start();
    my $model_slice_end       = $model_slice->end();
    my $model_slice_strand    = $model_slice->strand();

#a list of Bio::EnsEMBL::Mapper::Gap and Bio::EnsEMBL::Mapper::Coordinate objects
    my @model_genomic_coords
        = $trmapper->pep2genomic($amino_acid_start, $amino_acid_end);

    #not sure how to handle more complicated coordinates, so skip
    if (scalar(@model_genomic_coords) > 1) {
        return;
    }

    my $model_aligned_aa_start_genomic = $model_genomic_coords[0]->start;
    my $model_aligned_aa_end_genomic   = $model_genomic_coords[0]->end;
    my $model_aligned_strand_genomic   = $model_genomic_coords[0]->strand;

#check that the slice start and end is consistent with the position obtained using pep2genomic
    if ((!( $model_slice_start <= $model_aligned_aa_start_genomic))
        || (!($model_slice_end >= $model_aligned_aa_end_genomic)))
    {
        return;
    }

    my %model_genomic_cords = (
        chr    => $model_chromosome_name,
        start  => $model_aligned_aa_start_genomic,
        end    => $model_aligned_aa_end_genomic,
        strand => $model_aligned_strand_genomic
    );
    return \%model_genomic_cords;

}


sub get_phenotypes_from_variant_annotations{
#from a list of VariantAnnotation objects, returns lists of formatted strings with
#phenotype information for (@phenotypes, @cosmic_phenotypes, @hgmd_public_phenotypes)

    my $vas = shift;
    
    my %pheno; #key = phenotype_description[ (phenotype_name)]; value = list of pheno_sources
    my %pheno_source = (); #key = source; value = list of source_name variation_names[,p=p_value][,risk allele=risk_allele]
    my @phenotypes = ();
    my @phenotypes_short = (); #phenotype_description[ (phenotype_name)] only
    my %cosmic = (); #key = tumor site from phenotype_description; value = list of variation_names
    my @cosmic_phenotypes = ();
    my @cosmic_phenotypes_short = (); #tumor_site only
    my @hgmd_public_phenotypes = (); # list of variation_names
    my @hgmd_public_phenotypes_short = (); # empty

    foreach my $va (@{$vas}) {

        if (defined($va->source_name()) && $va->source_name() eq 'COSMIC'){
            if (defined($va->phenotype_description()) && defined( $va->variation_names())) {
                push(@{$cosmic{lstrip(lstrip(trim($va->phenotype_description()),'COSMIC:'),'tumour_site:')}}, trim($va->variation_names()));
                push(@cosmic_phenotypes_short, lstrip(lstrip(trim($va->phenotype_description()),'COSMIC:'),'tumour_site:'));
            }
        }
        elsif (defined($va->source_name()) && $va->source_name() eq 'HGMD-PUBLIC'){
            if (defined( $va->variation_names())) {
                push(@hgmd_public_phenotypes, trim($va->variation_names()));
                push(@hgmd_public_phenotypes_short, trim($va->variation_names()));
            }
        }
        else{
            my $key = undef;
            my @value = ();
            if (defined($va->phenotype_description()) && !defined($va->phenotype_name())) { #phenotype_name always null? at least in human
                $key = trim($va->phenotype_description());
            }
            elsif(defined($va->phenotype_description()) && defined($va->phenotype_name())){
                $key = trim($va->phenotype_description()).' ('.trim($va->phenotype_name()).')';
            }
            elsif (defined($va->phenotype_name())) {
                $key = trim($va->phenotype_name());
            }
            if(defined($key)){
                push(@phenotypes_short, $key);
                my $source_name = defined($va->source_name()) ? trim($va->source_name()) : '?';
                if (defined( $va->variation_names())) {
                    push(@value,  trim($va->variation_names()));
                }
                if (defined($va->p_value)){
                    push(@value, 'p='.trim($va->p_value()));
                }
                if (defined($va->associated_variant_risk_allele)){
                    my @risk_allele_array = split(/\-/, $va->associated_variant_risk_allele);
                    if (scalar(@risk_allele_array)>1){
                        push(@value, 'risk_allele='.trim($risk_allele_array[1]));
                    }
                }
                push(@{$pheno_source{$source_name}}, join(' ', @value));
                $pheno{$key} = \%pheno_source;
            }
        }
    }
    map {push(@cosmic_phenotypes, $_.'['.join(',',@{get_unique(\@{$cosmic{$_}})}).']')} keys %cosmic;
    foreach my $description(keys %pheno){
        my @p;
        foreach my $source(keys $pheno{$description}){
            push(@p, $source.'('.join(',',@{get_unique($pheno{$description}->{$source})}).')');
        }
        push(@phenotypes,$description.'['.join(';', @p).']');
    }
    #uncomment to return lengthy details
    #return (\@{get_unique(\@phenotypes)}, \@{get_unique(\@cosmic_phenotypes)}, \@{get_unique(\@hgmd_public_phenotypes)});    
    return (\@{get_unique(\@phenotypes_short)}, \@{get_unique(\@cosmic_phenotypes_short)}, \@{get_unique(\@hgmd_public_phenotypes_short)});    
}

sub get_phenotypes_from_phenotype_features{
#from a list of PhenotypeFeature objects, returns lists of formatted strings with
#phenotype information for (@phenotypes, @cosmic_phenotypes, @hgmd_public_phenotypes)

    my $pfs = shift;
    
    my %pheno; #key = phenotype_description[ (phenotype_name)]; value = list of pheno_sources
    my %pheno_source = (); #key = source; value = list of source_name variation_names[,p=p_value][,risk allele=risk_allele]
    my @phenotypes = ();
    my @phenotypes_short = (); #phenotype_description[ (phenotype_name)] only
    my %cosmic = (); #key = tumor site from phenotype_description; value = _object_id (COSMnnnn)
    my @cosmic_phenotypes = ();
    my @cosmic_phenotypes_short = (); #tumor_site only
    my @hgmd_public_phenotypes = (); # list of variation_names
    my @hgmd_public_phenotypes_short = (); # empty

    foreach my $pf (@{$pfs}) {
        my $phenotype = $pf->phenotype;
        if (defined($pf->source()) && $pf->source() eq 'COSMIC'){
            if (defined($phenotype->description()) && defined( $pf->{_object_id})) {
                push(@{$cosmic{lstrip(lstrip(trim($phenotype->description()),'COSMIC:'),'tumour_site:')}}, trim($pf->{_object_id}));
                push(@cosmic_phenotypes_short, lstrip(lstrip(trim($phenotype->description()),'COSMIC:'),'tumour_site:'));
            }
        }
        elsif (defined($pf->source()) && $pf->source() eq 'HGMD-PUBLIC'){
            if (defined( $pf->{_object_id})) {
                push(@hgmd_public_phenotypes, trim($pf->{_object_id}));
                push(@hgmd_public_phenotypes_short, trim($pf->{_object_id}));
            }
        }
        else{
            my $key = undef;
            my @value = ();
            if (defined($phenotype->description()) && !defined($phenotype->name())) { #phenotype_name always undef? at least in human
                $key = trim($phenotype->description());
            }
            elsif(defined($phenotype->description()) && defined($phenotype->name())){
                $key = trim($phenotype->description()).' ('.trim($phenotype->name()).')';
            }
            elsif (defined($phenotype->name())) {
                $key = trim($phenotype->name());
            }
            if(defined($key)){
                push(@phenotypes_short, $key);
                my $source_name = defined($pf->source()) ? trim($pf->source()) : '?';
                #variation_names are rs numbers
                if (defined( $pf->variation_names())) {
                    push(@value,  trim($pf->variation_names()));
                }
                if (defined($pf->p_value)){
                    push(@value, 'p='.trim($pf->p_value()));
                }
                if (defined($pf->risk_allele())){
                    my @risk_allele_array = split(/\-/, $pf->risk_allele());
                    if (scalar(@risk_allele_array)>1){
                        push(@value, 'risk_allele='.trim($risk_allele_array[1]));
                    }
                }
                push(@{$pheno_source{$source_name}}, join(' ', @value));
                $pheno{$key} = \%pheno_source;
            }
        }
    }
    map {push(@cosmic_phenotypes, $_.'['.join(',',@{get_unique(\@{$cosmic{$_}})}).']')} keys %cosmic;
    foreach my $description(keys %pheno){
        my @p;
        foreach my $source(keys $pheno{$description}){
            push(@p, $source.'('.join(',',@{get_unique($pheno{$description}->{$source})}).')');
        }
        push(@phenotypes,$description.'['.join(';', @p).']');
    }
    #uncomment to return lengthy details
    #return (\@{get_unique(\@phenotypes)}, \@{get_unique(\@cosmic_phenotypes)}, \@{get_unique(\@hgmd_public_phenotypes)});    
    return (\@{get_unique(\@phenotypes_short)}, \@{get_unique(\@cosmic_phenotypes_short)}, \@{get_unique(\@hgmd_public_phenotypes_short)});    
    #return (\@{get_unique(\@phenotypes_short)}, \@{get_unique(\@cosmic_phenotypes)}, \@{get_unique(\@hgmd_public_phenotypes_short)});    
}

1;

