=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      
 This is a work derived from the NGS-SNP annotate_SNPs.pl script, 
 Grant JR, Arantes AS, Liao X, Stothard P., In-depth annotation of SNPs arising from resequencing projects using NGS-SNP, Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.                                                                      
 and is licensed as a whole at no charge to all third parties under the terms of the original GNU GENERAL PUBLIC LICENSE.

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Alignment

=head1 SYNOPSIS

 mv Alignment.pm vw.pm blosum62.mat ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin Alignment

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP)
 that adds the following new columns:
 Alignment_Score_Change, C_blosum, Context_Conservation, Amino_Acids_In_Orthologues, Orthologue_Species.
 
 WARNING: this is a database-intensive, slow-runnig plugin. Depending on the versions
 of the Ensembl API and BioPerl, it may generate many (probably harmless) warnings.
 
    Requires the vw.pm plugin in the Plugins directory

    Requires the VAX.pm module in the Plugins directory
    
    Requires Bio::Matrix::IO (Bio/Matrix/IO.pm)
        This BioPerl module is NOT available in BioPerL 1.2.3, which is claimed to be necessary for the ensembl API.
        Suggestion: install a current version of BioPerl and add to PERL5LIB with lower priority than bioperl 1.2.3
        
    Requires blosum62.mat or other scoring matrix in the same directory as this plugin
    
    Contents of blosum62.mat:
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 

 References:

 (1) Grant JR, Arantes AS, Liao X, Stothard P. 
     In-depth annotation of SNPs arising from resequencing projects using NGS-SNP
     Bioinformatics. 2011 Aug 15;27(16):2300-1. doi: 10.1093/bioinformatics/btr372. Epub 2011 Jun 22.
     
 (2) Kowarsch A, Fuchs A, Frishman D et al.
     Correlated mutations: a hallmark of phenotypic amino acid substitutions
     PLoS Comput Biol 2010;6.

=head1 PARAMETERS

    parameters are in the form key=value, comma-separated
    'flanking=i', #amount of flanking genomic sequence to write on each side of SNPs when -f option is used (Optional; default is 100 bases)
    'flanking_output=s', #write genomic flanking sequence for each SNP to this file (Optional)
    'comparison_species=s', #names of species (semicolon-separated) to use when assessing sequence conservation (, default all
    'compara_db=s', #the name of the compara database to be used, default 'Multi'
    'model=s', #the model species to use when filling in the 'Model_Annotations' column in the output,default value is 'homo_sapiens'
    'scoring_matrix=s', #the scoring matrix file to use for amino acid comparisons, default blosum62.mat
    'flanking_for_context_conservation=i', #amount of flanking protein sequence on either side of SNP-affected residue to use for determining conservation of region containing coding SNP, default 10
    'degenerate_flanking', #indicate known SNP sites within flanking sequence as lowercase IUPAC DNA bases (Optional; default is to return unmodified reference sequence)
    'is_fake_vfa', #
    
=cut

package Alignment;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info get_consequence_info get_taxon_name space_to_underscore is_gap);
use Bio::Matrix::IO;
use List::Util qw(first);
use File::Basename;
use LWP::Simple;
use HTML::TokeParser;

    #When Ensembl adds new species, add abbreviations to this table
    my %species_short = (
        Aedes_aegypti => 'Aa',
        Ailuropoda_melanoleuca => 'Am',
        Anas_platyrhynchos => 'Ap',
        Anolis_carolinensis => 'Ac',
        Anopheles_gambiae => 'Ag',
        Apis_mellifera => 'Ame',
        Bos_taurus => 'Bt',
        Caenorhabditis_elegans => 'Ce',
        Callithrix_jacchus => 'Cj',
        Canis_familiaris => 'Cf',
        Canis_lupus_familiaris => 'Clf',
        Cavia_porcellus => 'Cp',
        Choloepus_hoffmanni => 'Ch',
        Ciona_intestinalis => 'Ci',
        Ciona_savignyi => 'Cs',
        Culex_quinquefasciatus => 'Cq',
        Danio_rerio => 'Dr',
        Dasypus_novemcinctus => 'Dn',
        Dipodomys_ordii => 'Do',
        Drosophila_melanogaster => 'Dm',
        Echinops_telfairi => 'Et',
        Equus_caballus => 'Ec',
        Erinaceus_europaeus => 'Ee',
        Felis_catus => 'Fc',
        Ficedula_albicollis => 'Fa',
        Gadus_morhua => 'Gm',
        Gallus_gallus => 'Gga',
        Gasterosteus_aculeatus => 'Ga',
        Gorilla_gorilla => 'Gg',
        Gorilla_gorilla_gorilla => 'Ggg',
        Homo_sapiens => 'Hs',
        Ictidomys_tridecemlineatus => 'It',
        Latimeria_chalumnae => 'Lc',
        Loxodonta_africana => 'La',
        Macaca_mulatta => 'Mm',
        Macropus_eugenii => 'Me',
        Meleagris_gallopavo => 'Mg',
        Microcebus_murinus => 'Mmur',
        Monodelphis_domestica => 'Md',
        Mus_musculus => 'Mmus',
        Mustela_putorius_furo => 'Mpf',
        Myotis_lucifugus => 'Ml',
        Nomascus_leucogenys => 'Nl',
        Ochotona_princeps => 'Op',
        Oreochromis_niloticus => 'On',
        Ornithorhynchus_anatinus => 'Oa',
        Oryctolagus_cuniculus => 'Oc',
        Oryzias_latipes => 'Ol',
        Otolemur_garnettii => 'Og',
        Ovis_aries => 'Oar',
        Pan_troglodytes => 'Pt',
        Papio_hamadryas => 'Ph',
        Pelodiscus_sinensis => 'Ps',
        Petromyzon_marinus => 'Pm',
        Pongo_abelii => 'Pa',
        Procavia_capensis => 'Pc',
        Pteropus_vampyrus => 'Pv',
        Rattus_norvegicus => 'Rn',
        Saccharomyces_cerevisiae => 'Sc',
        Sarcophilus_harrisii => 'Sh',
        Sorex_araneus => 'Sa',
        Spermophilus_tridecemlineatus => 'St',
        Sus_scrofa => 'Ss',
        Taeniopygia_guttata => 'Tg',
        Takifugu_rubripes => 'Tr',
        Tarsius_syrichta => 'Ts',
        Tetraodon_nigroviridis => 'Tn',
        Tupaia_belangeri => 'Tb',
        Tursiops_truncatus => 'Tt',
        Vicugna_pacos => 'Vp',
        Xenopus_tropicalis => 'Xt',
        Xiphophorus_maculatus => 'Xm',
    );

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->{config};
    my %params;
    foreach my $param(@{$self->{params}}){
        my ($key,$value) = split(/=/,$param,1);
        $params{$key} = $value;
    }
    $params{flanking} ||= 100;
    $params{flanking_output} ||= undef;
    $params{compara_db} ||= 'Multi';
    $params{model} ||= 'homo_sapiens';
    $params{scoring_matrix} ||= File::Spec->catfile(dirname(__FILE__), 'blosum62.mat');
    $params{flanking_for_context_conservation} ||= 10;
    $params{degenerate_flanking} ||= undef;
    $params{is_fake_vfa} ||= 0;
    if (defined($params{comparison_species})){
        my %comparison_species;
        map {$comparison_species{$_} = ''} split(/;/,$params{comparison_species});
        $params{comparison_species} = \%comparison_species;
    }
    $self->{params} = \%params;
    
    #get scoring matrix object
    my $parser = Bio::Matrix::IO->new(
        -format => 'scoring',
        -file   => $params{scoring_matrix}
    );
    $self->{matrix} = $parser->next_matrix;
    #determine maximum possible alignment score change value, for normalization
    #defaults to 15
    $self->{max_alignment_score_change} = get_max_alignment_score_change($self->{matrix});

    #adaptors used by NGS-SNP
    #direct access to Member is deprecated
    #$self->{ma} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Member');
    $self->{sma} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'SeqMember');
    $self->{gma} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'GeneMember');
    $self->{ha} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Homology');
    #$self->{fa} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Family');
    $self->{goa} = $config->{reg}->get_adaptor($params{compara_db}, 'Ontology', 'GOTerm');
    #$self->{translation_adaptor} = $config->{reg}->get_adaptor($config->{species}, 'core', 'translation');
    #these adaptors are used for the Model_Annotations field.
    #if they cannot be created then this field is not filled in.
    if (defined($params{model})) {
        $self->{model_translation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'translation');
        #$self->{model_transcript_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'transcript');
        #$self->{model_gene_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'gene');
        #$self->{model_variation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variation');
        #$self->{model_variationfeature_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variationfeature');
        #$self->{model_slice_adaptor} = $config->{reg}->get_adaptor($params{model}, 'core', 'slice');
        #$self->{model_variationannotation_adaptor} = $config->{reg}->get_adaptor($params{model}, 'variation', 'variationannotation');

        #if (   ( !defined( $config->{model_translation_adaptor} ) )
        ##    || ( !defined( $config->{model_transcript_adaptor} ) )
        ##    || ( !defined( $config->{model_gene_adaptor} ) )
        ##    || ( !defined( $config->{model_variation_adaptor} ) )
        ##    || ( !defined( $config->{model_variationfeature_adaptor} ) )
        ##    || ( !defined( $config->{model_slice_adaptor} ) )
        ##    || ( !defined( $config->{model_variationannotation_adaptor} ) ) )
        #{
        #    message( $config->{verbose}, $config->{log_file},
        #        "Unable to get adaptors for the species '$config->{model}' specified using the '-model' option.\n"
        #    );
        #}
    }

    $self->{species_short} = \%species_short;

    #ZILA
    # dynamically get species common and latin names
    my ($speciesListCommon_ptr,$speciesListLatin_ptr) = &getEnsemblSpecies();
    $self->{speciesListCommon_ptr} = $speciesListCommon_ptr;
    $self->{speciesListLatin_ptr} = $speciesListLatin_ptr;
    
    # generate distance array from phylogenetic tree
    my $distHash_ptr = &generatePhyloTreeDistArray($config->{species}, $speciesListCommon_ptr, $speciesListLatin_ptr);
    my %distHash = %{$distHash_ptr};
    $self->{distHash_ptr} = $distHash_ptr;
    #print map{"$distHash{$_}\t$_\n"} keys %distHash;
    #ZILA
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
        Alignment_Score_Change
        C_blosum
        Context_Conservation
        Amino_Acids_In_Orthologues
        Orthologue_Species
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Alignment_Score_Change => "the alignment score for the variant amino acid vs. the orthologous amino acids minus the alignment score for the reference amino acid vs. the orthologous amino acids. When there are multiple variant amino acids, the most extreme difference is given. A positive value indicates that the variant amino acid better resembles the orthologues than does the reference amino acid, whereas a negative value indicates that the reference amino acid better resembles the orthologues than does the variant amino acid. The value is scaled to between -1 and 1 (NGS-SNP)",
        C_blosum => "a measure of the conservation of the reference amino acid with the aligned amino acids in orthologous sequences. The alignment score for the reference amino acid vs. the orthologous amino acids is divided by the alignment score that would be obtained if all the orthologous residues matched the reference. Higher values tend to be associated with changes to the amino acid having a greater functional consequence. The formula used is equivalent to the C_blosum formula given in Kowarsch A et al. (2010 PLoS Comput Biol 6(9): e1000923), except that in NGS-SNP scoring matrices other than BLOSUM62 can be used (NGS-SNP)",
        Context_Conservation => "the average percent identity obtained when the region of the reference protein containing the SNP-affected residue is aligned with the orthologous region from other species. The size of the region examined is determined by the -cf option, which specifies how much sequence on either side of the SNP to examine. For example, if '-cf 10' is used, the size of the region is 10 + 1 + 10 = 21 (NGS-SNP)",
        Amino_Acids_In_Orthologues => "the amino acids aligned with the reference amino acid in orthologous sequences (NGS-SNP)",
        Orthologue_Species => "the species from which sequences were obtained to generate the 'Amino_Acid_In_Orthologues', 'Alignment_Score_Change', 'C_blosum', and 'Context_Conservation' values. The order of the species matches the order used to generate the 'Amino_Acid_In_Orthologues' value. Abbreviations - ".join(',',sort(map {"$species_short{$_}=$_"} keys %species_short))." (NGS-SNP)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my $config = $self->{config};
    my %params = %{$self->{params}};
    my %_vax_alignment;
    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $ci = get_consequence_info(@_);
    my %consequences_info;
    if (defined($ci)){
        %consequences_info = %{get_consequence_info(@_)};
    }
    
    # we cache the score on the BaseVariationFeature so we don't have to
    # fetch it multiple times if this variant overlaps multiple Features

    unless (exists $bvf->{_vax_alignment}) {
        if(defined($bvfoa_info{translation})
           && defined($bvfoa_info{ensg})
           && defined($bvfoa_info{ensp})
           && defined($bvfoa_info{amino_acid_reference})
           && defined($bvfoa_info{amino_acid_variant})
           && %consequences_info
           ){
    
            my ($aligned_homolog_residues, $homolog_species) =
                determine_aligned_homolog_residues(
                    $self,
                    $bvfoa_info{changes_protein},
                    $bvfoa_info{ensg},
                    $bvfoa_info{ensp},
                    $bvfoa_info{altered_aa_start},
                    $bvfoa_info{altered_aa_end},
                    $bvfoa_info{amino_acid_reference},
                    $consequences_info{consequences_ranks_so}, #$consequence
                    $params{comparison_species} #$comparison_species
                    );
            my $alignment_score_change =
                determine_alignment_score_change(
                    $self->{matrix},
                    $self->{max_alignment_score_change},
                    $aligned_homolog_residues,
                    $bvfoa_info{amino_acid_reference},
                    $bvfoa_info{amino_acid_variant}
                    );
            $_vax_alignment{Alignment_Score_Change} = $alignment_score_change;
    
            #C_blosom
            my $reference_amino_acid_conservation =
                determine_reference_amino_acid_conservation(
                    $aligned_homolog_residues,
                    $bvfoa_info{amino_acid_reference}
                    );
            my $reference_amino_acid_score =
                determine_reference_amino_acid_score(
                    $self->{matrix},
                    $aligned_homolog_residues,
                    $bvfoa_info{amino_acid_reference}
                    );
            $_vax_alignment{C_blosum} = $reference_amino_acid_score;
    
            #Context_Conservation
            my ($context_alignments, $homolog_species_context) =
                determine_context_alignments(
                    $self,
                    $bvfoa_info{changes_protein},
                    $bvfoa_info{altered_aa_start},
                    $bvfoa_info{protein_sequence},
                    $bvfoa_info{ensg},
                    $bvfoa_info{ensp}
                    );
            my $context_conservation = #context_average_percent_identity
                determine_context_average_percent_identity(
                    $context_alignments
                    );
            $_vax_alignment{Context_Conservation} = $context_conservation;
    
            #Amino_Acids_In_Orthologues & Orthologue_Species
            # order the orthologues according to genetic distance from my species
            if ((defined($aligned_homolog_residues))
                && (scalar( @{$aligned_homolog_residues}) > 0)
                && (defined($homolog_species))
                && (scalar(@{$homolog_species}) > 0))
            {
                ($aligned_homolog_residues, $homolog_species) =
                  &orderOrthologues($aligned_homolog_residues, $homolog_species, $self->{distHash_ptr});
            }
            #Information about protein site conservation
            my $first_long;
            if ((defined($aligned_homolog_residues))
                && (scalar(@{$aligned_homolog_residues}) > 0))
            {
                #separate AAs by semicolons if any are more than one AA
                $first_long = first {length($_)>1} @{$aligned_homolog_residues};
                $_vax_alignment{Amino_Acids_In_Orthologues} =
                    join(defined($first_long) ? ',' : '', @{$aligned_homolog_residues});
             }
            if ((defined($homolog_species))
                && (scalar(@{$homolog_species}) > 0))
            {
                my @species_short;
                for my $sp (@{$homolog_species}){
                    push @species_short, exists($self->{species_short}->{$sp}) ?
                        $self->{species_short}->{$sp} : $sp;
                }
                #my @foo = map {exists($self->{species_short}->{$_}) ?
                #    $self->{species_short}->{$_} : $_} @{$homolog_species};
                $_vax_alignment{Orthologue_Species} =
                    join(',', @{$homolog_species});
                $_vax_alignment{Orthologue_Species} =
                    join(defined($first_long) ? ',' : '', map {exists($self->
                        {species_short}->{$_}) ? $self->{species_short}->{$_} : $_}
                         @{$homolog_species});
            }
            $bvf->{_vax_alignment} = \%_vax_alignment;
        } #if(defined($translation))
    } #unless (exists $bvf->{_conservation_score})

    if (defined $bvf->{_vax_alignment}) {
        $line_hash->{Alignment_Score_Change} = defined($bvf->{_vax_alignment}->{Alignment_Score_Change}) ? $bvf->{_vax_alignment}->{Alignment_Score_Change} : '';
        $line_hash->{C_blosum} = defined($bvf->{_vax_alignment}->{C_blosum}) ? $bvf->{_vax_alignment}->{C_blosum} : '';
        $line_hash->{Context_Conservation} = defined($bvf->{_vax_alignment}->{Context_Conservation}) ? $bvf->{_vax_alignment}->{Context_Conservation} : '';
        $line_hash->{Amino_Acids_In_Orthologues} = defined($bvf->{_vax_alignment}->{Amino_Acids_In_Orthologues}) ? $bvf->{_vax_alignment}->{Amino_Acids_In_Orthologues} : '';
        $line_hash->{Orthologue_Species} = defined($bvf->{_vax_alignment}->{Orthologue_Species}) ? $bvf->{_vax_alignment}->{Orthologue_Species} : '';
    }

    return {};
}

sub get_max_alignment_score_change {
#determines maximum possible alignment score change value
    my $matrix   = shift;
    my @residues = (
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
        'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'
    );
    my $max_score = undef;
    foreach my $row_residue (@residues) {
        my $row_max = undef;
        my $row_min = undef;

        foreach my $column_residue (@residues) {
            my $score = $matrix->get_entry( $row_residue, $column_residue );
            if ( defined($score) ) {
                if ( ( !defined($row_max) ) || ( $score > $row_max ) ) {
                    $row_max = $score;
                }
                if ( ( !defined($row_min) ) || ( $score < $row_min ) ) {
                    $row_min = $score;
                }
            }
        }

        if (   ( !defined($max_score) )
            || ( ( $row_max - $row_min ) > $max_score ) )
        {
            $max_score = $row_max - $row_min;
        }
    }
    return $max_score;
}


sub get_model_orthologous_residues {
    my ($self, $gene_id, $protein_id, $altered_aa_start, $altered_aa_end) = @_;
    if ((!defined($gene_id))
        || (!defined($protein_id))
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
    #direct access to member is deprecated
    #my $member = $self->{ma}
    #    ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    my $member;
    my $gene_member = $self->{gma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    if (defined($gene_member)){
        $member = $gene_member
    }
    else{
        my $seq_member = $self->{sma}
            ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
        if (defined($seq_member)){
            $member = $seq_member
        }
    }

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
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
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
                $protein_id,
                $altered_aa_start
            );
        };
        if ($@) {
            next;
        }
        eval {
            $col_end = $aln->column_from_residue_number(
                $protein_id,
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
        #if (!defined($self->{config}->{verbose})){
            $SIG{'__WARN__'} = sub {return 1};
        #}
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
        if ($seq1->display_id eq $protein_id) {

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

sub determine_aligned_homolog_residues {
    my ($self, $changes_protein, $gene_id, $protein_id, $altered_aa_start, $altered_aa_end, $amino_acid_reference, $consequences_ranks_so, $comparison_species) = @_;
    if ( (!$changes_protein)
        || (!defined($gene_id))
        || (!defined($protein_id))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end))
        || (!defined($amino_acid_reference)))
    {
        return;
    }
    #return values
    my @aligned_homolog_residues;
    my @homolog_species;
    
    if (!($changes_protein)) {
        return;
    }

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    #direct access to member is deprecated
    #my $member = $self->{ma}
    #    ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    my $member;
    my $gene_member = $self->{gma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    if (defined($gene_member)){
        $member = $gene_member
    }
    else{
        my $seq_member = $self->{sma}
            ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
        if (defined($seq_member)){
            $member = $seq_member
        }
    }

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];
    if (defined($self->{params}->{comparison_species})) {

        #obtain orthologues from species of interest
        foreach my $species (keys( %{$self->{params}->{comparison_species}})) {
            push(
                @{$homologies},
                @{$self->{ha}
                    ->fetch_all_by_Member_paired_species($member,
                    $species, ['ENSEMBL_ORTHOLOGUES'])
                }
            );
        }
    }
    else {

        #obtain orthologues from all species
        #push(
        #    @{$homologies},
        #    @{$self->{ha}->fetch_all_by_Member_method_link_type($member,
        #            'ENSEMBL_ORTHOLOGUES')}
        #);
        push(
            @{$homologies},
            @{$self->{ha}->fetch_all_by_Member($member,
                    -METHOD_LINK_TYPE => 'ENSEMBL_ORTHOLOGUES')}
        );
    }

    foreach my $homology (@{$homologies}) {

        #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
        #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef; #Bio::SimpleAlign object

        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if (!(defined($aln))) {
            next;
        }

        #confirm that reference protein is in the alignment
        #Gets an array of Seq objects from the
        #alignment, the contents being those sequences
        #with the given name (there may be more than one)
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
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
                $protein_id,
                $altered_aa_start
            );
        };
        if ($@) {
            next;
        }
        eval {
            $col_end = $aln->column_from_residue_number(
                $protein_id,
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
        #2012-08-18 we don't want them even if verbose is on
        #if (!defined($self->{config}->{verbose})){
            $SIG{'__WARN__'} = sub {return 1};
        #}
        eval { $sub_align = $aln->slice($col, $col_end); };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        if (!(scalar($sub_align->each_seq()) == 2)) {
            next;
        }

        my $seq1 = $sub_align->get_seq_by_pos(1);
        my $seq2 = $sub_align->get_seq_by_pos(2);

        #in all the examples I've examined, $seq1 is the reference sequence
        #counter-example found
        if ($seq1->display_id eq $protein_id || $seq2->display_id eq $protein_id) {
            #swap sequences and taxons if $seq1 is not the reference
            if ($seq1->display_id ne $protein_id){
                my $seqtemp = $seq1;
                $seq1 = $seq2;
                $seq2 = $seqtemp;
                my $taxontemp = $taxon1;
                $taxon1 = $taxon2;
                $taxon2 = $taxontemp;
            }
            if (uc($seq1->seq()) ne
                uc($amino_acid_reference))
            {
#TO: Do we need to deal with this side effect somewhere?
                #if the residue in the alignment is selenocysteine (U), change the reference to U
                #as U in reference seems to be given as * in string returned by $con->pep_allele_string()
                if ((uc($seq1->seq()) eq 'U')
                    && ($amino_acid_reference eq '*'))
                {
                    $amino_acid_reference = 'U';

                    #assume base change in U-coding residue is NON_SYNONYMOUS_CODING
                    #$con->consequence_type seems to be incorrect for sites encoding 'U'
                    
                    #if ($consequence eq 'STOP_LOST') {
                    if (exists $consequences_ranks_so->{stop_lost}){
                        #$consequence = 'NON_SYNONYMOUS_CODING';
                        delete $consequences_ranks_so->{stop_lost};
#TODO: correct rank
                        $consequences_ranks_so->{missense_variant} = 0;
                    }
                }
                else {
#TODO: something other than ignore things like this
#VCF: 17	39340795	splice_acceptor_variant	ACGGCAGCAGCTGGACATACCACAGCTGGGGTGGCAGGTGGTCTGACAGCAGAGTGGG	A	.	.	.
#protein_id: ENSP00000381489
#seq1->display_id: ENSP00000381489
#amino_acid_reference: RPLCCQTTCHPSCGMSSCCR
#seq1->seq: RPLCCQTTCHP---SCGMSSCCR
                    #die("The alignment '"
                    #        . $seq1->seq()
                    #        . "' and reference '$amino_acid_reference' amino acids do not match for sequence $protein_id."
                    #);
                    next;
                }
            }
            push(
                @aligned_homolog_residues,
                $seq2->seq()
            );
            push(
                @homolog_species,
                space_to_underscore(get_taxon_name($taxon2))
            );
        }
        else {
            die("Unexpected sequence ordering in alignment.");
        }
    }
    return (\@aligned_homolog_residues, \@homolog_species)
}

sub determine_alignment_score_change {
    my ($matrix, $max_alignment_score_change, $aligned_homolog_residues, $amino_acid_reference, $amino_acid_variant) = @_;

    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_score = get_average_similarity_score(
        $matrix,
        $amino_acid_reference,
        $aligned_homolog_residues
    );
    my $variant_score = get_average_similarity_score(
        $matrix,
        $amino_acid_variant,
        $aligned_homolog_residues
    );

    my $alignment_score_change = sprintf("%.3f",
        ($variant_score - $reference_score) / $max_alignment_score_change);
    return $alignment_score_change;
}

sub get_average_similarity_score {
#compares $residue to each amino acid in $array_of_aligned
#and determines average similarity score
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum   = 0;
    my $score_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
            $score_count++;
        }
    }

    my $score = undef;
    if ( $score_count == 0 ) {
        $score = undef;
    }
    else {
        $score = $score_sum / $score_count;
        $score = sprintf( "%.2f", $score );
    }
    return $score;
}

sub determine_reference_amino_acid_conservation {
    my ($aligned_homolog_residues, $amino_acid_reference) = @_;

    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_amino_acid_conservation = get_percent_identity_column(
        $amino_acid_reference,
        $aligned_homolog_residues
    );

    if (defined($reference_amino_acid_conservation)) {
        my $reference_amino_acid_conservation
            = sprintf("%.1f", $reference_amino_acid_conservation);
    }
    return $reference_amino_acid_conservation;
}

sub get_percent_identity_column {
#compares $residue to each amino acid in $array_of_aligned
#and determines percentage of residues in $array_of_aligned
#that are identical.
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give 50
#G vs -G will give 100
#
#Is equivalent to "Cident" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $match_count   = 0;
    my $checked_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }
        $checked_count++;
        if ( uc($residue) eq uc($aligned) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;
}

sub determine_reference_amino_acid_score {
    my ($matrix, $aligned_homolog_residues, $amino_acid_reference) = @_;
 
    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_amino_acid_score = get_score_column(
        $matrix,
        $amino_acid_reference,
        $aligned_homolog_residues
    );

    if (defined($reference_amino_acid_score)) {
        $reference_amino_acid_score
            = sprintf("%.1f", $reference_amino_acid_score);
    }
    return $reference_amino_acid_score;
}

sub get_score_column {
#compares $residue to each amino acid in $array_of_aligned
#and determines sum of scores using scoring matrix.
#this value is then divided by the score obtained if all
#residues in $array_of_aligned were to match $residue
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give score(G,G) + score(G,A) / 2 * score(G,G)
#G vs -G will give score(G,G) / score(G,G)
#
#Is equivalent to "Cblosum" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum     = 0;
    my $max_score_sum = 0;

    my $max_score = $matrix->get_entry( uc($residue), uc($residue) );
    if ( !defined($max_score) ) {
        return undef;
    }

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }

        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
        }

        $max_score_sum = $max_score_sum + $max_score;

    }
    if ( $max_score_sum != 0 ) {
        my $score = sprintf( "%.1f", $score_sum / $max_score_sum );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "column score: $score\n";

        return $score;
    }
    return undef;

}

sub determine_context_alignments {
    my ($self, $changes_protein, $altered_aa_start, $protein_sequence, $gene_id, $protein_id) = @_;

    if (! $changes_protein) {
        return;
    }
    
    #return values
    my (@context_alignments, @homolog_species_context);
    
    my $context_start = $altered_aa_start
        - $self->{params}->{flanking_for_context_conservation};
    my $context_end = $altered_aa_start
        + $self->{params}->{flanking_for_context_conservation};

    if ( $context_start < 1 ) {
        $context_start = 1;
    }
    if ($context_end > length($protein_sequence)) {
        $context_end = length($protein_sequence);
    }

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    #direct access to member is deprecated
    #my $member = $self->{ma}
    #    ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    my $member;
    my $gene_member = $self->{gma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
    if (defined($gene_member)){
        $member = $gene_member
    }
    else{
        my $seq_member = $self->{sma}
            ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);
        if (defined($seq_member)){
            $member = $seq_member
        }
    }

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];

    if (defined($self->{params}->{comparison_species})) {

        #obtain orthologues from species of interest
        foreach my $species (keys(%{$self->{params}->{comparison_species}})) {
            push(
                @{$homologies},
                @{$self->{ha}->fetch_all_by_Member_paired_species($member,
                    $species, ['ENSEMBL_ORTHOLOGUES'])}
            );
        }
    }
    else {

        #obtain orthologues from all species
        #push(
        #    @{$homologies},
        #    @{$self->{ha}->fetch_all_by_Member_method_link_type($member,
        #            'ENSEMBL_ORTHOLOGUES')}
        #);
        push(
            @{$homologies},
            @{$self->{ha}->fetch_all_by_Member($member,
                    -METHOD_LINK_TYPE => 'ENSEMBL_ORTHOLOGUES')}
        );
    }

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

        if (!(defined($aln))) {
            next;
        }

        #confirm that reference protein is in the alignment
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
        {
            next;
        }

        #get an alignment containing the SNP-affected column and flanking sequence
        my $sub_align_context = undef;

        eval {
            my $col_context_start
                = $aln->column_from_residue_number(
                $protein_id,
                $context_start);

            my $col_context_end = $aln->column_from_residue_number(
                $protein_id, $context_end);

            $sub_align_context
                = $aln->slice($col_context_start, $col_context_end);
        };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        my $context_string = undef;
        if ((defined($sub_align_context))
            && (scalar($sub_align_context->each_seq()) == 2))
        {
            my $seq1_context = $sub_align_context->get_seq_by_pos(1);
            my $seq2_context = $sub_align_context->get_seq_by_pos(2);

            my %seq_hash = (reference => undef, homolog => undef);

           #in all the examples I've examined, $seq1 is the reference sequence
            #counter-example found
            if ($seq1_context->display_id eq $protein_id || $seq2_context->display_id eq $protein_id) {
                #swap sequences and taxons if $seq1 is not the reference
                if ($seq1_context->display_id ne $protein_id){
                    my $seq_contexttemp = $seq1_context;
                    $seq1_context = $seq2_context;
                    $seq2_context = $seq_contexttemp;
                    my $taxontemp = $taxon1;
                    $taxon1 = $taxon2;
                    $taxon2 = $taxontemp;
                }

                $seq_hash{reference} = $seq1_context->seq();
                $seq_hash{homolog}   = $seq2_context->seq();

                push( @context_alignments, \%seq_hash );
                push(
                    @homolog_species_context,
                    space_to_underscore( get_taxon_name($taxon2) )
                );
            }
            else {
                die("Unexpected sequence ordering in alignment.");
            }
        }
    }
    return (\@context_alignments, \@homolog_species_context);
}

sub determine_context_average_percent_identity {
    my ($context_alignments) = @_;

    if ((!defined($context_alignments))
        || (scalar(@{$context_alignments}) == 0))
    {
        return;
    }
    
    #return value
    my $context_conservation;

    my $identity_sum   = 0;
    my $identity_count = 0;
    foreach my $subsequence_alignment (@{$context_alignments})
    {
        my $percent_identity = get_percent_identity_between_two_seqs(
            $subsequence_alignment->{reference},
            $subsequence_alignment->{homolog}
        );
        if (defined($percent_identity)) {
            $identity_sum = $identity_sum + $percent_identity;
            $identity_count++;
        }
    }

    if ( $identity_count > 0 ) {
        $context_conservation
            = sprintf( "%.1f", ( $identity_sum / $identity_count ) );
    }
    return $context_conservation;
}
sub get_percent_identity_between_two_seqs {
#determines percent identity between two sequences
#by looking for character matches.
#
#Aligning gaps are ignored
#
#--G
#--G will give 100 match
#
#-GA
#-GG will give 50 match

    my $seq1 = shift;
    my $seq2 = shift;

    my @seq1 = split( //, $seq1 );
    my @seq2 = split( //, $seq2 );

    my $match_count   = 0;
    my $checked_count = 0;
    for ( my $i = 0; $i < scalar(@seq1); $i++ ) {
        if ( ( is_gap( $seq1[$i] ) ) || ( is_gap( $seq2[$i] ) ) ) {
            next;
        }

        $checked_count++;
        if ( uc( $seq1[$i] ) eq uc( $seq2[$i] ) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "seq1: $seq1\n";
        #   print "seq2: $seq2\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;
}

#ZILA
sub generatePhyloTreeDistArray {
    my $species = shift;
    my $speciesListCommon_ptr = shift;
    my @speciesListCommon = @{$speciesListCommon_ptr};
    my $speciesListLatin_ptr = shift;
    my @speciesListLatin = @{$speciesListLatin_ptr};
    #return will be ptr to hash of species:distance
    my %distanceHash;
  
	# get index of our species, default to Homo_sapiens
 	my $speciesIndex = first { lc($speciesListLatin[$_]) eq lc($species) } 0..$#speciesListLatin;
    if (! defined($speciesIndex)){
        $species = 'Homo_sapiens';
        $speciesIndex = first { lc($speciesListLatin[$_]) eq lc($species) } 0..$#speciesListLatin;
    }
	# check defined species exists
	die("ERROR: Could not find species \"", $species, "\" in Ensembl species list\n") unless defined $speciesIndex;

	# default phylogenetic tree; branch lengths for Ensembl species phylogenetic tree obtained from http://tinyurl.com/ensembltree on 2011-08-20
    #my $default_phyloTree = lc("(((((((((((((((((((((((Homo_sapiens:0.0067,Pan_troglodytes:0.006667):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,(Macaca_mulatta:0.007853,?Papio_hamadryas:0.007637):0.029618):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Cavia_porcellus:0.225629):0.01015,Spermophilus_tridecemlineatus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((Vicugna_pacos:0.107275,(Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,?Ovis_aries:0.061796):0.061796):0.025153):0.0201675,Sus_scrofa:0.079):0.0201675,((Equus_caballus:0.109397,(Felis_catus:0.098612,(Ailuropoda_melanoleuca:0.051229,Canis_familiaris:0.051229):0.051229):0.049845):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(Monodelphis_domestica:0.125686,Macropus_eugenii:0.122008):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,Taeniopygia_guttata:0.171542):0.199223,Anolis_carolinensis:0.489241):0.105143):0.172371,Xenopus_tropicalis:0.855573):0.311354,(((Tetraodon_nigroviridis:0.224159,Takifugu_rubripes:0.203847):0.195181,(Gasterosteus_aculeatus:0.316413,Oryzias_latipes:0.48197):0.05915):0.32564,Danio_rerio:0.730752):0.147949):0.526688,?Petromyzon_marinus:0.526688),(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(?Apis_mellifera:0.9,(((?Aedes_aegypti:0.25,?Culex_quinquefasciatus:0.25):0.25,?Anopheles_gambiae:0.5)Culicinae:0.2,Drosophila_melanogaster:0.8)Diptera:0.1)Endopterygota:0.7)Coelomata:0.1,Caenorhabditis_elegans:1.7)Bilateria:0.3,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group:0.3)");
	# default phylogenetic tree; branch lengths for Ensembl species phylogenetic tree obtained from http://tinyurl.com/ensembltree on 2012-01-18
    my $default_phyloTree = lc("((((((((((((((((((((((((Homo_sapiens:0.0067,Pan_troglodytes:0.006667):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,(Macaca_mulatta:0.007853,?Papio_hamadryas:0.007637):0.029618):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Cavia_porcellus:0.225629):0.01015,Spermophilus_tridecemlineatus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((Vicugna_pacos:0.107275,(Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,?Ovis_aries:0.061796):0.061796):0.025153):0.0201675,Sus_scrofa:0.079):0.0201675,((Equus_caballus:0.109397,(Felis_catus:0.098612,(Ailuropoda_melanoleuca:0.051229,Canis_familiaris:0.051229):0.051229):0.049845):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(Monodelphis_domestica:0.125686,(Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,Taeniopygia_guttata:0.171542):0.199223,Anolis_carolinensis:0.489241):0.105143):0.172371,Xenopus_tropicalis:0.855573):0.155677,Latimeria_chalumnae:0.155677):0.155677,((((Tetraodon_nigroviridis:0.224159,Takifugu_rubripes:0.203847):0.195181,(Gasterosteus_aculeatus:0.316413,Oryzias_latipes:0.48197):0.05915):0.16282,Gadus_morhua:0.16282):0.16282,Danio_rerio:0.730752):0.147949):0.526688,Petromyzon_marinus:0.526688):0.526688,(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(?Apis_mellifera:0.9,(((?Aedes_aegypti:0.25,?Culex_quinquefasciatus:0.25):0.25,?Anopheles_gambiae:0.5)Culicinae:0.2,Drosophila_melanogaster:0.8)Diptera:0.1)Endopterygota:0.7)Coelomata:0.1,Caenorhabditis_elegans:1.7)Bilateria:0.3,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group:0.3);");
    # dynamically get tree with branch lengths
    my $url = 'http://tinyurl.com/ensembltree';
    my $content = get $url;
    my $phyloTree = defined $content ? lc($content) : $default_phyloTree;

    my @distanceToSpecies = ((-1) x scalar(@speciesListLatin));
    for (my $i = 0; $i < scalar(@speciesListLatin); $i++){
        $distanceHash{$speciesListLatin[$i]} = -1;
    }
	my $specieslc;
	my $currspecieslc;

	# trim branches of tree that aren't in our species list
	my @trimmedtree = split(/([\(\):,])/, $phyloTree);
	for(my $i = 0; $i < scalar(@trimmedtree); $i++) {
		my $item = $trimmedtree[$i];
		$item =~ s/\?//;
		my $itemIndex = first { lc($speciesListLatin[$_]) eq $item } 0..$#speciesListLatin;
		if (($item =~ m/[a-z_]+/)  && !(defined $itemIndex)) {
			#print "Item $item not found. Removing from tree.\n";
			my @replacement = ("","","");
			splice(@trimmedtree, $i, 3, @replacement);			
		}
	}
	$phyloTree = join("", @trimmedtree);
	#print $phyloTree."\n";

	# calculate distances between our species and every other species
	for(my $i = 0; $i < scalar(@speciesListLatin); $i++) {
		# set distance to our own species to zero
		if ( $i == $speciesIndex) {
			$distanceToSpecies[$i] = 0;
            $distanceHash{$speciesListLatin[$i]} = 0;
		}
		else {
			# does our species and the current species exist in our phylogenetic tree?
			$specieslc = lc($species);
			$currspecieslc = lc($speciesListLatin[$i]);
			if (($phyloTree =~ m/$specieslc/) && ($phyloTree =~ m/$currspecieslc/)) {
				#print "Found ".$specieslc." and ".$currspecieslc."\n";

				# make a copy of our tree
				# e.g. ((((T:0.7,F:0.3):0.1,(G:0.4,O:0.5):0.6):0.2;D:0.8):0.9;A:1.1)
				my $tree = $phyloTree;

				# trim all leaves of the tree that aren't our species or the current species
				# e.g. if our species are F and D --> ((((,F:0.3):0.1,(,):0.6):0.2;D:0.8):0.9,)
				for(my $j = 0; $j < scalar(@speciesListLatin); $j++) {
					my $jspecies = lc($speciesListLatin[$j]);
					if (($jspecies ne $specieslc) && ($jspecies ne $currspecieslc)) {
						$tree =~ s/$jspecies:\d.[\d]+//g;
					}
				}

				# trim all branches of the tree that aren't on the path between the two species of interest
				# e.g. if our species are F and D --> (((,F:0.3):0.1,):0.2;D:0.8)
				while ($tree =~ m/\(,\):\d.[\d]+/) {
					$tree =~ s/\(,\):\d.[\d]+//g;
				}
				while ($tree =~ m/\(,\)/) {
					$tree =~ s/\(,\)//g;
				}
				# find set of parentheses that encloses our species and discard everything outside
				# code adapted from http://www.perlmonks.org/?node_id=660316
				my @queue = ( $tree );			
				my $regex = qr/
					(			# start of bracket 1
					\(			# match an opening parentheses
						(?:               
						[^\(\)]++	# one or more non parentheses, non backtracking
							|                  
						(?1)		# recurse to parentheses 1
						)*                 
					\)			# match a closing parentheses
					)			# end of parentheses 1
					/x;

				$" = "\n\t";

				my @potentials;
				while( @queue )
				{
					my $string = shift @queue;
					my @groups = $string =~ m/$regex/g;
					#print "Found:\n\t@groups\n\n" if @groups;
					if (($string =~ m/$specieslc/) && ($string =~ m/$currspecieslc/)) {
						push @potentials, $string;
					}
					unshift @queue, map { s/^\(//; s/\)$//; $_ } @groups;
				}
				#print "Potentials: ".join("\n", @potentials)."\n";
				my $shortest = $potentials[0];
				foreach my $potential (@potentials) {
					$shortest = $potential if length($potential) < length($shortest);
				    }
				#print "Shortest: (".$shortest.")\n";

				# add together all remaining numbers in our string
				# e.g. 0.3 + 0.1 + 0.2 + 0.8 --> 1.4
				my @splittree = split(/([\(\):,])/, $shortest);
				my $sum = 0;
				foreach my $item (@splittree) {
					if ($item =~ /\d.[\d]+/) {	# if it's a number
						$sum += $item;	# branch lengths to scale
						#$sum += 0.1;	# fixed branch lengths
					}
				}
				#print "Distance between ".$speciesListCommon->[$speciesIndex]." and ".$speciesListCommon->[$i]." is: ".$sum."\n";
				$distanceToSpecies[$i] = $sum;
                $distanceHash{$speciesListLatin[$i]} = $sum;
			}
		}
	}
    for (my $i = 0; $i < scalar(@speciesListLatin); $i++){
      $distanceHash{$speciesListLatin[$i]} = $distanceToSpecies[$i];
    }
	return \%distanceHash;
}

sub orderOrthologues {
#rearrange orthologues by evolutionary distance of species
    my $AAsInOrthologues_ptr = shift;
    my @AAsInOrthologues = @{$AAsInOrthologues_ptr};
    my $orthologues_ptr = shift;
    my @orthologues = @{$orthologues_ptr};
    #species:distance
    my $distHash_ptr = shift;
    my %distHash = %{$distHash_ptr};
    
    #uniqueify species/AAs
    my %distinctOrthologues;
    for (my $i=0;$i<scalar(@orthologues);$i++){
      $distinctOrthologues{$orthologues[$i]}->{$AAsInOrthologues[$i]}++;
    }
    
    #sort species by evolutionary distance
    my %positives = map {$_=>$distHash{$_}} grep {$distHash{$_}>=0} keys %distHash;
    my @sortedSpecies = sort {$positives{$a} <=> $positives{$b}} keys %positives;
    #append any species with no or negative distance
    while(my($species,$aas) = each(%distinctOrthologues)){
      push @sortedSpecies, $species unless exists($positives{$species});
    }
    
    #push species and AAs to new arrays ordered by evolutionary distance
    my @sortedAAsInOrthologues = ();
    my @sortedOrthologues = ();
    for my $thisSpecies(@sortedSpecies){
        while (my($aa, $count) = each(%{$distinctOrthologues{$thisSpecies}})){
            push @sortedOrthologues, $thisSpecies;
            push @sortedAAsInOrthologues, $aa;
        }
    }
    return \(@sortedAAsInOrthologues, @sortedOrthologues);
}

sub getEnsemblSpecies {
    # assign arrays of species common and latin names
    #defaults downloaded from http://www.ensembl.org/info/about/species.html on 2011-08-20
    my @speciesListCommon = qw(Alpaca Gorilla Pig Anole_Lizard Guinea_Pig Pika Armadillo Hedgehog Platypus Horse Rabbit Bushbaby Human Rat C.elegans Hyrax S.cerevisiae C.intestinalis Kangaroo_rat C.savignyi Shrew Cat Lesser_hedgehog_tenrec Sloth Chicken Macaque Squirrel Chimpanzee Marmoset Stickleback Cow Medaka Tarsier Dog Megabat Tetraodon Dolphin Microbat Tree_Shrew Mouse Turkey Elephant Mouse_Lemur Wallaby Fruitfly Opossum X.tropicalis Fugu Orangutan Zebra_Finch Gibbon Panda Zebrafish);
    my @speciesListLatin = qw(Vicugna_pacos Gorilla_gorilla Sus_scrofa Anolis_carolinensis Cavia_porcellus Ochotona_princeps Dasypus_novemcinctus Erinaceus_europaeus Ornithorhynchus_anatinus Equus_caballus Oryctolagus_cuniculus Otolemur_garnettii Homo_sapiens Rattus_norvegicus Caenorhabditis_elegans Procavia_capensis Saccharomyces_cerevisiae Ciona_intestinalis Dipodomys_ordii Ciona_savignyi Sorex_araneus Felis_catus Echinops_telfairi Choloepus_hoffmanni Gallus_gallus Macaca_mulatta Spermophilus_tridecemlineatus Pan_troglodytes Callithrix_jacchus Gasterosteus_aculeatus Bos_taurus Oryzias_latipes Tarsius_syrichta Canis_familiaris Pteropus_vampyrus Tetraodon_nigroviridis Tursiops_truncatus Myotis_lucifugus Tupaia_belangeri Mus_musculus Meleagris_gallopavo Loxodonta_africana Microcebus_murinus Macropus_eugenii Drosophila_melanogaster Monodelphis_domestica Xenopus_tropicalis Takifugu_rubripes Pongo_abelii Taeniopygia_guttata Nomascus_leucogenys Ailuropoda_melanoleuca Danio_rerio);
    my $html = get("http://uswest.ensembl.org/info/about/species.html");
    if (defined $html){
      $html =~ m{Ensembl [Ss]pecies[</a-z0-9>\s]*<table(.*?)</table>}s;
      if (defined $1){
          @speciesListCommon = ();
          @speciesListLatin = ();
          my $speciestable = $1;
          while ($speciestable =~ m{<a href=.*?>(.*?)</a>[)]?<br />(.*?)</td>}g) {
            if (!($1 =~ m{preview })) { 	# skip species that are only preview
              my $commonname = $1;
              $commonname =~ tr[ ][_];
              my $latinname;
              my $temp = $2;
              if ($temp =~ m{<i>(.*?)</i>}) {
                $latinname = $1;
                $latinname =~ tr[ ][_];
                if ($commonname eq $latinname) {
                  $commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
                  $commonname = substr($1, 0, 1).".".$2;
                }
              }
              else {
                $latinname = $commonname;
                $commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
                $commonname = substr($1, 0, 1).".".$2;
              }
              push(@speciesListCommon, $commonname);
              push(@speciesListLatin, $latinname);
            }		
          }
      }
    }
    #print '"'.join('", "',@speciesListCommon).'"'."\n";
    #print '"'.join('", "',@speciesListLatin).'"'."\n";
    return \(@speciesListCommon,@speciesListLatin);
}
#ZILA


1;

