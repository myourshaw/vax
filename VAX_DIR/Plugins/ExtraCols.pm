=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 ExtraCols

=head1 SYNOPSIS

 mv ExtraCols.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf [--plugin Condel ...] --plugin ExtraCols
 
 Plugins, such as Condel, must be installed and selectd to run for this plugin to use their data.
 
 In the config file or command line, the ExtraCols plugin must be specified after
 plugins, such as Condel, that generate data consumed by ExtraCols

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the VEP's Extra data to additional output columns.
 By default, it adds columns for the built-in VEP extra headers
 and for the downloadable Carol, Condel, Conservation, dbNSFP, FATHMM, and Grantham plugins.
 
 It expands the PolyPhen and SIFT values to prediction and score:
 PolyPhen_prediction, PolyPhen_score, SIFT_prediction, SIFT_score
 and similarly expands the Carol, Condel, and FATHMM plugin values.

 Aditional columns may be specified as parameters, to accomodate other plugin values
 Without needing to modify this code.
 
 A default EXTRA KEY can be suppressed with a negated parameter, e.g., -MAF.
 
 All parameters are case-sensitive.

=head1 PARAMETERS

    comma-separated list of EXTRAS keys (case-sensitive)

=cut

package ExtraCols;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my %config = %{$self->{config}};
    my @plugins = @{$config{plugins}};
    my @params = @{$self->{params}};
    
    #parameters that add new or remove default columns
    my (%add_cols, %subtract_cols);
    while (my ($index,$param) = each @params) {
        if (substr($param,0,1) ne '-') {
            $add_cols{$param} = 1;
        }
        else {
            $subtract_cols{substr($param,1)} = 1;
        }
        
    }
        
    #keys that will be added to EXTRAS
    my %extra_keys;
    
    #VEP built-in EXTRA columns
    #from variant_effect_predictor.pl
    # define headers that would normally go in the extra field
    # keyed on the config parameter used to turn it on
    my %extra_headers = (
        allele_number   => ['ALLELE_NUM'],
        biotype         => ['BIOTYPE'],
        canonical       => ['CANONICAL'],
        ccds            => ['CCDS'],
        cell_type       => ['CELL_TYPE'],
        check_existing  => ['CLIN_SIG'],
        check_frequency => ['FREQS'],
        check_svs       => ['SV'],
        domains         => ['DOMAINS'],
        gmaf            => ['GMAF'],
        hgvs            => ['HGVSc','HGVSp'],
        individual      => ['IND','ZYG'],
        maf_1kg         => ['AFR_MAF','AMR_MAF','ASN_MAF','EUR_MAF'],
        maf_esp         => ['AA_MAF','EA_MAF'],
        numbers         => ['EXON','INTRON'],
        polyphen        => ['PolyPhen'],
        protein         => ['ENSP'],
        pubmed          => ['PUBMED'],
        regulatory      => ['MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE'],
        sift            => ['SIFT'],
        symbol          => ['SYMBOL','SYMBOL_SOURCE'],
        user            => ['DISTANCE'],
        xref_refseq     => ['RefSeq'],
    );
    
    #add keys for each built-in extra that has been activated with a config option
    while (my ($option, $headers) = each %extra_headers) {
        if(exists($config{$option})) {
            while (my ($index,$col) = each @{$headers}) {
                $extra_keys{$col} = 1 unless exists($subtract_cols{$col});
            }
        }
    }

    #plugin EXTRA keys
    #key=plugin name; value=EXTRA key
    #
    my %plugin_headers = (
        Carol => ['CAROL'],
        Condel => ['Condel'],
        Conservation => ['Conservation'],
        FATHMM => ['FATHMM'],
        GO => ['GO'], #this may need to be put in the GO_VEP column to avoid conflict with the UniProt GO
        Grantham => ['Grantham'],
        dbNSFP => [],
    );
    
    #add keys for each plugin in %plugin_headers
    foreach my $plugin (@plugins) {
        my $plugin_type = Scalar::Util::blessed($plugin);
        if (exists($plugin_headers{$plugin_type})) {
            if (scalar(@{$plugin_headers{$plugin_type}} == 0)) {
                if ($plugin_type eq 'dbNSFP') {
                    @{$plugin_headers{$plugin_type}} = keys %{$plugin->{cols}};
                }
            }
            map {$extra_keys{$_} = 1 unless exists($subtract_cols{$_})} @{$plugin_headers{$plugin_type}};
        }
        
    }
    
    #add user specified columns to hash
    map { $extra_keys{$_} = 1 } keys %add_cols;
    
    $self->{extra_keys} = \%extra_keys;
    
    #list of keys in EXTRAS that will be expanded to prediction/score
    my %prediction_score_keys = (PolyPhen => 1, SIFT => 1, Condel => 1, CAROL => 1, FATHMM => 1,);
    $self->{prediction_score_keys} = \%prediction_score_keys;
    
    #extra column names that need to be renamed
    my %column_translation = (
        FOO => "BAR",
    );
    $self->{column_translation} = \%column_translation;

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
    my $class = shift;
    my %extra_keys = %{$class->{extra_keys}};
    my %prediction_score_keys = %{$class->{prediction_score_keys}};
    my %column_translation = %{$class->{column_translation}};
    
    #add extra columns for prediction/score
    my %extra_columns = %extra_keys;
    my @types = ('prediction', 'score');
    while (my ($key,$value) = each %prediction_score_keys) {
        if (exists($extra_columns{$key})) {
            #uncomment next line to exclude combined 'prediction (score)' column
            #delete($extra_columns{$col});
            while (my ($index,$col_type) = each @types) {
                $extra_columns{$key.'_'.$col_type} = 1;
            }
        }
    }
    $class->{extra_columns} = \%extra_columns;
    
    my @new_output_cols = map {exists($column_translation{$_}) ? $column_translation{$_} : $_} keys %extra_columns;
    
    #define ordering of new extra output columns
    my @Extra_cols_order= qw(
        SYMBOL
        SYMBOL_SOURCE
        CANONICAL
        ENSP
        CCDS
        HGVSc
        HGVSp
        Conservation
        PolyPhen
        PolyPhen_prediction
        PolyPhen_score
        SIFT
        SIFT_prediction
        SIFT_score
        Condel
        Condel_prediction
        Condel_score
        CAROL
        CAROL_prediction
        CAROL_score
        FATHMM
        FATHMM_prediction
        FATHMM_score
        Grantham
        EXON
        INTRON
        DOMAINS
        MOTIF_NAME
        MOTIF_POS
        HIGH_INF_POS
        MOTIF_SCORE_CHANGE
        CELL_TYPE
        IND
        ZYG
        RefSeq
        SV
        FREQS
        GMAF
        AFR_MAF
        AMR_MAF
        ASN_MAF
        EUR_MAF
        AA_MAF
        EA_MAF
        PUBMED
        DISTANCE
        CLIN_SIG
        BIOTYPE
        ALLELE_NUM
        GO
        Interpro_domain
        SLR_test_statistic
        SIFT_score_converted
        LRT_score
        LRT_score_converted
        LRT_pred
        MutationTaster_score
        MutationTaster_score_converted
        MutationTaster_pred
        MutationAssessor_score
        MutationAssessor_score_converted
        MutationAssessor_pred
        FATHMM_score_converted
        FATHMM_pred
        GERP++_NR
        GERP++_RS
        phyloP
        29way_pi
        29way_logOdds
        LRT_Omega
        1000Gp1_AC
        1000Gp1_AF
        1000Gp1_AFR_AC
        1000Gp1_AFR_AF
        1000Gp1_EUR_AC
        1000Gp1_EUR_AF
        1000Gp1_AMR_AC
        1000Gp1_AMR_AF
        1000Gp1_ASN_AC
        1000Gp1_ASN_AF
        ESP6500_AA_AF
        ESP6500_EA_AF
    );
    
    my %extra_cols_ordinals;
    while (my ($index,$col) = each @Extra_cols_order) {
        $extra_cols_ordinals{$col} = $index;
    }
    my %new_output_cols_ordinals;
    my $i = scalar(@Extra_cols_order);
     while (my($index,$col) = each @new_output_cols) {
        $new_output_cols_ordinals{$col} = exists($extra_cols_ordinals{$col}) ? $extra_cols_ordinals{$col} : $i++;
    }
    my @sorted_new_output_cols = map {$_} sort { $new_output_cols_ordinals{$a} <=> $new_output_cols_ordinals{$b} or $a cmp $b } keys %new_output_cols_ordinals;
    
    @OUTPUT_COLS = (@OUTPUT_COLS, @sorted_new_output_cols);

    #add header definitions for any expanded prediction/score columns
    my %header_items;
    while (my ($key,$value) = each %prediction_score_keys) {
        while (my ($index,$col_type) = each @types) {
            if (exists($new_output_cols_ordinals{$key.'_'.$col_type})) {
                $header_items{$key.'_'.$col_type} = "$key $col_type";
            }
        }
    }

    return \%header_items;
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my %extra_columns = %{$self->{extra_columns}};
    my %prediction_score_keys = %{$self->{prediction_score_keys}};
    my %column_translation = %{$self->{column_translation}};

    #add values from EXTRAS to individual columns if there is a column for them
    my %line_hash_extra = %{$line_hash->{Extra}};
    while (my($key,$value) = each %line_hash_extra ) {
        #expand prediction/score columns if this is an EXTRA key that should be expanded
        #and add to output line
        if (exists($prediction_score_keys{$key})) {
            #if ($value =~ /([^(]+)\(([0-9.-]+)/){
            if ($value =~ /([^(]+)\(([^)]+)/){
                if (exists($extra_columns{$key.'_prediction'})) {
                    $line_hash->{$key.'_prediction'} = $1;
                }
                if (exists($extra_columns{$key.'_score'})) {
                    $line_hash->{$key.'_score'} = $2;
                }
            }
        }
        #add EXTRA key (including unexpanded keys of those that were expanded) to output
        if (exists($extra_columns{$key})) {
            if (exists($column_translation{$key})) {
                $line_hash->{$column_translation{$key}} = $value;
            }
            else{
                $line_hash->{$key} = $value;
            }
        }
    }
    
    return {};
}


1;

