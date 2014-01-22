=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 GeneIDs

=head1 SYNOPSIS

 mv GeneIDs.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin GeneIDs

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 strand, ENSG, Gene_Description, RefSeq_summary, Entrez_Gene_Name, UniProtKB_AC, UniProt_ID, MT.
 
 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
 
 Requires that the VAX.pm module be in the Plugins directory

 References:
    (1) Pagliarini DJ, Calvo SE, Chang B et al.
        A mitochondrial protein compendium elucidates complex I disease biology
        Cell 2008;134:112-123
        
    (2) Pruitt KD, Tatusova T, Brown GR et al.
        NCBI Reference Sequences (RefSeq): current status, new features and genome annotation policy
        Nucleic Acids Res 2012;40:D130-135

=cut

package GeneIDs;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info get_unique);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my $goa = $self->{config}->{reg}->get_adaptor('Multi', 'Ontology', 'GOTerm');
    $self->{goa} = $goa;

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
        strand
        ENSG
        Gene_Description
        RefSeq_summary
        Entrez_Gene_Name
        UniProt_ID
        MT
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        strand => "Strand of transcript (Ensembl)",
        ENSG => "Gene stable ID (Ensembl)",
        Gene_Description => "A short description of the gene (Ensembl)",
        RefSeq_summary => "Gene summary (RefSeq)",
        Entrez_Gene_Name => "The Entrez Gene name of the relevant gene (Ensembl)",
        UniProt_ID => "The UniProt ID of the relevant protein (Ensembl)",
        MT => "1 if the gene is annotated as a mitochondrial gene by MitoCarta (MitoCarta)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    if (defined ($bvfoa_info{transcript})){
        my $transcript = $bvfoa_info{transcript};
        $line_hash->{strand} = $transcript->strand;
    }
    
    my $gene = $bvfoa_info{gene};
    if (defined $gene){
        $line_hash->{ENSG} = $bvfoa_info{ensg};

        if (defined $gene->{description}) {
            $line_hash->{Gene_Description} = $gene->{description}
        }
        
        my $xrefs = $gene->get_all_DBLinks(); 
        foreach my $xref (@{$xrefs}) {
            if ($xref->dbname() eq 'EntrezGene') {
                $line_hash->{Entrez_Gene_Name} = $xref->display_id();
#                $line_hash->{Entrez_Gene_ID} = $xref->primary_id();
                last;
            }
        }
        
        $xrefs = defined($bvfoa_info{transcript}) ?
            $bvfoa_info{transcript}->get_all_DBLinks() :
            $gene->get_all_DBLinks(); 
        my @uniprotkb_ac;
        my @uniprot_id;
        my @omim_id;
        my @hgmd_id;
        my @go_names = ();
        foreach my $xref (@{$xrefs}) {
            if ($xref->dbname() eq 'Uniprot/SWISSPROT') {
                push(
                    @uniprotkb_ac,
                    $xref->primary_id()
                );
                push(
                    @uniprot_id,
                    $xref->display_id()
                );
            }
            if ($xref->dbname() =~/^MIM/) {
                push(
                    @omim_id,
                    $xref->display_id()
                );
            }
            #if ($xref->dbname() =~ /HGMD/) {
            #    push(
            #        @hgmd_id,
            #        $xref->display_id()
            #    );
            #}
            #elsif ($xref->dbname() eq 'goslim_goa') {
            #    my $go_id = $xref->display_id();
            #    my $go_term;
            #    my $go_name;
            #    my $go_definition;
            #    if (defined($self->{goa})) {
            #        $go_term = $self->{goa}->fetch_by_accession($go_id);
            #        $go_name = $go_term->name();
            #        $go_definition = $go_term->definition();
            #    }
            #    if (defined($go_name)) {
            #        push(@go_names, "$go_name [$go_id]");
            #    }
            #    else {
            #        push(@go_names, "$go_id");
            #    }
            #}
        }
        @uniprotkb_ac = @{get_unique(\@uniprotkb_ac)};
        $line_hash->{UniProtKB_AC} = $uniprotkb_ac[0];
        @uniprot_id = @{get_unique(\@uniprot_id)};
        $line_hash->{UniProt_ID} = join(',',@uniprot_id);
        #@omim_id = @{get_unique(\@omim_id)};
        #$line_hash->{OMIM_ID} = join(',',@omim_id);
        #@hgmd_id = @{get_unique(\@hgmd_id)};
        #$line_hash->{HGMD_ID} = join(',',@hgmd_id);
        #my $gene_ontology = join('|', @{get_unique(\@go_names)});
        #$line_hash->{Gene_Ontology} = $gene_ontology;

        #use hgnc (ensg preferred when a table is created that has it)
        if (defined $bvfoa_info{hgnc}){
            my @data;
            my $query = "CALL $vw::vw_database.hgnc2refseq_summary('$bvfoa_info{hgnc}')";
            my $qh = $vw::vw_conn->prepare($query);
            $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
            while (my @row = $qh->fetchrow_array()){
                push @data, $row[0];
            }
            $line_hash->{RefSeq_summary} = @data ? join('|', @data) : '';
        }
    }
    if (defined $bvfoa_info{hgnc}){
        my $query = "CALL $vw::vw_database.get_mitocarta_gene('$bvfoa_info{hgnc}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my @row = $qh->fetchrow_array();
        if(defined($row[0]) && $row[0] ne ''){
            $line_hash->{MT} = '1';
        }
        else{
            $line_hash->{MT} = '';
        }
    }
    return {};
}


1;

