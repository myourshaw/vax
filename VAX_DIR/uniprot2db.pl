#!perl
use strict;
#note: no -w flag; suppress warnings because Swissknife_1.68/lib/SWISS/BaseClass.pm line 210 doesn't deal properly with undefs

#convert downloaded data into database friendly tab-delimited files

use Getopt::Long;
#
# SWISS::Entry is part of Swissknife
# Available from http://swissknife.sourceforge.net/
# See: http://swissknife.sourceforge.net/docs/
#
use SWISS::Entry;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#mysql import
#WORKDIR="/Volumes/scratch/uniprot";
#table=uniprot_human_protein;
#COLLIST=`head -1 $WORKDIR/$table.txt`
#mysql -u $USER -p$PASSWORD --execute "TRUNCATE TABLE vw.$table";
#mysql -u $USER -p$PASSWORD --execute "LOAD DATA LOCAL INFILE '$WORKDIR/$table.txt' INTO TABLE vw.$table IGNORE 1 LINES; SHOW WARNINGS;" > $table.output;
#table=uniprot_human_xref;
#COLLIST=`head -1 $WORKDIR/$table.txt`
#mysql -u $USER -p$PASSWORD --execute "TRUNCATE TABLE vw.$table";
#mysql -u $USER -p$PASSWORD --execute "LOAD DATA LOCAL INFILE '$WORKDIR/$table.txt' INTO TABLE vw.$table IGNORE 1 LINES; SHOW WARNINGS;" > $table.output;
#table=uniprot_human_feature;
#COLLIST=`head -1 $WORKDIR/$table.txt`
#mysql -u $USER -p$PASSWORD --execute "TRUNCATE TABLE vw.$table";
#mysql -u $USER -p$PASSWORD --execute "LOAD DATA LOCAL INFILE '$WORKDIR/$table.txt' INTO TABLE vw.$table IGNORE 1 LINES; SHOW WARNINGS;" > $table.output;

#-dir /scratch1/tmp/myourshaw/resources/uniprot/20131228


# configure from command line opts
my $config = &configure(scalar @ARGV);

&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;
    
    my $create_db_tables = 0;
    my $vwconn;
    # DBI CONFIG VARIABLES
    my $uniprot_entry_table = "uniprot_entry";
    my $uniprot_protein_table = "uniprot_protein";
    my $uniprot_feature_table = "uniprot_feature";
    my $uniprot_xref_table = "uniprot_xref";
    if($create_db_tables){
        use DBI;
        use DBD::mysql;
        my $platform = "mysql";
        my $host = "myourshaw-dev.genome.ucla.edu";
        my $port = "3306";
        my $database = "vw";
        my $user = "USER";
        my $pw = "PASSWORD";
        #print "user:";
        #my $user = <>;
        #print "password:";
        #my $pw = <>;
        my $dsn = "dbi:$platform:$database:$host:$port";
        $vwconn = DBI->connect($dsn, $user, $pw)
            or die "Unable to connect: $DBI::errstr\n";
    }
    
    # output columns
    #protein
    my @ENTRY_COLS = qw(
        UniProtKB_AC
        Entry
    );
    my @PROTEIN_COLS_COLS = qw(
        UniProtKB_AC
        Gene
        RecName
        ACs
        ID
        PE
        SQ
        DE
        GeneNames
        KW
        ALTERNATIVE_PRODUCTS
        ALLERGEN
        BIOPHYSICOCHEMICAL_PROPERTIES
        BIOTECHNOLOGY
        CATALYTIC_ACTIVITY
        CAUTION
        COFACTOR
        DEVELOPMENTAL_STAGE
        DISEASE
        DISRUPTION_PHENOTYPE
        DOMAIN
        ENZYME_REGULATION
        FUNCTION
        INDUCTION
        INTERACTION
        MASS_SPECTROMETRY
        MISCELLANEOUS
        PATHWAY
        PHARMACEUTICAL
        POLYMORPHISM
        PTM
        RNA_EDITING
        SEQUENCE_CAUTION
        SIMILARITY
        SUBCELLULAR_LOCATION
        SUBUNIT
        TISSUE_SPECIFICITY
        TOXIC_DOSE
        WEB_RESOURCE
        HGNC
        GO_term
        MIM_gene
        MIM_phenotype
        Reactome
        Pathway_Interaction
    );
     
    my @PROTEIN_KV_COLS = qw(
        UniProtKB_AC
        topic
        value
    );
    my @FEATURE_COLS = qw(
        UniProtKB_AC
        feature
        aaStart
        aaEnd
        description
    );
    my @XREF_COLS = qw(
        UniProtKB_AC
        RESOURCE_ABBREVIATION
        RESOURCE_IDENTIFIER
        OPTIONAL_INFORMATION_1
        OPTIONAL_INFORMATION_2
        OPTIONAL_INFORMATION_3
    );
    my $base = $config->{dir}."/";
    my $protein_cols_out = $base."uniprot_human_protein_cols.txt";
    my $protein_kv_out = $base."uniprot_human_protein.txt";
    my $feature_out = $base."uniprot_human_protein_feature.txt";
    my $xref_out = $base."uniprot_human_xref.txt";
    open PROTEIN_COLS,">",$protein_cols_out or die "can't open $protein_cols_out";
    print PROTEIN_COLS '#' . (join "\t", @PROTEIN_COLS_COLS) . "\n";
    open PROTEIN_KV,">",$protein_kv_out or die "can't open $protein_kv_out";
    print PROTEIN_KV '#' . (join "\t", @PROTEIN_KV_COLS) . "\n";
    open FEATURE,">",$feature_out or die "can't open $feature_out";
    print FEATURE '#' . (join "\t", @FEATURE_COLS) . "\n";
    open XREF,">",$xref_out or die "can't open $xref_out";
    print XREF '#' . (join "\t", @XREF_COLS) . "\n";
    
    foreach my $file($base.'uniprot_sprot_human.dat', $base.'uniprot_trembl_human.dat'){
        if(`gzip -t $file` == 0){
            open(DATA, "gunzip -c $file |") || die "can't open pipe to $file";
        }
        else {
            open(DATA, $file) || die "can't open $file";
        }
        #my $protein_out = $data.".protein.txt";
        #if ($create_db_tables){
        #    while(<DATA>){
        #        my $entry = SWISS::Entry->fromText($_);
        #        my $UniProtKB_AC = $entry->AC;
        #        my $Entry = $entry->toText();
        #        my $query = "INSERT INTO $uniprot_entry_table (UniProtKB_AC, Entry) VALUES ($UniProtKB_AC, $Entry)";
        #        my $query_handle = $vwconn->prepare($query);
        #        $query_handle->execute();
        #    }
        #    
        #    exit;
        #}
        #open PROTEIN,">",$protein_out or die "can't open $protein_out";
        #print PROTEIN '#' . (join "\t", @PROTEIN_COLS) . "\n";
        #my $feature_out = $data.".feature.txt";
        #open FEATURE,">",$feature_out or die "can't open $feature_out";
        #print FEATURE '#' . (join "\t", @FEATURE_COLS) . "\n";
        #my $xref_out = $data.".xref.txt";
        #open XREF,">",$xref_out or die "can't open $xref_out";
        #print XREF '#' . (join "\t", @XREF_COLS) . "\n";
        
        # Change the line termination string so we read an entire entry at a time
        local $/ = "\n//\n";
        my $fullParse=0; #lazy loading
        my $null = '';
        my $line_count = 0;
        
        while(<DATA>){
            $line_count++;
            my $entry = SWISS::Entry->fromText($_, $fullParse);
            my $primary_accession_number;
            my %line;
            my @tmp;
            
            #basic info
            $primary_accession_number = $entry->AC;
            $line{UniProtKB_AC} = $primary_accession_number;
            undef(@tmp);
            foreach my $ac ($entry->ACs->elements) {
                push @tmp, $ac;
            }
            $line{ACs} = join(";", @tmp);
            $line{ID} = $entry->ID;
            undef(@tmp);
            foreach my $id ($entry->IDs->elements) {
                push @tmp, $id;
            }
            #$line{IDs} = join(";", @tmp);
            $line{SQ} = $entry->SQ;
        
            #protein existence
            if (defined($entry->PE->text) && length(trim($entry->PE->text))>0){
                $line{PE} = trim($entry->PE->text);
                $line{PE} =~ s/: /:/;
                $line{PE} =~ s/;$//;
            }
            #my $name = $entry->DE;
            #description
            undef(@tmp);
            foreach my $de ($entry->DEs->elements) {
                if ($de->category eq "RecName" && $de->type eq "Full"){
                    $line{RecName} = $de->text;
                }
                push @tmp, $de->text;
            }
            $line{DE} = clean_join(\@tmp);
            
            #genes
            undef(@tmp);
            $line{Gene} = $entry->GNs->getFirst();
            foreach my $gns ($entry->GNs) {
                foreach my $geneGroups ($gns->list){
                    foreach my $geneGroup (@{$geneGroups}){
                        foreach my $names ($geneGroup->Names->list){
                            foreach my $name (@{$names}){
                                my $g = $name->text;
                                push @tmp, $g;
                            }
                       }
                    }
                }
            }
            $line{GeneNames} = clean_join(\@tmp);
            
            #xref
            my @ENST;
            my @ENSP;
            my @ENSG;
            my @KEGG;
            my @UCSC;
            my @RefSeq_NP;
            my @RefSeq_NM;
            my @HGNC;
            my @GO;
            my @GO_term;
            my @MIM_gene;
            my @MIM_phenotype;
            my @Reactome;
            my @Pathway_Interaction;
            undef(@tmp);
            foreach my $dr ($entry->DRs->elements){
                my ($Database_identifier, $primary_key, $secondary_key, $tertiary_key) = @{$dr}[0..3];
                #if (defined($tertiary_key)){
                #    print "$Database_identifier, $primary_key, $secondary_key, $tertiary_key\n";
                #}
                #else{
                #    print "$Database_identifier, $primary_key, $secondary_key\n";
                #}
                if ($Database_identifier eq "Ensembl"){
                    push @ENST, $primary_key;
                    push @ENSP, $secondary_key;
                    push @ENSG, $tertiary_key;
                }
                elsif ($Database_identifier eq "KEGG"){
                    push @KEGG, $primary_key;
                }
                elsif ($Database_identifier eq "UCSC"){
                    push @UCSC, $primary_key;
                }
                elsif ($Database_identifier eq "RefSeq"){
                    push @RefSeq_NP, $primary_key;
                    push @RefSeq_NM, $secondary_key;
                }
                elsif ($Database_identifier eq "HGNC"){
                    push @HGNC, $secondary_key; #official gene name
                }
                elsif ($Database_identifier eq "GO"){
                    push @GO, $primary_key;
                    push @GO_term, $secondary_key;
                }
                elsif ($Database_identifier eq "MIM"){
                    if ($secondary_key eq "gene"){
                        push @MIM_gene, $primary_key;
                    }
                    elsif ($secondary_key eq "phenotype"){
                        push @MIM_phenotype, $primary_key;
                    }
                    else{
                        push @MIM_gene, $primary_key;
                        push @MIM_phenotype, $primary_key;
                    }
                }
                elsif ($Database_identifier eq "Reactome"){
                    push @Reactome, $secondary_key;
                }
                elsif ($Database_identifier eq "Pathway_Interaction_DB"){
                    push @Pathway_Interaction, $secondary_key;
                }
                push @tmp, join(',', @{$dr});
            }
            $line{ENST} = clean_join(\@ENST);
            $line{ENSP} = clean_join(\@ENSP);
            $line{ENSG} = clean_join(\@ENSG);
            $line{KEGG} = clean_join(\@KEGG);
            $line{UCSC} = clean_join(\@UCSC);
            $line{RefSeq_NP} = clean_join(\@RefSeq_NP);
            $line{RefSeq_NM} = clean_join(\@RefSeq_NM);
            $line{HGNC} = clean_join(\@HGNC);
            $line{GO} = clean_join(\@GO);
            $line{GO_term} = clean_join(\@GO_term);
            $line{MIM_gene} = clean_join(\@MIM_gene);
            $line{MIM_phenotype} = clean_join(\@MIM_phenotype);
            $line{Reactome} = clean_join(\@Reactome);
            $line{Pathway_Interaction} = clean_join(\@Pathway_Interaction);
            $line{DRs} = clean_join(\@tmp);
        
            #keywords
            undef(@tmp);
            foreach my $kw ($entry->KWs->elements) {
                push @tmp, $kw->text;
            }
            $line{KW} = clean_join(\@tmp);
            
            #comments
            undef(@tmp);
            my %comments;
            for my $CC ($entry->CCs->elements()) {
                my $key = $CC->topic;
                if ($key eq 'Copyright'){
                    next;
                }
                my $comment = $CC->topic eq 'RNA EDITING' ?
                    $CC->note :
                    $CC->comment;        
                push(@{$comments{$key}}, $comment);
            }
            for my $key (keys %comments){
                $line{space_to_underscore($key)} = clean_join(\@{$comments{$key}})
            }
        
            my $output = join "\t", map { $line{$_} || $null } @PROTEIN_COLS_COLS;
            print PROTEIN_COLS "$output\n";
            
            for my $key(keys %line){
                if ($key eq 'UniProtKB_AC'){
                    next;
                }
                else{
                    print PROTEIN_KV "$primary_accession_number\t$key\t$line{$key}\n";
                }
            }
            
            #features
            my @FTs = $entry->FTs->elements();
            for my $ft(@FTs){
                my %line;
                $line{UniProtKB_AC} = $primary_accession_number;
                ($line{feature}, $line{aaStart}, $line{aaEnd}, $line{description}) = @{$ft}[0..3];
                my $output = join "\t", map { $line{$_} || $null } @FEATURE_COLS;
                print FEATURE "$output\n";
            }
            
            #xrefs
            foreach my $dr ($entry->DRs->elements){
                my %line;
                $line{UniProtKB_AC} = $primary_accession_number;
                ($line{RESOURCE_ABBREVIATION}, $line{RESOURCE_IDENTIFIER},
                    $line{OPTIONAL_INFORMATION_1}, $line{OPTIONAL_INFORMATION_2},
                    $line{OPTIONAL_INFORMATION_3}) = @{$dr};
                my $output = join "\t", map { $line{$_} || $null } @XREF_COLS;
                print XREF "$output\n";
            }
        } #while(<DATA>)
    print "$line_count records processed from $file\n";
    } #foreach my $file
    print "done\n";
} #main

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    GetOptions(
        $config,
        'host=s',
        'port=i',
        'user=s',
        'pass:s',
        'dir=s',
    );
    $config->{host} ||= "cortex.local";
    $config->{port} ||= 3306;
    $config->{user} ||= "vw";
    $config->{pass} ||= "vw";
    
    ## connect to databases
    #$config->{reg} = &connect_to_dbs($config);
    #
    #&get_adaptors($config);
   
    return $config;
}

#remove tabs newlines and pipes
sub clean($) {
    my $tmp = shift;
    $tmp =~ s/[\n\t\|]//g;
    return trim($tmp);
}
#clean up strings and join with pipe
sub clean_join(@) {
    my $tmp = shift;
    my @clean;
    for my $t (@{$tmp}){
        push @clean, clean($t)
    }
    return join("|", @clean);
}
# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($) {
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
sub space_to_underscore {
    my $text = shift;
    $text =~ s/\s{1,}/_/g;
    return $text;
}
sub get_unique {
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}
