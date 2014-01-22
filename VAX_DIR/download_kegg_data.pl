#!/usr/bin/perl -w
use strict;

use File::Path qw(make_path);
use File::Spec qw(splitpath);
use LWP::Simple;
use DBI;
use DBD::mysql;
use Getopt::Long;
Getopt::Long::Configure('no_ignore_case');

#-o /share/apps/myourshaw/resources/kegg/20131228/kegg_gene_pathways.txt
#-o /share/apps/aliz/resources/kegg/20131118/kegg_gene_pathways.txt
 
sub unquote($){
    my $str = shift;
    $str =~ s/'/\\'/g;
    return $str;
}
 
my %options;
GetOptions(
    'o|output=s' => \$options{output},
    'host' => \$options{host},
    'port' => \$options{port},
    'u|user=s' => \$options{user},
    'p|password=s' => \$options{password},
    'platform=s' => \$options{platform},
    'database=s' => \$options{database},
    'table=s' => \$options{table},
);
$options{host} ||= "cortex.local";
$options{port} ||= "3306";
$options{platform} ||= "mysql";
$options{database} ||= "vw";
$options{table} ||= "kegg_gene_pathway";
 
my $create_tables = 0;
if (defined($options{host}) && defined($options{user}) && defined($options{password})){
    $create_tables = 1;
}
 
my $dsn = "dbi:$options{platform}:$options{database}:$options{host}:$options{port}";
my $insert = "INSERT INTO $options{table} (gene_id, path_id, pathway) VALUES ('%s', '%s', '%s')";
my $truncate = "TRUNCATE TABLE $options{table}";
my $conn;
if ($create_tables){
    $conn = DBI->connect($dsn, $options{user}, $options{password})
        or die "Unable to connect: $DBI::errstr\n";
    my $qh = $conn->prepare($truncate);
    $qh->execute() or die "Unable to execute $truncate: $DBI::errstr\n";
}
 
my $write_files = 0;
if (defined($options{output})){
    $write_files = 1;
}
if ($write_files){
    my ($volume,$directories,$file) = File::Spec->splitpath($options{output});
    make_path($directories);
    open OUT, ">", $options{output}
        or die "can't open $options{output}\n";
    print OUT "#gene_id\tpath_id\tpathway\n";
}
 
my $url = 'http://rest.kegg.jp';
my $pathslist = get("$url/list/pathway/hsa");
my @paths = split("\n", $pathslist);
foreach my $path (@paths) {
    my($path_id, $pathway) = split(/\t/, $path, 2);
    $path_id =~ s/path://;
    $pathway =~ s/ - Homo sapiens \(human\)$//;
    my $geneslist = get("$url/link/genes/$path_id");
    my @genes = split("\n", $geneslist);
    foreach my $gene (@genes){
	my($path_id, $gene_id) = split(/\t/, $gene, 2);
        if ($create_tables){
            my $query = sprintf($insert, ($gene_id, $path_id, unquote($pathway)));
            my $qh = $conn->prepare($query);
            $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        }
        if ($write_files){
            print OUT "$gene_id\t$path_id\t$pathway\n";
        }
    }
}
