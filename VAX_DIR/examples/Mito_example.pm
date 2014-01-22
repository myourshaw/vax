package Mito;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);
use vw;

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub version { return '74'; }

sub feature_types { return ['Transcript']; }

sub get_header_info {
    my @new_output_cols = qw( MT );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);
    return { MT => "annotated as in mitochondrion by MitoCarta", };
}

sub run {
    my ($self, $tva, $line_hash) = @_;
    my $config = $self->{config};
	my $hgnc = $tva->transcript->{_gene_hgnc};	
    if (defined $hgnc){
        my $query = "CALL $vw::vw_database.get_mitocarta_gene('$hgnc')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my @row = $qh->fetchrow_array();
        if( defined($row[0]) && $row[0] ne '' ) {
            $line_hash->{MT} = '1';
        }
        else {
            $line_hash->{MT} = '';
        }
    }
    return {};
}

1;

# SQL stored procedure 
#     CREATE DEFINER=`sa`@`%` PROCEDURE `get_mitocarta_gene`(hgnc varchar(15))
#     BEGIN
#     SELECT `mitocarta_gene`.`mito_gene`
#     FROM `vw`.`mitocarta_gene`
#     WHERE `mitocarta_gene`.`mito_gene` = hgnc;
#     END
