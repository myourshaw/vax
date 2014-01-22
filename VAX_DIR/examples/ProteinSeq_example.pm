package ProteinSeq;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub version { return '74'; }

sub feature_types { return ['Transcript']; }

sub get_header_info {
    return { ProteinSeq => "amino acid sequence of transcript's translated protein", };
}

sub run {
    my ($self, $tva) = @_;
    if ( defined $tva->transcript->translation ) {
	    return { ProteinSeq => $tva->transcript->translation->seq() };
    }
    return {};
}

1;

