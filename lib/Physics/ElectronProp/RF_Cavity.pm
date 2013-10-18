package Physics::ElectronProp::RF_Cavity;

use parent EM_lens;
use PDL::GSLSF::BESSEL;

sub _init {
  my $self = shift;
  $self->front_radius{$self->radius};
  $self->back_radius{$self->radius};
  $self->B_Field(pdl [0,0,0]);
  $self->E_Field(pdl [0,0,0]);
}

## RF Cavity Lens Extra Parameters ##

sub E_0		{ $_[0]->{E_0		} }
sub mode	{ $_[0]->{mode		} }
sub epsilon	{ $_[0]->{epsilon	} }
sub mu		{ $_[0]->{mu		} }
sub radius	{ $_[0]->{radius	} }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my $mode = $self->mode
  my ($r,$theta,$z) = list $pos;
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z) = list $pos;
}

## Field Functions ##

sub TM_Fields {
  my $self = shift;
  my $pos  = shift;
  my $mode_nums = shift; 
  my ($R, $d) = map { $self->$_ } qw/radius lens_length/
  my $gamma = 
  my ($r,$theta,$z) = list $pos;
  my ($m, $n, $p) = list $mode_nums;
}

my %fields = {
  TM => \&TM_Fields,
  TE => \&TE_Fields,
}

sub gamma_2 {
  my $self = shift;
  return 
}

sub omega {
  my $self = shift;
  my ($m, $n, $p) = $self->mode;
  return ( sqrt( ( ()**2 + ()**2 ) / () ) ) 
}
1;
