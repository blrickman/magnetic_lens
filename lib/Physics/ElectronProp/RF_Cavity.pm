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
sub epsilon	{ $_[0]->{mu		} }
sub mu		{ $_[0]->{epsilon	} }
sub radius	{ $_[0]->{radius	} }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z) = list $pos;
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z) = list $pos;
}

1;
