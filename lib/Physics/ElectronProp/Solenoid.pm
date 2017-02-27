package Physics::ElectronProp::Solenoid;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use PDL;
use PDL::GSLSF::ELLINT;

sub _init {
  my $self = shift;
  $self->{num_loops} -= 1;
}

## Solenoid Lens Extra Parameters ##

sub lens_type	   { 'solenoid' 	}
sub num_loops      { $_[0]->{num_loops     } }
sub current	   { $_[0]->{current       }=$_[1] if defined $_[1]; $_[0]->{current       } }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  my $Btot;
  for my $n (0..$self->num_loops) {
    $Btot += $self->Bloop($r,$z,$n/$self->num_loops);
  }
  return $Btot;
}

sub E_Field {
  return zeros(3)
}

## Setup of solenoid ##

sub Bloop {
  my $self = shift;
  my ($r,$z,$n) = @_;
  my $Br = B_r($r, $z - $self->loop_step($n) - $self->front_pos, $self->lens_shape->($n));
  my $Bz = B_z($r, $z - $self->loop_step($n) - $self->front_pos, $self->lens_shape->($n));
  return pdl ($Br,0,$Bz)*$self->current;
}

sub loop_step {
  my $self = shift;
  my $n = shift;
  return 0 if $self->num_loops == 0;
  return $self->lens_length * $n ;
}

## Magnetic Field Functions ##

sub B_r {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  my $Br = mu_0 * (-$K +($a**2 + $r**2 + $z**2)/((1 - k($r,$z,$a)**2)*Q($r,$z,$a)) * $E) / (2 * pi * sqrt(Q($r,$z,$a)))*$z;
  unless (dims($Br)) {
    return ($Br==0) ? $Br : $Br/$r;  # I don't remember why this is necessary, but I am too afraid to change it...
  }
  $Br/=$r;
  $Br->inplace->setnantobad;
  $Br->inplace->setbadtoval(0);
  return $Br;
}

sub B_z {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a)); 
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  return mu_0 * ($K + ($a**2 - $r**2 - $z**2) * $E /((1 - k($r,$z,$a)**2)*Q($r,$z,$a))) / (2 * pi * sqrt(Q($r,$z,$a)));
}

sub Q {
  my ($r,$z,$a) = @_;
  return ($a + $r)**2 + $z**2;
}

sub k {
  my ($r,$z,$a) = @_;
  return sqrt(4*$a*$r/Q($r,$z,$a));
}

1;
