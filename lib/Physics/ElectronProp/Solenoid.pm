package Physics::ElectronProp::Solenoid;

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;
use Physics::ElectronProp::Auxiliary ':constants';

sub new {
  my $class = shift;
  my $self = {@_};
  bless($self, $class);
  $self->_init;
  return $self;
}

sub _init {
  my $self = shift;
}

## Magnetic Lens Parameters ##

sub front_diameter { $_[0]->{front_diameter}=$_[1] if defined $_[1]; $_[0]->{front_diameter} }
sub back_diameter  { $_[0]->{back_diameter }=$_[1] if defined $_[1]; $_[0]->{back_diameter } }
sub num_loops      { $_[0]->{num_loops     }=$_[1] if defined $_[1]; $_[0]->{num_loops     } }
sub sol_length     { $_[0]->{sol_length    }=$_[1] if defined $_[1]; $_[0]->{sol_length    } }
sub sol_name	   { $_[0]->{sol_name      }=$_[1] if defined $_[1]; $_[0]->{sol_name      } }
sub sol_shape	   { $_[0]->{sol_shape     }=$_[1] if defined $_[1]; $_[0]->{sol_shape     } }
sub current	   { $_[0]->{current       }=$_[1] if defined $_[1]; $_[0]->{current       } }
sub front_pos	   { $_[0]->{front_pos     }=$_[1] if defined $_[1]; $_[0]->{front_pos     } }

sub mag_field	   {$_[0]->{mag_field      }=$_[1] if defined $_[1]; $_[0]->{mag_field     } }

## Setup of solenoid ##

sub Bloop {
  my $self = shift;
  my ($r,$z,$n) = @_;
  my $Br = B_r($r, $z - $self->loop_step($n), $self->sol_shape($n));
  my $Bz = B_z($r, $z - $self->loop_step($n), $self->sol_shape($n));
  return pdl ($Br,$Bz)*$self->current;
}

sub loop_step {
  my $self = shift;
  my $n = shift;
  return $self->sol_length * $n / $self->num_loops;
}

sub mag_tot {
  my $self = shift;
  my ($r,$z) = @_;
  my $Btot;
  for my $n (0..$self->num_loops) {
    $Btot += $self->Bloop($r,$z,$n);
  }
  return [list $Btot];
}

## Magnetic Field Functions ##

sub B_r {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  my $Br = mu_0 * (-$K +($a**2 + $r**2 + $z**2)/((1 - k($r,$z,$a)**2)*Q($r,$z,$a)) * $E) / (2 * pi * sqrt(Q($r,$z,$a)))*$z;
  unless (dims($Br)) {
    return ($Br==0) ? $Br : $Br/$r;
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
