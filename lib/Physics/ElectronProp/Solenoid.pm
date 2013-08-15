package Physics::ElectronProp::Solenoid;

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;

sub new {
  my $class = shift;
  my $self = {@_};
  bless($self, $class);
  return $self;
}

## Magnetic Lens Parameters ##

sub front_diameter { $_[0]->{front_diameter}=$_[1] if defined $_[1]; $_[0]->{front_diameter} }
sub back_diameter  { $_[0]->{back_diameter }=$_[1] if defined $_[1]; $_[0]->{back_diameter } }
sub num_loops      { $_[0]->{num_loops     }=$_[1] if defined $_[1]; $_[0]->{num_loops     } }
sub sol_length     { $_[0]->{sol_length    }=$_[1] if defined $_[1]; $_[0]->{sol_length    } }
sub current	   { $_[0]->{current       }=$_[1] if defined $_[1]; $_[0]->{current       } }
sub front_pos	   { $_[0]->{front_pos     }=$_[1] if defined $_[1]; $_[0]->{front_pos     } }
sub shape	   { $_[0]->{shape         }=$_[1] if defined $_[1]; $_[0]->{shape         } }

## Setup of solenoid ##

sub Bloop {
  my ($pos,$n) = @_;
  my ($r,$x,$y,$z) = (sqrt(($pos->index(0))**2+($pos->index(1))**2),$pos->index(0),$pos->index(1),$pos->index(2));
  my $theta = atan2($y,$x);
  my $Br = Br($r, $z - zstep($n), a($n));
  my $Bz = Bz($r, $z - zstep($n), a($n));
  my $Bloop = pdl ($Br*cos($theta), $Br*sin($theta), $Bz);
  return $Bloop;
}

sub a {
  my $n = shift;
  return $afront - ($afront - $aback) * $n / $nloop;
}

sub zstep {
  my $n = shift;
  return $length * $n / $nloop;
}

sub Btot {
  my $pos = shift;
  my $Btot;
  for my $n (0..$nloop) {
    $Btot += Bloop($pos,$n);
  }
  return $Btot;
}

## Magnetic Field Functions ##

sub Br {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  my $Br = $mu * $Icur * (-$K +($a**2 + $r**2 + $z**2)/((1 - k($r,$z,$a)**2)*Q($r,$z,$a)) * $E) / (2 * $pi * sqrt(Q($r,$z,$a)))*$z;
  unless (dims($Br)) {
    return ($Br==0) ? $Br : $Br/$r;
  }
  $Br/=$r;
  $Br->inplace->setnantobad;
  $Br->inplace->setbadtoval(0);
  return $Br;
}

sub Bz {
  my ($r,$z,$a) = @_;
  #print k($r,$z,$a) .", ".Q($r,$z,$a)."\n";
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a)); 
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  #print "$K, \t$E\n";
  return $mu * $Icur * ($K + ($a**2 - $r**2 - $z**2) * $E /((1 - k($r,$z,$a)**2)*Q($r,$z,$a))) / (2 * $pi * sqrt(Q($r,$z,$a)));
}

sub Q {
  my ($r,$z,$a) = @_;
  return ($a + $r)**2 + $z**2;
}

sub k {
  my ($r,$z,$a) = @_;
  return sqrt(4*$a*$r/Q($r,$z,$a));
}

