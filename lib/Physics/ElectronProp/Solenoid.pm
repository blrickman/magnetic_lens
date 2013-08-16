package Physics::ElectronProp::Solenoid;

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;
use Physics::ElectronProp::Auxiliary ':constants';

use Digest::SHA1 'sha1_hex';
use PDL::IO::Dumper;
use File::chdir;

my $cache_data = 1;

sub new {
  my $class = shift;
  my $self = {@_};
  bless($self, $class);
  $self->_init;
  return $self;
}

sub _init {
  my $self = shift;
  my $cache_filename = sha1_hex(join('-', 
    $self->front_diameter,
    $self->back_diameter,
    $self->num_loops,
    $self->sol_length,
    $self->current,
    $self->sol_name,
    $self->step_length, 
    $self->step_radial, 
    $self->mag_z_start, 
    $self->mag_z_end, 
    $self->mag_r_end,
  )) . ".dat";
  local $CWD = 'cache';
  if (-e $cache_filename & $cache_data) {
    print "Found cached magnetic field! Loading from $cache_filename\n";
    $self->mag_field(frestore($cache_filename)); 
  } else {
    $self->mag_field($self->generate_field());
    fdump($self->mag_field,$cache_filename);
    print "Magnetic field cached as $cache_filename\n";
  }
  $self->mag_cur_pos_adjust();
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

## Magnetic Field Mesh Creation ##

sub step_length    { $_[0]->{step_length}}
sub step_radial    { $_[0]->{step_radial}}
sub mag_z_start    { $_[0]->{mag_z_start}}
sub mag_z_end	   { $_[0]->{mag_z_end  }}
sub mag_r_end	   { $_[0]->{mag_r_end  }}

sub generate_field {
  my $self = shift;
  my $zstepsize = ($self->mag_z_end - $self->mag_z_start)/$self->step_length;
  my $rstepsize = ($self->mag_r_end)/$self->step_radial;
  my @mag_field;

  for my $z (0..$self->step_length) {
    my @z_col;
    $z = $z*$zstepsize + $self->mag_z_start;
    for my $r (0..$self->step_radial) {
      $r = $r*$rstepsize;
      push @z_col, [$r,$z,$self->mag_tot($r,$z)->[0],$self->mag_tot($r,$z)->[1]];
    }
    push @mag_field, \@z_col;
  }
  \@mag_field;
}

sub mag_cur_pos_adjust {
  my $self = shift;
  my $mag_field = $self->mag_field;

}

## Setup of solenoid ##

sub Bloop {
  my $self = shift;
  my ($r,$z,$n) = @_;
  my $Br = B_r($r, $z - $self->zstep($n), $self->sol_shape($n));
  my $Bz = B_z($r, $z - $self->zstep($n), $self->sol_shape($n));
  return pdl ($Br,$Bz)*$self->current;
}

sub zstep {
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
  return [$Btot->index(0),$Btot->index(1)];
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
