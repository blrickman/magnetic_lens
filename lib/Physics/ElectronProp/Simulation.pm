package Physics::ElectronProp::Simulation;

use warnings;
use strict;
use v5.10;
use PDL;

use Physics::ElectronProp::Solenoid;
use Physics::ElectronProp::Electron;
use Physics::ElectronProp::Auxiliary ':constants';

use Digest::SHA1 'sha1_hex';
use PDL::IO::Dumper;
use File::chdir;

sub new {
  my $class = shift;
  my $self = {@_};
  bless($self, $class);
  $self->_init;
  return $self;
}

sub _init{
  my $self = shift;
  $self->{time_step} = ($self->{z_end} - $self->{z_start}) / $self->{steps} / @{$self->electrons}[0]->{velocity}->index(2) unless defined $self->{time_step}; 
}

## Sim Components ##

sub electrons 	{$_[0]->{electrons }}
sub solenoids 	{$_[0]->{solenoids }}

## Sim Parameters ##

sub z_start 	{ $_[0]->{z_start  }}
sub z_end   	{ $_[0]->{z_end    }}
sub r_end   	{ $_[0]->{r_end    }}
sub steps    	{ $_[0]->{steps	   }}
sub time_step	{ $_[0]->{time_step}}


sub evolve {
  my $self = shift;
  my $electron = shift;

  my $dt = $self->time_step;
  my $r  = $electron->position;
  my $v  = $electron->velocity;
  my $qe = $electron->charge;

  $electron->history('pos_hist',$r);
  $electron->history('vel_hist',$v);

  my @force = map { $self->mag_force($v,$_) * $qe } map { $_->mag_tot($r) } @{ $self->solenoids };
  $force[0] += pop @force while @force > 1; 
  my $acc = $force[0] / ($electron->mass);

  $electron->velocity( $v + $acc * $dt );
  $electron->position( $r + $v * $dt + $acc * $dt**2 / 2);
}

sub run {
  my $self = shift;
  while (@{$self->electrons}[0]->{position}->index(2) < $self->z_end) {
    $self->evolve( $_ ) for @{ $self->electrons };
  }
}

sub mag_force {
  my $self = shift;
  my $v = shift;
  my $B = shift;
  my $force = crossp($v, $B);
  return $force;
}

1;

__END__
## Mesh Stuff ##


sub step_length    { $_[0]->{step_length}}
sub step_radial    { $_[0]->{step_radial}}

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
  if (-e $cache_filename && $cache_data) {
    print "Found cached magnetic field! Loading from $cache_filename\n";
    $self->mag_field(frestore($cache_filename)); 
  } else {
    $self->mag_field($self->generate_field());
    fdump($self->mag_field,$cache_filename);
    print "Magnetic field cached as $cache_filename\n";
  }
  $self->mag_cur_pos_adjust();
}

sub generate_field {
  my $self = shift;
  my $zstepsize = ($self->z_end - $self->z_start)/$self->step_length;
  my $rstepsize = ($self->r_end)/$self->step_radial;
  my @mag_field;

  for my $z (0..$self->step_length) {
    my @z_col;
    $z = $z*$zstepsize + $self->z_start;
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
