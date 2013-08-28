package Physics::ElectronProp::Simulation_Cart;

use warnings;
use strict;
use v5.10;
use PDL;

use Physics::ElectronProp::Solenoid;
use Physics::ElectronProp::Electron;
use Physics::ElectronProp::Auxiliary ':constants';

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
  $_->history('time',-$self->time_step) for @{ $self->electrons };
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
  my $pos= $electron->position;
  my $vel= $electron->velocity;
  my $qe = $electron->charge;

  $electron->history('pos_hist',$pos);
  $electron->history('vel_hist',$vel);
  $electron->history('time',$dt + $electron->history('time')->[-1]);

  my $r = sqrt($pos->index(0)**2 + $pos->index(1)**2);
  my $theta = atan2($pos->index(1),$pos->index(0));
  my $z = $pos->index(2); 
  $electron->min_rad([$r,$z]) if ($electron->min_rad->[0] > $r);
  
  my @field = map { $_->mag_tot(pdl ($r,0,$z)) } @{ $self->solenoids };
  $field[0] += pop @field while @field > 1;
  my ($Br,undef,$Bz) = dog($field[0]);
  my ($Bx, $By) = ($Br * cos($theta), $Br * sin($theta));

  my $acc = $self->mag_force($vel, pdl ($Bx, $By, $Bz) ) * $qe / ($electron->mass);
  $electron->velocity( $vel + $acc * $dt );
  $electron->position( $pos + $vel * $dt + $acc * $dt**2 / 2 );
}

sub run {
  my $self = shift;
  while (@{$self->electrons}[0]->{position}->index(2) < $self->z_end) {
    $self->evolve( $_ ) for @{ $self->electrons };
  }
  shift @{$_->history('time')} for @{ $self->electrons };
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
