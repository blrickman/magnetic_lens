package Physics::ElectronProp::Simulation_theta;

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
  for (@{ $self->electrons }) {
    my $electron = $_;
    $electron->previous_force(0);
    for (@{ $self->solenoids }) {
      $electron->near_lens($_,0);
    }
  }
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
sub time_step	{ $_[0]->{time_step} = $_[1] if defined $_[1]; $_[0]->{time_step} }
sub func_time_step  { $_[0]->{func_time_step}}


sub evolve {
  my $self = shift;
  my $electron = shift;

  my $dt    = $self->time_step;
  my $r     = $electron->position;
  my $v	    = $electron->velocity;

  $electron->history('pos_hist',$r);
  $electron->history('vel_hist',$v);
  
  my @force = map { crossp($v,$_) * $electron->charge } map { $_->mag_tot($r) } @{ $self->solenoids };
  $force[0] += pop @force while @force > 1; 
  my $acc = $force[0] / ($electron->mass);

  $acc->slice(0) += $v->slice(1)**2 / $r->slice(0);
  $acc->slice(1) -= 2 * $v->slice(0) * $v->slice(1);

  $self->check_negative_rad($electron) if $self->func_time_step->[0];
  $self->going_2_0($electron) if $self->func_time_step->[1];
  $self->near_lens($electron) if $self->func_time_step->[2];
  $self->check_force_magnitude($electron,$force[0]) if $self->func_time_step->[3];

  $electron->velocity( $v + $acc * $dt );
  $electron->position( $r + $v * $dt + $acc * $dt**2 / 2);
}

sub run {
  my $self = shift;
  while (@{$self->electrons}[0]->{position}->index(2) <= $self->z_end) {
    $self->evolve( $_ ) for @{ $self->electrons };
    #say ${$_->history('pos_hist')}[-10] for @{ $self->electrons };
  }
}

sub check_negative_rad {
  my $self = shift;
  my $electron = shift;
  if ($electron->position->slice(0) < 0) {
    my $backstep = 20;
    my $dec_time = 50;
    print "Running into negative r values, decreasing time step by $dec_time and taking $backstep steps back \n";
    $self->{time_step} /= $dec_time;
    $electron->velocity( ${$_->history('vel_hist')}[- $backstep]);
    $electron->position( ${$_->history('pos_hist')}[- $backstep]);
  }
}

sub going_2_0 {
  my $self = shift;
  my $electron = shift;
  if ($electron->position->slice(0) < $electron->min_rad) {
    $electron->{min_rad} /= 30;
    my $dec_time = 20;
    $self->{time_step} /= $dec_time;
    print "The electron is approaching zero, lowering time step \nElectron at z = ". $electron->position->slice(2)*100 . "cm\n";
  } elsif ($electron->position->slice(0) > $electron->min_rad * 50) {
    $electron->{min_rad} *= 50;
    my $inc_time = 20;
    $self->{time_step} *= $inc_time;
    print "The electron is moving away from zero, increasing time step \nElectron at z = ". $electron->position->slice(2)*100 . "cm\n";
  }
}

sub near_lens {
  my $self = shift;
  my $electron = shift;
  for (@{ $self->solenoids }) {
    if ($electron->position->slice(2) > $_->front_pos - $_->sol_length * .5 && ! $electron->near_lens($_) ) {
      print "Electron is approaching a lens, decreasing time step \n";
      $electron->near_lens($_,1);
      $self->{time_step} /= 10 ;
    } elsif ($electron->position->slice(2) > $_->front_pos + $_->sol_length * 1.5 && $electron->near_lens($_) == 1 ) {
      print "Electron is leaving a lens, increasing time step \n";
      $electron->near_lens($_,2);
      $self->{time_step} *= 10 ;
    }
  }
}

sub check_force_magnitude {
  my $self = shift;
  my $electron = shift;
  my @force = list shift; 
  @force = map {$_**2} @force;
  my $force;
  $force += pop @force while @force > 1;
  $force = $force**.5;
  if ($force > $electron->previous_force * 100) {
    print "Force on electron is growing stronger, decreasing time step \nElectron at z = ". $electron->position->slice(2)*100 . "cm\n";
    $electron->previous_force($force);
    $self->{time_step} /= 10 ;
  } elsif ($force < $electron->previous_force /100) {
    print "Force on electron is growing weaker, increasing time step\nElectron at z = ". $electron->position->slice(2)*100 . "cm\n";
    $electron->previous_force($force);
    $self->{time_step} *= 10 ;
  }
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
