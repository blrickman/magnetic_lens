package Physics::ElectronProp::Simulation_Cart;

use warnings;
use strict;
use v5.10;
use PDL;
use Term::ProgressBar;
use File::chdir;

use Physics::ElectronProp::Solenoid;
use Physics::ElectronProp::RF_Cavity;
use Physics::ElectronProp::Generic_Lens;
use Physics::ElectronProp::Mesh_Lens;
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
  $self->sim_time(-$self->time_step) unless defined $self->sim_time;
  $self->{start_time} = ($self->sim_time);
  $_->history('time',-$self->time_step) for @{ $self->electrons };
}

## Sim Components ##

sub electrons 	{$_[0]->{electrons }}
sub lens 	{$_[0]->{lens }}

## Sim Parameters ##

sub z_start 	{ $_[0]->{z_start  }}
sub z_end   	{ $_[0]->{z_end    }}
sub r_end   	{ $_[0]->{r_end    }}
sub steps    	{ $_[0]->{steps	   }}
sub time_step	{ $_[0]->{time_step}}
sub sim_time	{ $_[0]->{sim_time} = $_[1] if defined $_[1]; $_[0]->{sim_time  }}
sub start_time	{ $_[0]->{start_time}}
sub prog_silent { $_[0]->{prog_silent}}
sub dir		{ $_[0]->{dir}	= $_[1] if defined $_[1]; $_[0]->{dir	}}
sub fields	{ $_[0]->{fields} = $_[1] if defined $_[1]; $_[0]->{fields  }}

sub evolve {
  my $self = shift;
  my $electron = shift;

  my $dt = $self->time_step;
  my $pos= $electron->position;
  my $vel= $electron->velocity;
  my $qe = $electron->charge;

#  $electron->history('pos_hist',$pos);
#  $electron->history('vel_hist',$vel);
#  $electron->history('time',$dt + $electron->history('time')->[-1]);
  my $t = $self->sim_time($dt+$self->sim_time);

  my $r = sqrt($pos->index(0)**2 + $pos->index(1)**2);
  my $theta = atan2($pos->index(1),$pos->index(0));
  my $z = $pos->index(2); 
  #$electron->min_rad([$r,$z]) if ($electron->min_rad->[0] > $r);
  
  ## B-Field ##

  my @Bfield = map { $_->B_Field(pdl ($r,$theta,$z,$t)) } @{ $self->lens };
  $Bfield[0] += pop @Bfield while @Bfield > 1;
  my ($Br,$Bt,$Bz) = dog($Bfield[0]);
  my ($Bx, $By) = ($Br * cos($theta) - $Bt * sin($theta), $Br * sin($theta) + $Bt * cos($theta));

  ## E-Field ##

  my @Efield = map { $_->E_Field(pdl ($r,$theta,$z,$t)) } @{ $self->lens };
  $Efield[0] += pop @Efield while @Efield > 1;
  my ($Er,$Et,$Ez) = dog($Efield[0]);
  my ($Ex, $Ey) = ($Er * cos($theta) - $Et * sin($theta), $Er * sin($theta) + $Et * cos($theta));

  $self->fields([$Ex,$Ey,$Ez,$Bx,$By,$Bz,$t]);

  ## Lorentz Force Calculation ##

  my $acc = (crossp($vel, pdl($Bx, $By, $Bz)) + pdl($Ex, $Ey, $Ez)) * $qe / ($electron->mass);
  $electron->accel( $acc );
  $electron->velocity( $vel + $acc * $dt );
  $electron->position( $pos + $vel * $dt + $acc * $dt**2 / 2 );

  #$self->near_lens($electron);
}

sub run {
  my $self = shift;
  my $CWD = $self->{dir};
  my $progress = Term::ProgressBar->new({
    count => @{ $self->electrons } * $self->steps, 
    ETA => 'linear',
    name => 'Progress',
    silent => $self->prog_silent,
  });
  my $l = 0;
  my $first = 0;
  for my $electron (@{ $self->electrons }) { 
    open my $DATA_OUT, "> $CWD/e" . $electron->{id} . ".dat";
    open my $FIELD_OUT, "> $CWD/e" .  $electron->{id} . "_fields.dat" unless $first;
    while ($electron->{position}->index(2) < $self->z_end) {
      $self->evolve( $electron );
      $self->export($electron,$DATA_OUT);
      $self->export_field($FIELD_OUT) unless $first;
      $progress->update($l) unless $l++ % 100; 
    }
    close $FIELD_OUT unless $first++;
    close $DATA_OUT;
    $self->sim_time($self->start_time);
  }
  $progress->update(@{ $self->electrons } * $self->steps);
  print "\n";
#  shift @{$_->history('time')} for @{ $self->electrons };
}

sub near_lens {
  my $self = shift;
  my $electron = shift;
  for (@{ $self->solenoids }) {
    if ($electron->position->slice(2) > $_->front_pos - $_->sol_length * 1.5 && ! $electron->near_lens($_) ) {
      print "Electron is approaching a lens, decreasing time step \n";
      $electron->near_lens($_,1);
      $self->{time_step} /= 10 ;
    } elsif ($electron->position->slice(2) > $_->front_pos + $_->sol_length * 1.1 && $electron->near_lens($_) == 1 ) {
      print "Electron is leaving a lens, increasing time step \n";
      $electron->near_lens($_,2);
      $self->{time_step} *= 10 ;
    }
  }
}

sub export {
  my $self = shift;
  my ($electron,$FH) = @_;
  my @export = (list($electron->{position}), list($electron->{velocity}), list($electron->{accel}), $self->{sim_time} );
  print $FH join(' ', @export) . "\n";
}

sub export_field {
  my $self = shift;
  my $FH = shift;
  print $FH join(' ',@{$self->fields}) . "\n";
}

1;

__END__
