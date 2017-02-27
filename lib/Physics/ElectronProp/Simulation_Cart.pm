package Physics::ElectronProp::Simulation_Cart;

use warnings;
use strict;
use v5.10;
use PDL;
use Term::ProgressBar;
use File::chdir;

use Physics::ElectronProp::Solenoid;
use Physics::ElectronProp::RF_Cavity;
use Physics::ElectronProp::Dielectric_RF_Cavity;
use Physics::ElectronProp::Generic_Lens;
use Physics::ElectronProp::Mesh_Lens;
use Physics::ElectronProp::Aperture;
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
  if (defined $self->{time_step}) {
    print "The estimated number of steps per run is " . sprintf( "%.2e", ($self->{z_end} - $self->{z_start}) / $self->{time_step} / @{$self->electrons}[0]->{velocity}->index(2)) . "\n";
  } else {
    $self->{time_step} = ($self->{z_end} - $self->{z_start}) / $self->{steps} / @{$self->electrons}[0]->{velocity}->index(2)
  }
  $self->sim_time(-$self->time_step) unless defined $self->sim_time;
  $self->all_fields(0) unless defined $self->all_fields;
  $self->{start_time} = ($self->sim_time);
  $_->history('time',-$self->time_step) for @{ $self->electrons };
  $self->{debug} = 0 unless $self->{debug};
  if (defined $self->batch) {
    my ($batch,$num) = @{ $self->batch };
    my ($print,@keys) = @{ $self->bt_print };
    my ($subdir,@dkeys) = @{ $self->sub_dir };
    $self->{sub_dir} = sprintf($subdir,map { $batch->[$num]{$_} } @dkeys);
    $self->{bt_print} = sprintf($print,map { $batch->[$num]{$_} } @keys);
  }
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
sub batch	{ $_[0]->{batch    }}
sub bt_print	{ $_[0]->{bt_print }}
sub dir		{ $_[0]->{dir      }}
sub sub_dir	{ $_[0]->{sub_dir  }}
sub prog_silent { $_[0]->{prog_silent}}
sub fields	{ $_[0]->{fields} = $_[1] if defined $_[1]; $_[0]->{fields  }}
sub all_fields	{ $_[0]->{all_fields} = $_[1] if defined $_[1]; $_[0]->{all_fields  }}
sub debug	{ $_[0]->{debug    }}

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

  my $acc = (crossp($vel, pdl($Bx, $By, $Bz)) + pdl($Ex, $Ey, $Ez)) * $qe / ($electron->mass) / gamma($vel)**3;
  $electron->accel( $acc );
  $electron->velocity( $vel + $acc * $dt );
  $electron->position( $pos + $vel * $dt + $acc * $dt**2 / 2 );

  #$self->near_lens($electron);
}

sub run {
  my $self = shift;
  my $CWD = $self->{dir};
  my $progress = Term::ProgressBar->new({
    count => @{ $self->electrons } * 20, 
    ETA => 'linear',
    name => 'Progress',
    silent => $self->prog_silent,
  });
  my $l = 0;
  my $first = 0;
  for my $electron (@{ $self->electrons }) { 
    my $nd=1;
    open my $DATA_OUT, "> $CWD/e" . $electron->{id} . ".dat";
    open my $FIELD_OUT, "> $CWD/e" .  $electron->{id} . "_fields.dat" if !$first || $self->all_fields;
    while ($electron->{position}->index(2) < $self->z_end && $electron->{position}->index(2)  >=$self->z_start) {
      $self->evolve( $electron );
      $self->export($electron,$DATA_OUT);
      $self->export_field($electron,$FIELD_OUT) if !$first || $self->all_fields;
      my $dist = ($electron->{position}->index(2)-$self->z_start)/( $self->z_end-$self->z_start);
      $dist = $dist > 1 ? 1 : $dist;
      $progress->update($l + int (20*$dist)) if $progress->last_update < $l + int (20*$dist);
      if ($electron->{position}->index(2) > ($self->z_end - $self->z_start)*$nd/10 && $self->debug) {
        $progress->message( $electron->{position}->index(2) ); $nd++;
      }
    } 
    close $FIELD_OUT unless $first++ || $self->all_fields;
    close $DATA_OUT;
    $self->sim_time($self->start_time);
    $l += 20;
    $progress->message("Electron " . $l / 20 . " of " . @{ $self->electrons } . " is complete.");
  }
  $progress->update(@{ $self->electrons } * 20);
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
  my ($electron,$FH) = @_;
  print $FH join(' ',(list($electron->{position}),@{$self->fields})) . "\n";
}

sub gamma {
  my ($vel) = @_;
  return 1 / sqrt(1 - ($vel/vc)**2)
}

1;

__END__
