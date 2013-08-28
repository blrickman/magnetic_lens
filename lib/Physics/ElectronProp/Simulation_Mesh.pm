package Physics::ElectronProp:Simulation;

use warnings;
use strict;
use v5.10;
use PDL;

use Physics::Electron::Solenoid;
use Physics::Electron::Electron;
use Physics::Electron::Auxiliary ':constants';

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

sub step_length    { $_[0]->{step_length}}
sub step_radial    { $_[0]->{step_radial}}
sub z_start	   { $_[0]->{mag_z_start}}
sub z_end	   { $_[0]->{mag_z_end  }}
sub r_end	   { $_[0]->{mag_r_end  }}
my $step 	= 500;
my $tstep 	= ($zf - $zi) / $vi /$step; # s

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

sub sim {

  my @x_store	= $x_c;
  my @v_store	= $v_c;


  for (0..$step) {
    my $ao = mag_acc($x_c,$v_c);
  
    #print "$x_c; $v_c; ".Btot($x_c)."\n";

    $x_c = $x_c + $v_c * $tstep + $ao * $tstep**2 / 2;
    $v_c = $ao * $tstep + $v_c;

    push @x_store, $x_c;
    push @v_store, $v_c;
  }

  my $x_data = pdl \@x_store;
  my $v_data = pdl \@v_store;

  #print $x_data;
  #print $v_data;

  my $x = $x_data->slice('(0),');
  my $y = $x_data->slice('(1),');
  my $z = $x_data->slice('(2),');
  return ($x,$y,$z);
}

sub mag_acc {
  my ($x,$v) = @_;
  my $B = Btot($x); 
  my $a = crossp($v, $B);
  $a *= $charge / $mass;
  return $a;
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

## Magnetic Field Mesh Creation ##








my $xi		= 10**-2/40	; # m
my $yi		= 0		; # m
my $zi		= -5*10**-2	; # m
my $zf		=  5*10**-2	; # m




