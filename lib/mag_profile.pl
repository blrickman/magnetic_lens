#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Solenoid;
use PDL;
use File::chdir;

my $front_rad	= 10**-2;
my $back_rad 	=.25*10**-2;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * $n;};

my $solenoid = Physics::ElectronProp::Solenoid->new(
  sol_name	=> 'conical',
  front_radius 	=> $front_rad,
  back_radius	=> $back_rad,
  num_loops	=> 100,
  sol_length	=> 10**-2,
  current	=> 10,
  front_pos	=> 0,
  sol_shape	=> $shape,
  test 		=> 0,
);

my $column_names = "r, z, Br, Btheta, Bz";

my $steps = 200;


my $filename = join('_',
  $solenoid->sol_name,
  $solenoid->front_radius,
  $solenoid->back_radius,
  $solenoid->num_loops,
  $solenoid->sol_length,
  $solenoid->current,
  $steps,
);

generate(); 

sub generate {
  local $CWD = 'solenoid-data'; 
  mkdir $filename;
  local $CWD = $filename;

  for my $r (0..9) {
    $r *= $solenoid->back_radius / 10;
    open my $fh, '>', "r-". sprintf ("%.5f",$r) .".dat";
#    print $fh "$column_names\n";
    for my $z (0..$steps) {
      $z = 5*10**-2 * ($z/$steps - 0.5);
      print $fh "$r, $z, " . join(", ",list($solenoid->mag_tot(pdl [$r,0, $z]))) . "\n";
    }
    print "r = $r m is finished\n";
  }
}
