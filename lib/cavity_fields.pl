#!/usr/bin/perl

use warnings;
use strict;
use v5.10;
use PDL;
use File::chdir;
use File::Path 'make_path';

use Physics::ElectronProp::Auxiliary ':constants';
use Physics::ElectronProp::RF_Cavity;

## Set up Cavity and parameters ##

my $len = 1;
my $rad	= 1;
my $E0  = 1;
my $w	= 1*10**9;

for my $phase (0..3) {
$phase *= pi/2;

my $cavity = Physics::ElectronProp::RF_Cavity->new(
  sol_name	=> 'cav',
  radius 	=> $rad,
  lens_length	=> $len,
  front_pos	=> 0,
  mode		=> {field => 'TM', m => 0, n => 1, p => 0},
  epsilon	=> epsilon_0,
  mu		=> mu_0,
  omega		=> $w,
  E_0		=> $E0,
  phase		=> $phase,
);
$w = sprintf("%.0e", $w);

my $mode = $cavity->{mode}{field} . $cavity->{mode}{m} . $cavity->{mode}{n} . $cavity->{mode}{p};

## Set up simulation ##

my $t = 0;	#Time that removes temporal dependance of fields
my @r_heights = map {$_/10} (0..10);
my $z_step = 0.01;

my @fields = qw/ E_r E_t E_z B_r B_t B_z/;
my %FIELDFH;

my $time = $t == -99 ? "t-independent" : "phase-" . sprintf("%.1f", $cavity->phase / pi) . "pi";
my $dir = "fields/$mode/${w}w/$time";
make_path $dir;

for (@fields) {
  open $FIELDFH{$_}, "> $dir/$_.dat";
  print {$FIELDFH{$_}} <<END; 
# radius = $rad,
# length = $len,
# E-field = $E0,
# omega = $w,
# mode = $mode,
END
}
print {$FIELDFH{$_}} ("z, " . join(', ', @r_heights) . "\n") for @fields;

for my $z (0..$len/$z_step) {
  $z *= $z_step;
  print {$FIELDFH{$_}} ("$z, ") for @fields;
  for my $r (@r_heights) {
    my %field;
    my $pos = pdl [$r,0,$z,$t];
    my @field_val = list ( append($cavity->E_Field($pos),$cavity->B_Field($pos)));
    for (0..@fields-1) {
      $field{$fields[$_]} = $field_val[$_];
    }
    print {$FIELDFH{$_}} ($field{$_} . ", ") for @fields;
  }
  print {$FIELDFH{$_}} ("\n") for @fields;
}

## Set up Plotting ##

my $rs = @r_heights + 1;
my %y_range = (
  E_r => 1,
  E_t => 1,
  E_z => 1,
  B_r => 1,
  B_t => 1.3*10**-8,
  B_z => 1,
);

for my $field (@fields) {
  my $yrange = $y_range{$field};
  `gnuplot -e "dir = '$dir';
set output dir . '/../$field-$time.png';
datafile = dir . '/$field.dat';
set title '$field - $time';
height = $rs;
rangey = $yrange;" fields/fields.gp`;
}
}
