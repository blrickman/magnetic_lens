#!/usr/bin/perl

use warnings;
use strict;
use v5.10;
use Physics::ElectronProp::Electron;
use Physics::ElectronProp::Solenoid;
use PDL;
use PDL::Graphics::Prima::Simple;

use Data::Dumper;


my $solenoid = Physics::ElectronProp::Solenoid->new(
  front_diameter=> 10**-2,
  back_diameter	=> 10**-2,
  num_loops	=> 0,
  sol_length	=> 10**-2,
  sol_shape	=> sub {my $n = shift; 10**-2;},
  current	=> 10,
  front_pos	=> 0,
);

for my $r (0..100) {
  $r = ($r)/20000;
  for my $z (0..100) {
    $z = ($z - 50)/1000;
    say (join(', ', ($r,$z,list( $solenoid->mag_tot(pdl [$r,0,$z])) ) ));
  }
}

my @mag;

for (0..99) {
  $_ -= 50;
  $_ /= 1000;
  my $pos = pdl [5*10**-3,0,$_];
  push @mag, $solenoid->mag_tot($pos);
}

my $z = (sequence(100)-50)/1000;

my ($magr,undef, $magz) = dog(transpose( pdl \@mag));	


plot(
  -height1	=> ds::Pair($z,$magr,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -height2	=> ds::Pair($z,$magz,
    color 	=> cl::Blue,
    plotType 	=> ppair::Lines,
  ),
  x		=> {label => 'z (m)' },# , min => $zi, max => $zf},
  y		=> {label => 'B (T)' },# , min => -.00001, max => 1.1*$aback/10},
);
