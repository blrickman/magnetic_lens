#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Electron;
use Physics::ElectronProp::Solenoid;
use PDL;
use PDL::Graphics::Prima::Simple;

use Data::Dumper;

my $electron = Physics::ElectronProp::Electron->new(
  velocity => [0,0,1],
  position => [0,0,0],
);

for (1..10) {
  my $x = pdl [0,0,$_];
  $electron->history('pos_hist',$x);
  $electron->position($x);
}

print Dumper($electron->history('pos_hist'));
