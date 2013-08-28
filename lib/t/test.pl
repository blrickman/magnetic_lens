#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Simulation;
use PDL::Util 'export2d';
use PDL;
use PDL::Graphics::Prima::Simple [1000,500];

use Data::Dumper;

my $solenoid = Physics::ElectronProp::Solenoid->new(
  front_radius 	=> 10**-2,
  back_radius	=> .25*10**-2,
  num_loops	=> 100,
  sol_length	=> 10**-2,
  current	=> 10,
  front_pos	=> 0,
  sol_shape	=> sub {my $n = shift; 10**-2 - (10**-2 - .25*10**-2) * $n;},
  test 		=> 0,
);

my ($lensx,$lensy) = $solenoid->plot_lens;

my $electron = Physics::ElectronProp::Electron->new(
  energy 	=> 100*10**3,
  position 	=> [.025*10**-2,0,-5*10**-2],
);

my $electron2 = Physics::ElectronProp::Electron->new(
  energy 	=> 100*10**3,
  position 	=> [.025*10**-2,0,-5*10**-2],
);

my $sim = Physics::ElectronProp::Simulation->new(
  electrons	=> [$electron],
  solenoids	=> [$solenoid],
  z_start	=>  -5*10**-2,
  z_end		=>  10*10**-2,
  steps		=>  100,
  func_time_step=>  [0,1,1,1],
);

my $sim2 = Physics::ElectronProp::Simulation->new(
  electrons	=> [$electron2],
  solenoids	=> [$solenoid],
  z_start	=>  -5*10**-2,
  z_end		=>  10*10**-2,
  steps		=>  100,
  func_time_step=>  [0,1,1,2],
);

$sim->run();
my $position = pdl $electron->history('pos_hist');
my ($r,undef, $z) = dog(transpose($position));

open my $fh, '>', 'concical_1_025_1_100_10A-100it.csv';
export2d($position, $fh);

$sim2->run();
my $position2 = pdl $electron2->history('pos_hist');
my ($r2,undef, $z2) = dog(transpose($position2));

open my $fh, '>', 'concical_1_025_1_100_10A-100it.csv';
export2d($position, $fh);

print "End Position is: " . ${$electron->history('pos_hist')}[-1] ."\n";


plot(
  -height1	=> ds::Pair($z*100,$r*1000,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -height2	=> ds::Pair($z2*100,$r2*1000,
    color 	=> cl::Blue,
    plotType 	=> ppair::Lines,
  ),
  -lens	=> ds::Pair($lensx*100,$lensy*1000,
    color 	=> cl::Black,
    plotType 	=> ppair::Lines,
  ),
  x		=> {label => 'z (cm)' },# , min => $zi, max => $zf},
  y		=> {label => 'r (mm)' , min => -.0001, max => .26},
);

