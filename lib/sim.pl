#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Simulation_Cart;
use Physics::ElectronProp::Auxiliary ':constants';
use PDL::Util 'export2d';
use PDL;
use PDL::Graphics::Prima::Simple [1000,500];
use File::chdir;

my $no_export = shift;

my $front_rad	= .625*10**-2;
my $back_rad 	= .625*10**-2;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * $n;};

my $solenoid = Physics::ElectronProp::Solenoid->new(
  sol_name	=> 'cylindrical',
  front_radius 	=> $front_rad,
  back_radius	=> $back_rad,
  num_loops	=> 100,
  sol_length	=> 10**-2,
  current	=> 19.25,
  front_pos	=> 0,
  sol_shape	=> $shape,
  test 		=> 0,
);

my $electron1 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.025*10**-2,0,-5*10**-2],
);

my $electron2 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.025*10**-2*2/3,0,-5*10**-2],
);

my $electron3 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.025*10**-2/3,0,-5*10**-2],
);

my $sim = Physics::ElectronProp::Simulation_Cart->new(
  electrons	=> [$electron1],#,$electron2,$electron3],
  solenoids	=> [$solenoid],
  z_start	=>  -5*10**-2,
  z_end		=>  10*10**-2,
  steps		=>  100000,
);

my $fn_sol = join('_',
  $solenoid->sol_name,
  $solenoid->front_radius,
  $solenoid->back_radius,
  $solenoid->sol_length,
);
export(transpose(cat($solenoid->plot_lens)),$fn_sol . '.lns', 'data/lens_shape',0);

my $fn = join('_',
  $fn_sol,
  $solenoid->num_loops + 1 . 'lps',
  $solenoid->current . 'A',
  $electron1->energy . 'eV',
  $sim->steps . 'it',
  'Cart.dat',  
);
$fn = fn_exists($fn,"e1-", 'data');
print "Writing data files to $fn\n";


$sim->run();
my $i = 1;
for my $electron (@{ $sim->electrons }) {
  my $position = pdl $electron->history('pos_hist');
  my $velocity = pdl $electron->history('vel_hist');
  my $time     = transpose(pdl $electron->history('time'));
  export(append(append($position,$velocity),$time),"e" . $i++ . "-" . $fn, 'data',1) unless $no_export;

  my @final_v = list(transpose($velocity)->slice(-1));
  my $vf = $final_v[0]**2 + $final_v[1]**2  + $final_v[2]**2;
  print "End Position is: {" . join(', ',list (1000*${$electron->history('pos_hist')}[-1])) ."} mm\n";
  print "End speed is " . sprintf ("%.4f",$vf**.5 / vc) . "c\n";
  print "The focus is " . sprintf ("%.4f", ($electron->min_rad->[1] - $solenoid->sol_length) * 100) . "cm behind the lens\n"; 
  print "and gets focused to " . sprintf ("%.4f", $electron->min_rad->[0] * 10**6) . "microm\n\n"; 
}

sub fn_exists {
  my ($fn,$prefix) = (shift,shift);
  local $CWD = shift;
  if (-e $fn) {
    warn "$fn already exists. Overwrite? (y/n)\n";
    while (<>) {
      return $fn if /y|Y/;
      print "Define a new fn\n";
      $fn = <>;
      chomp $fn;
      return $fn;
    }
  }
  return $fn;
}

sub export {
  my ($data,$fn) = (shift,shift);
  local $CWD = shift; 
  my $override = shift;
  if (!(-e $fn) || $override) {
    open my $fh, '>', $fn;
    export2d($data, $fh);
  }
}
__END__
my ($r,undef, $z) = dog(transpose($position));
plot(
  -height1	=> ds::Pair($z*100,$r*1000,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -lens	=> ds::Pair($lensx*100,$lensy*1000,
    color 	=> cl::Black,
    plotType 	=> ppair::Lines,
  ),
  x		=> {label => 'z (cm)' },# , min => $zi, max => $zf},
  y		=> {label => 'r (mm)' , min => -.0001, max => .26},
);
