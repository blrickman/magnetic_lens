#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Simulation_Cart;
use Physics::ElectronProp::Auxiliary ':constants';
use PDL::Util 'export2d';
use PDL;
use PDL::Graphics::Prima::Simple [1000,500];
use File::chdir;
use Getopt::Long;

my $export = 1;
my $message = '';

GetOptions(
  'export!'    => \$export,
  'message=s'  => \$message,
);

unless ($message) {
  print "Please enter message: ";
  chomp($message = <>);
}

my $front_rad	= 0.25*10**-2;
my $back_rad 	= 1.00*10**-2;
my $current	= 17.238;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * ($n);};
my $shape2 	= sub {my $n = shift; $back_rad + ($front_rad - $back_rad) * ($n);};

my $solenoid1 = Physics::ElectronProp::Solenoid->new(
  sol_name	=> 'con-f',
  front_radius 	=> $front_rad,
  back_radius	=> $back_rad,
  num_loops	=> 100,
  sol_length	=> 10**-2,
  current	=> $current,
  front_pos	=> 0,
  sol_shape	=> $shape,
  test 		=> 0,
);

my $solenoid2 = Physics::ElectronProp::Solenoid->new(
  sol_name	=> 'con-2b',
  front_radius 	=> $back_rad,
  back_radius	=> $front_rad,
  num_loops	=> 100,
  sol_length	=> 10**-2,
  current	=> $current,
  front_pos	=> (1 + 0.4524242490*2)*10**-2,
  sol_shape	=> $shape2,
  test 		=> 0,
);

my $electron1 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.1*10**-3,0,-5*10**-2],
);

my $electron2 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.25*10**-3,0,-5*10**-2],
);

my $electron3 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.5*10**-3,0,-5*10**-2],
);

my $sim = Physics::ElectronProp::Simulation_Cart->new(
  electrons	=> [$electron2,$electron3],#$electron1],#,
  solenoids	=> [$solenoid1,$solenoid2],
  z_start	=>  -5*10**-2,
  z_end		=>  10*10**-2,
  steps		=>  100000,
);

for my $solenoid (@{ $sim->solenoids }) {
  my $fn_sol = join('_',
  $solenoid->sol_name,
  $solenoid->front_radius . 'm',
  $solenoid->back_radius . 'm',
  $solenoid->sol_length . 'm',
  $solenoid->front_pos . 'm',
  );
  export(transpose(cat($solenoid->plot_lens)),$fn_sol . '.lns', 'data/lens_shape',0);
}

my $dir = join('_',
  $solenoid1->sol_name,
  $solenoid2->sol_name,
  $solenoid1->front_radius,
  $solenoid1->back_radius,
  $solenoid1->sol_length,
  $solenoid1->num_loops + 1 . 'lps',
  $solenoid1->current . 'A',
  's2-' .   $solenoid2->front_pos,
  $electron1->energy . 'eV',
  $sim->steps . 'it',
  'Cart',  
);
$dir = dir_exists($dir, 'data');
print "Writing data files to:\n$dir\n" if $export;


$sim->run();
my $i = 2;
for my $electron (@{ $sim->electrons }) {
  my $position = pdl $electron->history('pos_hist');
  my $velocity = pdl $electron->history('vel_hist');
  my $time     = transpose(pdl $electron->history('time'));
  export(append(append($position,$velocity),$time),"e" . $i++ . '.dat', "data/$dir",1) if $export;

  my @final_v = list(transpose($velocity)->slice(-1));
  my $vf = $final_v[0]**2 + $final_v[1]**2  + $final_v[2]**2;
  print "End Position is: {" . join(', ',list (1000*${$electron->history('pos_hist')}[-1])) ."} mm\n";
  print "End speed is " . sprintf ("%.4f",$vf**.5 / vc) . "c\n";
  print "The focus is " . sprintf ("%.4f", ($electron->min_rad->[1] - $solenoid1->sol_length) * 100) . "cm behind the lens\n"; 
  print "and gets focused to " . sprintf ("%.4f", $electron->min_rad->[0] * 10**6) . "microm\n\n"; 
}

sub dir_exists {
  my $dir = shift;
  local $CWD = shift;
  if (-d $dir) {
    warn "$dir already exists. Overwrite? (y/n)\n";
    while (<>) {
      return $dir if /y|Y/;
      print "Define a new dir\n";
      chomp($dir = <>);
      mkdir $dir;
      return $dir;
    }
  }
  mkdir $dir;
  return $dir;
}

sub export {
  my ($data,$fn) = (shift,shift);
  local $CWD = shift; 
  my $override = shift;
  if (!(-e $fn) || $override) {
    open my $fh, '>', $fn;
    print $fh "#" . $message . "\n";
    export2d($data, $fh);
  }
}

print "Runtime: " . (time - $^T) . "s\n";
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
