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
require configure;

my $export = 1;
my $dir = '';

GetOptions(
  'export!'    => \$export,
  'directory=s'  => \$dir,
);

unless ($dir) {
  print "Please enter a directory name: ";
  chomp($dir = <>);
}

$dir = dir_exists($dir, 'data');
print "Writing data files to:\n$dir\n" if $export;

my $sim = $configure::sim;
$sim->run();

## Export Data to a new directory
system("cp configure.pm data/$dir/configure.pm");
my $i = 1;
for my $electron (@{ $sim->electrons }) {
  my $position = pdl $electron->history('pos_hist');
  my $velocity = pdl $electron->history('vel_hist');
  my $time     = transpose(pdl $electron->history('time'));
  export(append(append($position,$velocity),$time),"e" . $i++ . '.dat', "data/$dir",1) if $export;

  #my @final_v = list(transpose($velocity)->slice(-1));
  #my $vf = $final_v[0]**2 + $final_v[1]**2  + $final_v[2]**2;
}

## Subroutines
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
    export2d($data, $fh);
  }
}
__END__
## Export Lens Shapes
for my $lens (@{ $sim->lens }) {
  my $fn_sol = join('_',
  $lens->{sol_name},
  $lens->{front_radius} . 'm',
  $lens->{back_radius} . 'm',
  $lens->{sol_length} . 'm',
  $lens->{front_pos} . 'm',
  );
  export(transpose(cat($lens->plot_lens)),$fn_sol . '.lns', 'data/lens_shape',0);
}

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
