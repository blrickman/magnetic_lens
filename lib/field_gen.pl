#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Simulation_Cart;
use Physics::ElectronProp::Auxiliary ':constants';
use PDL;
use Term::ProgressBar;
use Audio::Beep;
use File::chdir;
use Getopt::Long;
use File::Path 'make_path';

GetOptions(
  'directory=s'  	=> \my $basedir,
);

die "Enter a directoy!" unless defined $basedir;


my ($R,$L,$N) = (10,10,20);
my $dir = $basedir . "R" . sprintf("%02d", $R) . "mm_L" . sprintf("%02d", $L) . "mm_N$N";
make_path $dir;

my $zf  = 10.0*10**-2;
my $rad = 1.00*10**-3;
my $zs  = $zf  / 100;
my $rs  = $rad / 100;

my $front_rad	= $R*10**-3;
my $back_rad 	= $R*10**-3;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * ($n);};

my $lens = Physics::ElectronProp::Solenoid->new(
  lens_name	=> 'con-f',
  front_radius 	=> $front_rad,
  back_radius	=> $back_rad,
  num_loops	=> $N,
  lens_length	=> ($zf-$L*10**-3)/2,
  current	=> 1,
  front_pos	=> 0.0475,
  lens_shape	=> $shape,
  test 		=> 0,
);

my $progress = Term::ProgressBar->new({
  count => $zf/$zs, 
  ETA => 'linear',
  name => 'Progress',
});

my $fn = join "_", (0,$zf,$zs,0,$rad,$rs,$lens->front_pos);
my ($BZ, $BR);
{
  local $CWD = $dir;
  open $BZ, "> bz_$fn.dat" or die $!;
  open $BR, "> br_$fn.dat" or die $!;
} 

for my $z (0..$zf/$zs) {
  $z *= $zs;
  for my $r (0..$rad/$rs) {
    $r *= $rs;
    my $field = $lens->B_Field(pdl ($r,0,$z,0));
    print $BR sprintf("%.16e", $field->index(0)) . ", ";
    print $BZ sprintf("%.16e", $field->index(2)) . ", ";
  }
  print $_ "\n" for ($BZ,$BR);
  $progress->update($z/$zs);
}

my $beep = Audio::Beep->new();
my $music = "\\bpm400 g'' fes des a ges c ges' c.. ";
$beep->play($music);
