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

## Simulation Ranges
my $zf  = 10.0*10**-2;
my $rad = 1.00*10**-3;
my $zs  = $zf  / 1000;
my $rs  = $rad / 1000;
my @fields = qw/ ez bt /;

## Lens Parameters
my ($R,$L,$N) = (10,10,20);
my $front_rad	= $R*10**-3;
my $back_rad 	= $R*10**-3;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * ($n);};
my @mode_value	= qw/ TM 0 1 0/;
my %mode;
@mode{ qw/ field m n p /} = @mode_value;

my $dir = $basedir . join('',@mode_value) . "_R" . sprintf("%02d", $R) . "mm_L" . sprintf("%02d", $L) . "mm";
make_path $dir;

my $lens2 = Physics::ElectronProp::Solenoid->new(
  front_radius 	=> $front_rad,
  back_radius	=> $back_rad,
  num_loops	=> $N,
  lens_length	=> $L*10**-3,
  current	=> 1,
  front_pos	=> ($zf-$L*10**-3)/2,
  lens_shape	=> $shape,
  test 		=> 0,
);

my $lens = Physics::ElectronProp::RF_Cavity->new(
  name		=> 'cav1',
  radius 	=> $front_rad,
  lens_length	=> $L*10**-3,
  E_0		=> 1,
  front_pos	=> ($zf-$L*10**-3)/2,
  mode		=> \%mode,
  epsilon	=> epsilon_0,
  mu		=> mu_0,
  phase		=> 0,
  test 		=> 1,
);

my $progress = Term::ProgressBar->new({
  count => $zf/$zs, 
  ETA => 'linear',
  name => 'Progress',
});

my $fn = join "_", (0,$zf,$zs,0,$rad,$rs,$lens->front_pos);
my %FILES;
{
  local $CWD = $dir;
  for ( @fields ) {
    open $FILES{$_}, "> ${_}_$fn.dat" or die $!;
  }
} 

for my $z (0..$zf/$zs) {
  $z *= $zs;
  for my $r (0..$rad/$rs) {
    $r *= $rs;
    my %field;
    @field{ qw/ er et ez br bt bz / } = (list( $lens->E_Field(pdl ($r,0,$z,0)) ),list($lens->B_Field( pdl ($r,0,$z,0)) ));
    for my $f (@fields) {
      my $FH = $FILES{$f}; 
      print $FH sprintf("%.16e", $field{$f}) . ", ";
    }
  }
  print $_ "\n" for map $FILES{$_}, @fields;
  $progress->update($z/$zs);
}

my $beep = Audio::Beep->new();
my $music = "\\bpm400 g'' fes des a ges c ges' c.. ";
$beep->play($music);
