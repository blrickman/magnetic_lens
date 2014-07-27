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
my ($R,$L,$N,$hl) = (10,10,20,1);
my $zf  = $L*10**-3;
my $rad = 3*$hl*10**-3;
my $zs  = $zf  / 100;
my $rs  = $rad / 100;
my @fields = qw/ ez er /;

## Lens Parameters
my $front_rad	= $R*10**-3;
my $back_rad 	= $R*10**-3;
my $shape 	= sub {my $n = shift; $front_rad - ($front_rad - $back_rad) * ($n);};
my $E0		= .85*10**8;
my @mode_value	= qw/ TM 0 1 0/;
my %mode;
@mode{ qw/ field m n p /} = @mode_value;

my $subdir = "Hole" . $hl . "mm_R" . sprintf("%02d", $R) . "mm_L" . sprintf("%02d", $L) . "mm";
#join('',@mode_value) . "_R" . sprintf("%02d", $R) . "mm_L" . sprintf("%02d", $L) . "mm";

my $dir = $basedir . $subdir;
make_path $dir;

my $lens = Physics::ElectronProp::Aperture->new(
  name		=> 'ap',
  radius 	=> $R/1000,
  lens_length	=> $R/1000,
  front_pos	=> -$L*10**-3/2,
  ap_radius	=> $hl/1000,
  E_0		=> -$E0,
  E_1		=> 0,
  omega		=> 1,
  phase		=> 0,
  test		=> 1,
);

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

my $lens1 = Physics::ElectronProp::RF_Cavity->new(
  name		=> 'cav1',
  radius 	=> $front_rad,
  lens_length	=> $L*10**-3,
  E_0		=> $E0,
  front_pos	=> -$L*10**-3/2,
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

my $fn = join "_", (-$zf,$zf,$zs,0,$rad,$rs,$lens->front_pos);
my %FILES;
{
  local $CWD = $dir;
  for ( @fields ) {
    open $FILES{$_}, "> ${_}_$fn.dat" or die $!;
  }
} 

my @lenses = ($lens,$lens1);

for my $z (-$zf/$zs..$zf/$zs) {
  $z *= $zs;
  for my $r (0..$rad/$rs) {
    $r *= $rs;
    my ($E,$B);
    for (@lenses) {
      $E += $_->E_Field(pdl ($r,0,$z,0));
      $B += $_->B_Field(pdl ($r,0,$z,0));
    }
    my %field;
    @field{ qw/ er et ez br bt bz / } = (list( $E ),list( $B ));
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
