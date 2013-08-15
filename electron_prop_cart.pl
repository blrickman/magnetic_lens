#!/usr/bin/perl

use warnings;
use strict;
use v5.10;
use PDL;
use PDL::GSLSF::ELLINT;
use PDL::Graphics::Prima::Simple [700,500];

## Constants ##

my $pi 		= 4 * atan2(1,1);
my $charge 	= -1; # -1.6*10**-19;
my $mass 	= 1;
my $mu	   	= 1; # 4*$pi*10**-7;

## Magnetic Lens Parameters ##

my $afront 	= 1;
my $aback  	= $afront;
my $nloop  	= 100;
my $length 	= 1;
my $Icur   	= 1;

## Initial Conditions of Electron ##

my $vi		= 50;
my $xi		= $aback/10;
my $yi		= 0;
my $zi		= -5 ;
my $zf		= 15 ;

my $x_c		= pdl ($xi,$yi,$zi);
my $v_c		= pdl (0,0,$vi);

#my $B = mag_acc([0,0,0],[0,0,1]); print join(', ',@$B) . "\n";

## Simulation ##

my $step 	= 500;
my $tstep 	= ($zf - $zi) / $vi /$step;

my ($x1,$y1,$z1) = sim();
my $r1 = sqrt($x1**2+$y1**2);

$xi		= $aback/20;
$x_c	  	 = pdl ($xi,$yi,$zi);
$v_c		 = pdl (0,0,$vi);
my ($x2,$y2,$z2) = sim();
my $r2 = sqrt($x2**2+$y2**2);

$x_c	  	 = pdl (0,0,$zi);
$v_c		 = pdl (0,0,$vi);
my ($x3,$y3,$z3) = sim();
my $r3 = sqrt($x3**2+$y3**2);

my $time = sequence($step+2)*$tstep;

plot(
  -height1	=> ds::Pair($z1,$r1,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -height2	=> ds::Pair($z2,$r2,
    color 	=> cl::Blue,
    plotType 	=> ppair::Lines,
  ),
  -height3	=> ds::Pair($z3,$r3,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
    lineStyle 	=> lp::Dash,
  ),
  x		=> {label => 'z', min => -5, max => 16},
  y		=> {label => 'r', min => -.001, max => .03},
);





########################### FUNCTIONS ###########################

sub sim {

  my @x_store	= $x_c;
  my @v_store	= $v_c;


  for (0..$step) {
    my $ao = mag_acc($x_c,$v_c);
  
    #print "$x_c; $v_c; ".Btot($x_c)."\n";

    $x_c = $x_c + $v_c * $tstep + $ao * $tstep**2 / 2;
    $v_c = $ao * $tstep + $v_c;

    push @x_store, $x_c;
    push @v_store, $v_c;
  }

  my $x_data = pdl \@x_store;
  my $v_data = pdl \@v_store;

  #print $x_data;
  #print $v_data;

  my $x = $x_data->slice('(0),');
  my $y = $x_data->slice('(1),');
  my $z = $x_data->slice('(2),');
  return ($x,$y,$z);
}

sub mag_acc {
  my ($x,$v) = @_;
  my $B = Btot($x); 
  my $a = crossp($v, $B);
  $a *= $charge / $mass;
  return $a;
}

## Setup of solenoid ##

sub Bloop {
  my ($pos,$n) = @_;
  my ($r,$x,$y,$z) = (sqrt(($pos->index(0))**2+($pos->index(1))**2),$pos->index(0),$pos->index(1),$pos->index(2));
  my $theta = atan2($y,$x);
  my $Br = Br($r, $z - zstep($n), a($n));
  my $Bz = Bz($r, $z - zstep($n), a($n));
  my $Bloop = pdl ($Br*cos($theta), $Br*sin($theta), $Bz);
  return $Bloop;
}

sub a {
  my $n = shift;
  return $afront - ($afront - $aback) * $n / $nloop;
}

sub zstep {
  my $n = shift;
  return $length * $n / $nloop;
}

sub Btot {
  my $pos = shift;
  my $Btot;
  for my $n (0..$nloop) {  ## CHANGE 1 to 0 FOR MULTIPLE LOOPS ##
    $Btot += Bloop($pos,$n);
  }
  return $Btot;
}

## Magnetic Field Functions ##

sub Br {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  my $Br = $mu * $Icur * (-$K +($a**2 + $r**2 + $z**2)/((1 - k($r,$z,$a)**2)*Q($r,$z,$a)) * $E) / (2 * $pi * sqrt(Q($r,$z,$a)))*$z;
  unless (dims($Br)) {
    return ($Br==0) ? $Br : $Br/$r;
  }
  $Br/=$r;
  $Br->inplace->setnantobad;
  $Br->inplace->setbadtoval(0);
  return $Br;
}

sub Bz {
  my ($r,$z,$a) = @_;
  #print k($r,$z,$a) .", ".Q($r,$z,$a)."\n";
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a)); 
  my ($E,) = gsl_sf_ellint_Ecomp(k($r,$z,$a));
  #print "$K, \t$E\n";
  return $mu * $Icur * ($K + ($a**2 - $r**2 - $z**2) * $E /((1 - k($r,$z,$a)**2)*Q($r,$z,$a))) / (2 * $pi * sqrt(Q($r,$z,$a)));
}

sub Q {
  my ($r,$z,$a) = @_;
  return ($a + $r)**2 + $z**2;
}

sub k {
  my ($r,$z,$a) = @_;
  return sqrt(4*$a*$r/Q($r,$z,$a));
}
