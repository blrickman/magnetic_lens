#!/usr/bin/perl

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;
use PDL::Graphics::Prima::Simple;

## Constants ##

my $pi 		= 4 * atan2(1,1);
my $charge 	= 1;
my $mass 	= 1;
my $mu	   	= 1;

## Magnetic Lens Parameters ##

my $afront 	= 1;
my $aback  	= $afront;
my $nloop  	= 1;
my $length 	= 0;
my $Icur   	= 1;

## Initial Conditions of Electron ##

my $vi		= 10;
my $xi		= $aback/2;
my $yi		= 0;
my $zi		= -10 ;
my $zf		=  10 ;

my $x_c		= pdl ($xi,$yi,$zi);
my $v_c		= pdl (0,0,$vi);
my @x_store	= $x_c;
my @v_store	= $v_c;

#my $B = mag_acc([0,0,0],[0,0,1]); print join(', ',@$B) . "\n";

## Simulation ##

my $step 	= 300;
my $tstep 	= ($zf - $zi) / $vi /$step;

my $z = sequence(100)/10 -5;
my $x = ones(100)*0;
my $y = -ones(100)*0;
my $pos = transpose pdl ($x,$y,$z);

my $B = Btot($pos);
my ($Bx,$By,$Bz) = dog $B;

plot(
  -by		=> ds::Pair($z,$By,
    color 	=> cl::Red,
  ),
  -bx		=> ds::Pair($z,$Bx,
    color 	=> cl::Blue,
  ),
  -bz		=> ds::Pair($z,$Bz,
    color 	=> cl::Green,
  ),
  x		=> {label => 'z',},# min => -1, max => -.75},
  y		=> {label => 'B',},# min => 0.024, max => .026},
);








########################### FUNCTIONS ###########################

sub sim {
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
  for my $n (1..$nloop) {  ## CHANGE 1 to 0 FOR MULTIPLE LOOPS ##
    $Btot += Bloop($pos,$n);
  }
  return $Btot;
}

## Magnetic Field Functions ##

sub Br {
  my ($r,$z,$a) = @_;
  #$r = ($r==0) ? 10**(-17) : $r;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
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
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a)); 
  my ($E,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  return $mu * $Icur * ($K + ($a**2 - $r**2 - $z**2) * $E /((1 - k($r,$z,$a)**2)*Q($r,$z,$a))) / (2 * $pi * sqrt(Q($r,$z,$a)));
}

sub Q {
  my ($r,$z,$a) = @_;
  return ($a + $r)**2 + $z**2;
}

sub k {
  my ($r,$z,$a) = @_;
  return sqrt( 4 * $a * $r / Q($r,$z,$a) );
}
