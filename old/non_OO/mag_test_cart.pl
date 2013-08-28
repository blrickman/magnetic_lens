#!/usr/bin/perl

use warnings;
use strict;
use v5.10;
use PDL;
use PDL::GSLSF::ELLINT;
use PDL::Graphics::Prima::Simple [1000,500];

## Constants ##

my $pi 		= 4 * atan2(1,1);
my $charge 	= -1.6*10**-19	; # C
my $mass 	= 9.11*10**-31	; # kg
my $mu	   	= 4*$pi*10**-7	; # 
my $c		= 3*10**8	; # m/s

## Magnetic Lens Parameters ##

my $afront 	= 10**-2	; # m
my $aback  	= $afront	; # m
my $nloop  	= 0		; 
my $length 	= 10**-2	; # m 
my $Icur   	= 10		; # A

## Initial Conditions of Electron ##

my $KE		= 100*10**3	; # eV
my $vi   	= $c*(1-($KE/($mass*$c**2)
		  +1)**-2)**.5  ; # m/s
$vi		= 0.55*$c	; # m/s
my $xi		= 10**-2/40	; # m
my $yi		= 0		; # m
my $zi		= -5*10**-2	; # m
my $zf		=  10*10**-2	; # m

my $x_c		= pdl ($xi,$yi,$zi);
my $v_c		= pdl (0,0,$vi);

#my $B = mag_acc([0,0,0],[0,0,1]); print join(', ',@$B) . "\n";

## Simulation ##


say Btot(pdl [.00025,0,-0.05]);


############# TEST ################

sub test {
for my $r (0..100) {
  $r = ($r )/20000;
  for my $z (0..100) {
    $z = ($z - 50)/1000;
    say (join(', ', ($r,$z,list( Btot(pdl [$r,0,$z])) ) ));
  }
}



my $z = (sequence(101)-50)/1000;

my (@magr,@magz);
for (0..100) {
  $_ = ($_ -50)/1000;
  my ($mx,$my,$mz) = list Btot(pdl [5*10**-3,0,$_]);
  push @magr, $mx;
  push @magz, $mz;	
}

my $magr = pdl \@magr;
my $magz = pdl \@magz;
plot(
  -height1	=> ds::Pair($z, $magr,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -height2	=> ds::Pair($z, $magz,
    color 	=> cl::Blue,
    plotType 	=> ppair::Lines,
  ),
  x		=> {label => 'z (m)' },# , min => $zi, max => $zf},
  y		=> {label => 'B (T)' },# , min => -.00001, max => 1.1*$aback/10},
);
}


########################### FUNCTIONS ###########################


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
  return $afront if $nloop == 0;
  return $afront - ($afront - $aback) * $n / $nloop;
}

sub zstep {
  my $n = shift;
  return 0 if $nloop == 0;
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
