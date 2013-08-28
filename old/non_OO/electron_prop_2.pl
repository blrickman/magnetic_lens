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
my $aback  	= $afront/4;
my $nloop  	= 1;
my $length 	= 0;
my $Icur   	= 1;

## Initial Conditions of Electron ##

my $vi		= 1;
my $xi		= $aback/10;
my $yi		= 0;
my $zi		= -1 ;
my $zf		=  2 ;

my $x_c		= pdl ($xi,$yi,$zi);
my $v_c		= pdl (0,0,$vi);
my @x_store	= $x_c;
my @v_store	= $v_c;

#my $B = mag_acc([0,0,0],[0,0,1]); print join(', ',@$B) . "\n";

## Simulation ##

my $step 	= 100;
my $tstep 	= ($zf - $zi) / $vi /$step;

for (0..$step) {
  my $ao = mag_acc($x_c,$v_c);
  
  #print "$x_c; $v_c\n";

  $x_c = $x_c + $v_c * $tstep + $ao * $tstep**2 / 2;
  $v_c = $ao * $tstep + $v_c;

  push @x_store, $x_c;
  push @v_store, $v_c;
}

my $x_data = pdl \@x_store;
my $v_data = pdl \@v_store;

print $x_data;

print $v_data;

my $x = $x_data->slice('(0),');
my $y = $x_data->slice('(1),');
my $z = $x_data->slice('(2),');

line_plot($x,$z);

sub mag_acc {
  my ($x,$v) = @_;
  my $B = Btot($x); 
  my $a = crossp $v, $B;
  $a *= $charge / $mass;
  return $a;
}

## Setup of solenoid ##

sub Bloop {
  my ($x,$n) = @_; 
  my $Br = Br($x->index(0), $x->index(2) - zstep($n), a($n));
  my $Bz = Bz($x->index(0), $x->index(2) - zstep($n), a($n));
  my $Bloop = pdl ($Br , 0, $Bz);
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
  my $x = shift;
  my $Btot;
  for my $n (1..$nloop) {  ## CHANGE FOR MULTIPLE LOOPS ##
    $Btot += Bloop($x,$n);
  }
  return $Btot;
}

## Magnetic Field Functions ##

sub Br {
  my ($r,$z,$a) = @_;
  my ($K,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  my ($E,) = gsl_sf_ellint_Kcomp(k($r,$z,$a));
  return $mu * $Icur * (-$K +($a**2 + $r**2 + $z**2)/((1 - k($r,$z,$a)**2)*Q($r,$z,$a)) * $E) / (2 * $pi * sqrt(Q($r,$z,$a)))
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
