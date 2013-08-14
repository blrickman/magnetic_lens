#!/usr/bin/perl

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;
use PDL::Graphics::Gnuplot;

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
my $Icur   	= .5;

## Initial Conditions of Electron ##

my $vi		= 1;
my $xi		= $aback/10;
my $yi		= 0;
my $zi		= -1 ;
my $zf		=  2 ;

my $x_c		= [$xi,$yi,$zi];
my $v_c		= [0,0,$vi];
my @x_store	= [$xi,$yi,$zi];
my @v_store	= [0,0,$vi];

#my $B = mag_acc([0,0,0],[0,0,1]); print join(', ',@$B) . "\n";

## Simulation ##

my $step 	= 10;
my $tstep 	= ($zf - $zi) / $vi /$step;

for (0..$step) {
  
  my $ao = mag_acc($x_c,$v_c);
  print join(", \t",@$x_c) . ";\t" . join(", \t",@$v_c) . "\n";
  for my $dim (0..2) {
    $$x_c[$dim] = $$x_c[$dim] + $$v_c[$dim] * $tstep + $$ao[$dim] * $tstep**2 / 2;
    $$v_c[$dim] = $$ao[$dim] * $tstep + $$v_c[$dim];
  }
  print $x_store[-2]->[0] . "\n\n";
  push @x_store, $x_c;
  push @v_store, $v_c;
}

my $x_data = pdl \@x_store;
my $v_data = pdl \@v_store;

for my $dim (0..2) {
  for my $ele (0..@x_store-1) {
print $x_store[$ele]->[$dim];
}}

my $x = $x_data->slice('(0),');
my $y = $x_data->slice('(1),');
my $z = $x_data->slice('(2),');

#gplot($x,$z);

sub mag_acc {
  my ($x,$v) = @_;
  my @i = (0,1,2,0,1);
  my @a;
  my $B = Btot($x); 
  for (0..2) {
    $a[$_] = $$v[$i[$_ + 1]] * $$B[$i[$_ + 2]] - $$v[$i[$_ + 2]] * $$B[$i[$_ + 1]];
    $a[$_] *= $charge / $mass;
  }
  return \@a;
}

## Setup of solenoid ##

sub Bloop {
  my ($x,$n) = @_;
  my @Bloop = (Br( $$x[0], $$x[2] - zstep($n), a($n)), 0, Bz( $$x[0], $$x[2] - zstep($n), a($n)));
  return \@Bloop;
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
  my @Btot;
  for my $n (0..$nloop) {
    my $Bloop = Bloop($x,$n);
    for my $dim (0..2) {
      $Btot[$dim] += $$Bloop[$dim];
    }
  }
  return \@Btot;
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
