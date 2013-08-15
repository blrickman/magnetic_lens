package Physics::ElectronProp:Simulation;

use warnings;
use strict;
use v5.10;
use PDL;


my $xi		= 10**-2/40	; # m
my $yi		= 0		; # m
my $zi		= -5*10**-2	; # m
my $zf		=  5*10**-2	; # m
my $step 	= 500;
my $tstep 	= ($zf - $zi) / $vi /$step; # s


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
