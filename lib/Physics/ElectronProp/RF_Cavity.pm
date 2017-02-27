package Physics::ElectronProp::RF_Cavity;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use PDL;
use YAML 'LoadFile';
use Math::Polynomial::Solve qw(:classical); 

sub _init {
  my $self = shift;
  $self->{front_radius}=$self->radius;
  $self->{back_radius}=$self->radius;
  if ($self->{cav_mode}{field} eq 'TM') {
    $self->rootx(LoadFile('J_m10-n10.dat'));
    $self->{Cavity_Field} = \&TM_Fields
  } elsif ($self->{cav_mode}{field} eq 'TE') {
    $self->rootx(LoadFile('Jprime_m10-n10.dat'));
    $self->{Cavity_Field} = \&TE_Fields
  } else {
    die "Cavity type is undefined: $!";
  }
  print $self->name . " has a frequency of " . sprintf('%.3e',$self->omega/(2*pi)) . "Hz.\n";
  $self->{off}{'B'} = 0 unless defined $self->{off}{'B'};
  $self->{off}{'E'} = 0 unless defined $self->{off}{'E'};
}

## RF Cavity Lens Extra Parameters ##

sub lens_type	{ 'rfcavity' 	}
sub name	{ $_[0]->{name		} }
sub E_0		{ $_[0]->{E_0		} }
sub cav_mode	{ $_[0]->{cav_mode	} }
sub epsilon	{ $_[0]->{epsilon	} }
sub mu		{ $_[0]->{mu		} }
sub radius	{ $_[0]->{radius	} }
sub phase	{ $_[0]->{phase		} }
sub test	{ $_[0]->{test		} }
sub off		{ $_[0]->{off		}{$_[1]} }

sub Cavity_Field{ $_[0]->{Cavity_Field		} }

sub gamma_mn	{ $_[0]->{rootx}->[$_[0]->{cav_mode}{m}][$_[0]->{cav_mode}{n}] / $_[0]->radius }
sub omega	{ sqrt(($_[0]->gamma_mn**2 + ($_[0]->{cav_mode}{p} * pi / $_[0]->lens_length )**2) / ($_[0]->mu * $_[0]->epsilon)) }

## Methods ##

sub phase_shift {
  my ($self,$electron, $type) = @_;
  my ($E0, $L, $omega) = map { $self->$_ } qw/E_0 lens_length omega/;
  my ($mass, $charge) = map { $electron->$_ } qw/mass charge/;
  my ($z0,$v0) = map { $electron->$_->index(2) } qw /position velocity/;
  $charge *= -1;

  my $shift = $self->omega*($self->front_pos+$L/2-$z0)/$v0;
  if ($type ==2) {
    my  @x3 = (cubic_roots(sprintf("%.16f", -2*$E0*$charge/(3*$mass*$omega**2)),0,sprintf("%.16f", 2*$v0/$omega),sprintf("%.16f", -$L)));
    my $lowroot = pi;
    for my $root (@x3) {
      next if $root =~ /i/;
      #next if $root < 0;
      $lowroot = $root < $lowroot ? $root : $lowroot;
    }
    $shift = $lowroot + $self->omega*$self->front_pos/$v0 if $type == 2;
  } elsif ($type == 3) {
    my $C1 = (-3 * $E0**2 * $L * $mass * $charge**2 * $omega**2 + sqrt($E0**3 * $mass**2 * $charge**3 * $omega**3 * (16 * $mass * $v0**3 + 9 * $E0 * $L**2 * $charge * $omega)))**(1/3)/2**(2/3);
    $shift = $mass*$v0*$omega/$C1 - $C1/$E0/$charge + $self->omega*$self->front_pos/$v0 if $type == 3;
  }
  #printf "Phase shift is : %.5f \n", $shift;
  return $shift;
}

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  return zeros(3) if $self->off('B');
  if ($r <= $self->radius && $z <= $self->lens_length + $self->front_pos && $z >= $self->front_pos) {
    return $self->Fields('B',($r,$theta,$z-$self->front_pos,$t))
  }
  return zeros(3)
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  return zeros(3) if $self->off('E');
  if ($r <= $self->radius && $z <= $self->lens_length + $self->front_pos && $z >= $self->front_pos) {
    return $self->Fields('E',($r,$theta,$z-$self->front_pos,$t))
  }
  return zeros(3)
}

## Field Functions ##

sub rootx { $_[0]->{rootx }=$_[1] if defined $_[1]; $_[0]->{rootx} }

sub Fields { 
  my $self = shift;
  my $cfield = $self->Cavity_Field;
  return $self->$cfield(@_);  
}

sub TM_Fields {
  my $self = shift;
  my $field  = shift;
  my ($r,$theta,$z,$t) = @_;

  my ($R, $d, $E0, $omega, $epsilon, $mu, $gamma_mn, $phi) = map { $self->$_ } qw/radius lens_length E_0 omega epsilon mu gamma_mn phase/;
  ($epsilon, $mu) = (epsilon_0,mu_0);
  my ($m, $n, $p) = map { $self->{cav_mode}{$_} } qw/m n p/;
  die "Cannot have a n=0 TM mode. $!" if $n == 0;
  #$omega = $self->test ? 1 : $omega;
  if ($field eq 'E') {
    my $Er = $p==0 ? 0 : -$E0 * $p * pi / (2 * $d * $gamma_mn) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * sin($p*pi*$z/$d) * cos($m * $theta);
    my $Et = $m==0 || $p ==0 ? 0 : $E0 * $m * $p * pi / ($d * $gamma_mn**2) * (bessjn($gamma_mn*$r,$m) / $r ) * sin($p*pi*$z/$d) * sin($m * $theta);
    my $Ez = $E0 * bessjn($gamma_mn*$r,$m) * cos($p*pi*$z/$d) * cos($m * $theta);
    my $E = pdl [$Er,$Et,$Ez];
    return $E if $self->test;
    return $E * pdl [cos($omega * $t - $phi),cos($omega * $t - $phi),cos($omega * $t - $phi)];
  } else {
    my $Br = $m==0 ? 0 : $E0 * $m * $mu * $epsilon * $omega / ($gamma_mn**2) * (bessjn($gamma_mn*$r,$m) / $r ) * cos($p*pi*$z/$d) * sin($m * $theta);
    my $Bt = $E0 * $mu *  $epsilon * $omega / (2 * $gamma_mn) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * cos($p*pi*$z/$d) * cos($m * $theta);
    my $Bz = 0;
    my $B = pdl [$Br,$Bt,$Bz];
    return $B if $self->test;
    return $B * pdl [sin($omega * $t - $phi),sin($omega * $t - $phi),sin($omega * $t - $phi)];
  }
}

sub TE_Fields {
  my $self = shift;
  my $field  = shift;
  my ($r,$theta,$z,$t) = @_;

  my ($R, $d, $E0, $omega, $epsilon, $mu, $gamma_mn, $phi) = map { $self->$_ } qw/radius lens_length E_0 omega epsilon mu gamma_mn phase/;
  ($epsilon, $mu) = (epsilon_0,mu_0);
  my ($m, $n, $p) = map { $self->{cav_mode}{$_} } qw/m n p/;
  die "Cannot have a n=0 or p=0 TE mode. $!" if $n == 0 || $p == 0;
  #$omega = $self->test ? 1 : $omega;
  if ($field eq 'E') {
    my $Er = $m==0 ? 0 : -$E0 * sqrt($epsilon*$mu) * $m * $omega / ($gamma_mn**2) * (bessjn($gamma_mn*$r,$m) / $r ) * sin($p*pi*$z/$d) * sin($m * $theta);
    my $Et = -$E0 * sqrt($epsilon*$mu) * $omega / (2 * $gamma_mn) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * sin($p*pi*$z/$d) * cos($m * $theta);
    my $Ez = 0;
    my $E = pdl [$Er,$Et,$Ez];
    return $E if $self->test;
    return $E * pdl [sin($omega * $t - $phi),sin($omega * $t - $phi),sin($omega * $t - $phi)];
  } else {
    my $Br = $E0 * sqrt($epsilon*$mu) * $p * pi / (2 * $d * $gamma_mn) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * cos($p*pi*$z/$d) * cos($m * $theta);
    my $Bt = $m==0 ? 0 : -$E0 * sqrt($epsilon*$mu) * $m * $p * pi / ($d * $gamma_mn**2) * (bessjn($gamma_mn*$r,$m) / $r ) * cos($p*pi*$z/$d) * sin($m * $theta);
    my $Bz = $E0 * sqrt($epsilon*$mu) * bessjn($gamma_mn*$r,$m) * sin($p*pi*$z/$d) * cos($m * $theta);
    my $B = pdl [$Br,$Bt,$Bz];
    return $B if $self->test;
    return $B * pdl [cos($omega * $t - $phi),cos($omega * $t - $phi),cos($omega * $t - $phi)];
  }
}
1;
__END__

my %fields = {
  TM => \&TM_Fields,
  TE => \&TE_Fields,
}


