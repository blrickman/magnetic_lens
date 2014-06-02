package Physics::ElectronProp::RF_Cavity;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use PDL;
use YAML 'LoadFile';

sub _init {
  my $self = shift;
  $self->{front_radius}=$self->radius;
  $self->{back_radius}=$self->radius;
  if ($self->{mode}{field} eq 'TM') {
    $self->rootx(LoadFile('J_m10-n10.dat'));
    $self->{Cavity_Field} = \&TM_Fields
  } elsif ($self->{mode}{field} eq 'TE') {
    $self->rootx(LoadFile('Jprime_m10-n10.dat'));
    $self->{Cavity_Field} = \&TE_Fields
  } else {
    die "Cavity type is undefined: $!";
  }
}

## RF Cavity Lens Extra Parameters ##

sub lens_type	{ 'rfcavity' 	}
sub E_0		{ $_[0]->{E_0		} }
sub mode	{ $_[0]->{mode		} }
sub epsilon	{ $_[0]->{epsilon	} }
sub mu		{ $_[0]->{mu		} }
sub radius	{ $_[0]->{radius	} }
sub omega	{ $_[0]->{omega		} }
sub phase	{ $_[0]->{phase		} }

sub Cavity_Field{ $_[0]->{Cavity_Field		} }
sub gamma_2	{ $_[0]->mu * $_[0]->epsilon * $_[0]->omega**2 - ($_[0]->{mode}{p} * pi / $_[0]->lens_length )**2 }
sub gamma_mn	{ $_[0]->{rootx}->[$_[0]->{mode}{m}][$_[0]->{mode}{n}] / $_[0]->radius }
sub omega_mnp	{ sqrt(($_[0]->gamma_mn**2 + ($_[0]->{mode}{p} * pi / $_[0]->lens_length )**2) / ($_[0]->mu * $_[0]->epsilon)) }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  if ($r <= $self->radius && $z <= $self->lens_length + $self->front_pos && $z >= $self->front_pos) {
    return $self->Fields('B',($r,$theta,$z-$self->front_pos,$t))
  }
  return zeros(3)
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
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

  my ($R, $d, $E0, $omega, $epsilon, $mu, $gamma2, $gamma_mn, $phi) = map { $self->$_ } qw/radius lens_length E_0 omega epsilon mu gamma_2 gamma_mn phase/;
  my ($m, $n, $p) = map { $self->{mode}{$_} } qw/m n p/;

  if ($field eq 'E') {
    my $Er = $p==0 ? 0 : -$E0 * $gamma_mn * $p * pi / (4 * $d * $gamma2) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * sin(2*$p*pi*$z/$d) * cos($omega * $t - $phi);
    my $Et = $m==0 || $p ==0 ? 0 : $E0 * $m * $p * pi / (2 * $d * $gamma2) * (bessjn($gamma_mn*$r,$m) / $r ) * sin(2*$p*pi*$z/$d) * sin($omega * $t - $phi);
    my $Ez = $E0 * bessjn($gamma_mn*$r,$m) * cos($p*pi*$z/$d) * cos($omega * $t - $phi);
    return pdl [$Er,$Et,$Ez]
  } else {
    my $Br = $m==0 ? 0 : $mu * $E0 * $m * $epsilon * $omega / ($gamma2) * (bessjn($gamma_mn*$r,$m) / $r ) * cos($p*pi*$z/$d)**2 * cos($omega * $t - $phi);
    my $Bt = $mu * $E0 * $gamma_mn * $epsilon * $omega / (2 * $gamma2) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * cos($p*pi*$z/$d)**2 * sin($omega * $t - $phi);
    my $Bz = 0;
    return pdl [$Br,$Bt,$Bz]
  }
}

sub TE_Fields {
  my $self = shift;
  my $field  = shift;
  my ($r,$theta,$z,$t) = @_;

  my ($R, $d, $E0, $omega, $epsilon, $mu, $gamma2, $gamma_mn, $phi) = map { $self->$_ } qw/radius lens_length E_0 omega epsilon mu gamma_2 gamma_mn phase/;
  my ($m, $n, $p) = map { $self->{mode}{$_} } qw/m n p/;

  if ($field eq 'E') {
    my $Er = $m==0 ? 0 : sqrt($epsilon*$mu) * $E0 * $m * $omega / ($gamma2) * (bessjn($gamma_mn*$r,$m) / $r ) * sin($p*pi*$z/$d)**2 * cos($omega * $t - $phi);
    my $Et = - sqrt($epsilon*$mu) * $E0 * $gamma_mn * $omega / (2 * $gamma2) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * sin($p*pi*$z/$d)**2 * sin($omega * $t - $phi);
    my $Ez = 0;
    return pdl [$Er,$Et,$Ez]
  } else {
    my $Br = $p==0 ? 0 : sqrt($epsilon/$mu) * $E0 * $gamma_mn * $p * pi / (4 * $d * $gamma2) * (bessjn($gamma_mn*$r,$m-1) - bessjn($gamma_mn*$r,$m+1)) * sin(2*$p*pi*$z/$d) * cos($omega * $t - $phi);
    my $Bt = $m==0 || $p ==0 ? 0 : sqrt($epsilon/$mu) * $E0 * $m * $p * pi / (2 * $d * $gamma2) * (bessjn($gamma_mn*$r,$m) / $r ) * sin(2*$p*pi*$z/$d) * sin($omega * $t - $phi);
    my $Bz = sqrt($epsilon/$mu) * $E0 * bessjn($gamma_mn*$r,$m) / $r * sin($p*pi*$z/$d) * cos($omega * $t - $phi);
    return pdl [$Br,$Bt,$Bz]
  }
}
1;
__END__

my %fields = {
  TM => \&TM_Fields,
  TE => \&TE_Fields,
}


