package Physics::ElectronProp::Electron;

use warnings;
use strict;
use v5.10;
use PDL;
use Physics::ElectronProp::Auxiliary ':constants';

sub new {
  my $class = shift;
  my $self = {@_};
  bless($self, $class);
  $self->_init;
  return $self;
}

sub _init {
  my $self = shift;
  $self->{charge} = -1 * qe unless defined $self->{charge};
  $self->{mass} = me unless defined $self->{mass};
  if (defined $self->energy && defined $self->velocity) {
    warn "energy and velocity given, possible conflict, using velocity parameter\n";
  } else{
    $self->velocity( $self->KE_to_vel( $self->energy * qe, $self->mass));
    print "Converting energy (" . $self->energy / 1000 ."keV) to velocity along z\nv = ". sprintf ("%.4f",$self->{velocity}[2] / vc ) . " c\n";
  }
  $self->min_rad([sqrt($self->position->[0]**2+$self->position->[1]**2),$self->position->[2]]);
  $self->small_rad($self->position->[0] * .1);
  $self->{position} = pdl $self->position;
  $self->{velocity} = pdl $self->velocity;
}

## Initial Conditions of Electron ##

sub energy   { $_[0]->{energy   }=$_[1] if defined $_[1]; $_[0]->{energy   } }
sub position { $_[0]->{position }=$_[1] if defined $_[1]; $_[0]->{position } }
sub velocity { $_[0]->{velocity }=$_[1] if defined $_[1]; $_[0]->{velocity } }

sub mass     { $_[0]->{mass     } }
sub charge   { $_[0]->{charge   } }

sub min_rad  { $_[0]->{min_rad  }=$_[1] if defined $_[1]; $_[0]->{min_rad  } }

sub small_rad  { $_[0]->{small_rad  }=$_[1] if defined $_[1]; $_[0]->{small_rad  } }
sub near_lens{ $_[0]->{near_lens}->{$_[1]} = $_[2] if defined $_[2]; $_[0]->{near_lens}->{$_[1]} }
sub previous_force { $_[0]->{previous_force} = $_[1] if defined $_[1]; $_[0]->{previous_force} }

## Storing and Conversion Functions ##

sub history {
  my $self  = shift;
  my $place = shift;
  my $store = shift;
  push @{$self->{$place}}, $store if defined $store;
  $self->{$place};
}

sub KE_to_vel { 
  my $self = shift;
  my ($energy,$mass) = @_;
  return [0,0,vc*(1-($energy/($mass*vc**2)+1)**-2)**.5];
}

1;
