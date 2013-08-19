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
    print "Converting energy (" . $self->energy."eV) to velocity along z\nv = ". $self->{velocity}[2] / vc . " c\n";
  }

  $self->{position} = pdl $self->position;
  $self->{velocity} = pdl $self->velocity;
}

## Initial Conditions of Electron ##

sub energy   { $_[0]->{energy   }=$_[1] if defined $_[1]; $_[0]->{energy   } }
sub position { $_[0]->{position }=$_[1] if defined $_[1]; $_[0]->{position } }
sub velocity { $_[0]->{velocity }=$_[1] if defined $_[1]; $_[0]->{velocity } }

sub mass     { $_[0]->{mass     }=$_[1] if defined $_[1]; $_[0]->{mass     } }
sub charge   { $_[0]->{charge   }=$_[1] if defined $_[1]; $_[0]->{charge   } }

## Storing and Conversion Functions ##

sub history {
  my $self  = shift;
  my $place = shift;
  my $store = shift;
  $store = [list($store)] if defined $store;
  for (0..2) {
    push @{$self->{$place}[$_]}, $$store[$_] if defined $store;
  }
  $self->{$place};
}

sub KE_to_vel { 
  my $self = shift;
  my ($energy,$mass) = @_;
  return [0,0,vc*(1-($energy/($mass*vc**2)+1)**-2)**.5];
}

1;
