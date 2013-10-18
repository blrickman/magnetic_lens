package Physics::ElectronProp::EM_lens;

use warnings;
use strict;
use PDL;
use PDL::GSLSF::ELLINT;
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
  $self->B_Field(pdl [0,0,0]);
  $self->E_Field(pdl [0,0,0]);
}

## Export Lens Diagram ##

sub plot_lens {
  my $self = shift;
  my $length = sequence(50)/49;
  my $height = $length;
  my $f_l = $height * $self->front_radius;
  my $t_l = $self->lens_shape->($length);
  my $b_l = $height->slice('-1:0') * $self->back_radius;
  my $lens_x = pdl (list($self->front_pos * ones(50)),list($length * $self->lens_length + $self->front_pos),list(($self->front_pos + $self->lens_length) * ones(50)));
  my $lens_y = pdl (list($f_l), list($t_l), list($b_l));
  return ($lens_x, $lens_y);
}

## Lens Parameters ##

sub lens_name  		{ $_[0]->{lens_name     } }
sub lens_shape  	{ $_[0]->{lens_shape    } }
sub front_radius   	{ $_[0]->{front_radius  } }
sub back_radius    	{ $_[0]->{back_radius   } }
sub lens_length     	{ $_[0]->{lens_length   } }
sub front_pos	   	{ $_[0]->{front_pos     } }

## E and B Fields to be defined by child classes ##

1;
