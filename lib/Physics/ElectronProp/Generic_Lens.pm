package Physics::ElectronProp::Generic_Lens;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use PDL;

sub _init {

}

## RF Cavity Lens Extra Parameters ##

sub lens_type	{ 'generic' 	}
sub coeff       { $_[0]->{coeff		} }
sub radius	{ $_[0]->{radius	} }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  return zeros(3)
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  if ($r <= $self->radius && $z <= $self->lens_length + $self->front_pos && $z >= $self->front_pos) {
    my $Er = 0;
    my @el = @{$self->{coeff}};
    for (0..@el-1) {
      $Er += $r**$_ * $el[$_];
    } 
    return pdl [$Er,0,0]
  }
  return zeros(3)
}

1;
