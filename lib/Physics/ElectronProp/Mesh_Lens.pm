package Physics::ElectronProp::Mesh_Lens;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use Tie::File;
use PDL;

sub _init {
  my $self = shift;
  my @field_files = keys $self->{mesh_file};
  for (@field_files) {
    tie my @data, 'Tie::File', $self->mesh_file($_) or die "Can't open " . $self->mesh_file($_) . ": $!";
    $self->{mesh_array}{$_} = \@data;
  }
  my (undef,$zi,$zf,$zs,$ri,$rf,$rs) = split '_', $self->mesh_file($field_files[0]);
  $rs =~ s/\.\w+$//;
  $self->{mesh_length} = $zf - $zi;
  $self->{mesh_radius} = ($rf - $ri);#/2;
  $self->{mesh_step  } = {'z' => $zs, 'r' => $rs};
  $self->{mesh_start } = {'z' => $zi, 'r' => $ri};
  $self->{phase} = {'B',0,'E',0} unless defined $self->{phase};
  $self->{omega} = 0 unless defined $self->{omega};
  if ($self->mesh_lens_pos > $zf || $self->mesh_lens_pos < $zi) {
    die "Check the lens position within the mesh :$!";
  }
  $self->{mesh_pos} = $self->front_pos - $self->mesh_lens_pos;
}

## Mesh Lens Extra Parameters ##

sub lens_type	{ 'mesh' 	}
sub mesh_file   { $_[0]->{mesh_file	}{$_[1]} }
sub mesh_step 	{ $_[0]->{mesh_step 	}{$_[1]} }
sub mesh_start 	{ $_[0]->{mesh_start 	}{$_[1]} }
sub mesh_radius	{ $_[0]->{mesh_radius	} }
sub mesh_length	{ $_[0]->{mesh_length	} }
sub mesh_lens_pos { $_[0]->{mesh_lens_pos } }
sub mesh_pos	{ $_[0]->{mesh_pos	} }
sub mesh_array	{ $_[0]->{mesh_array	}{$_[1]}[$_[2]] }
sub phase 	{ $_[0]->{phase 	}{$_[1]} }
sub omega	{ $_[0]->{omega		} }

## E and B Fields ##

sub B_Field {
  my ($self,$pos) = @_;
  $self->Field($pos,'B')
}

sub E_Field {
  my ($self,$pos) = @_;
  $self->Field($pos,'E')
}

sub Field {
  my ($self,$pos,$field) = @_;
  my ($r,$theta,$z,$t) = list $pos;
  my $omega = $self->omega;
  my $phase = $self->phase($field);
  if ($r <= $self->mesh_radius && $z <= $self->mesh_length + $self->mesh_pos && $z >= $self->mesh_pos) {
    return $self->get_fields($field,($r,$z-$self->mesh_pos)) * cos($omega*$t - $phase);
  }
  return zeros(3)
}


# Field Retrieval from mesh

sub get_fields {
  my $self = shift;
  my ($field,$r,$z) = @_;
  my @field;
  for ('z','r') {
    my $f_val = 0;
    if (defined($self->mesh_file($field . $_))) {
      my $pos = pdl (1,$z,$r,$z*$r);
      my $box = $self->box_zr($field . $_,$z,$r);
      $f_val = inner($pos,$self->fit_4pt($box));
    }
    push @field, $f_val;
  }
  return pdl ($field[1], 0, $field[0]);
}

sub box_zr {
  my $self = shift;
  my ($fn_key,$z,$r) = @_;
  my ($i,$j) = $self->zr_ij(($z,$r));
  my @box;
  $i = $i == $self->mesh_length / $self->mesh_step('z') ? $i-1 : $i;
  $j = $j == $self->mesh_radius / $self->mesh_step('r') ? $j-1 : $j;
  for my $ii ($i,$i+1) {
    for my $jj ($j,$j+1) {
      push @box, [$self->ij_zr($ii,$jj),$self->ij_element($fn_key,$ii,$jj)];
    }
  }
  return \@box;
}

sub fit_4pt {
  my ($self,$pts) = @_;
  my ($x0,$y0,$z00) = @{$$pts[0]};
  my ($x1,$y1,$z11) = @{$$pts[3]};
  my ($z01,$z10) = ($$pts[1][2],$$pts[2][2]);
  my $c0 = ($x1*$y1*$z00 - $x1*$y0*$z01 - $x0*$y1*$z10 + $x0*$y0*$z11)/(($x0 - $x1)*($y0 - $y1)); 
  my $c1 = ($y1*(-$z00 + $z10) + $y0*($z01 - $z11))/(($x0 - $x1)*($y0 - $y1)); 
  my $c2 = ($x1*(-$z00 + $z01) + $x0*($z10 - $z11))/(($x0 - $x1)*($y0 - $y1));
  my $c3 = ($z00 - $z01 - $z10 + $z11)/(($x0 - $x1)*($y0 - $y1));
  pdl ($c0,$c1,$c2,$c3);
}

sub ij_element {
  my $self = shift;
  my ($fn_key,$i,$j) = @_;
  my @row = split ', ', $self->mesh_array($fn_key,$i);
  return $row[$j];
}

sub zr_ij {
  my $self = shift;
  my ($z, $r) = @_;
  return (
    int ($z - $self->mesh_start('z'))/$self->mesh_step('z'), 
    int ($r - $self->mesh_start('r'))/$self->mesh_step('r')
  );
}

sub ij_zr {
  my $self = shift;
  my ($i, $j) = @_;
  return (
    $self->mesh_start('z') + $i*$self->mesh_step('z'), 
    $self->mesh_start('r') + $j*$self->mesh_step('r')
  );
}


1;
__END__


