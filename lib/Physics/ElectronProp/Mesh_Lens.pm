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
  my (undef,$zi,$zf,$zs,$yi,$yf,$ys) = split '_', $self->mesh_file($field_files[0]);
  $ys =~ s/\.\w+$//;
  $self->{mesh_length} = $zf - $zi;
  $self->{mesh_radius} = $yf - $yi;
  $self->{mesh_step  } = {'z' => $zs, 'y' => $ys};
  $self->{phase} = {'B',0,'E',0} unless defined $self->{phase};
  $self->{omega} = 0 unless defined $self->{omega};
}

## Mesh Lens Extra Parameters ##

sub lens_type	{ 'mesh' 	}
sub mesh_file   { $_[0]->{mesh_file	}{$_[1]} }
sub mesh_step 	{ $_[0]->{mesh_step 	}{$_[1]} }
sub mesh_radius	{ $_[0]->{mesh_radius	} }
sub mesh_length	{ $_[0]->{mesh_length	} }
sub mesh_lineup { $_[0]->{mesh_lineup	} }
sub mesh_array	{ $_[0]->{mesh_array	}{$_[1]}[$_[2]] }
sub phase 	{ $_[0]->{phase 	}{$_[1]} }
sub omega	{ $_[0]->{omega		} }

## E and B Fields ##

sub B_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  my $omega = $self->omega;
  my $phase = $self->phase(B);
  if ($r <= $self->mesh_radius && $z <= $self->mesh_length + $self->front_pos && $z >= $self->front_pos) {
    return $self->get_fields('B',($r,$z-$self->front_pos)) * cos($omega*$t - $phi);
  }
  return zeros(3)
}

sub E_Field {
  my $self = shift;
  my $pos  = shift;
  my ($r,$theta,$z,$t) = list $pos;
  my $omega = $self->omega;
  my $phase = $self->phase(E);
  if ($r <= $self->radius && $z <= $self->lens_length + $self->front_pos && $z >= $self->front_pos) {
    return $self->get_fields('E',($r,$z-$self->front_pos)) * cos($omega*$t - $phi);
  }
  return zeros(3)
}

# Field Retrieval from mesh

sub get_efields {
  my $self = shift;
  my ($field,$r,$z) = @_;
  
  my @field_2d;
  for my $dir (0..1) {
    my ($ez, $er) = ( arr_grab($pos->index(0),\@z), arr_grab($pos->index(1),\@r) );
    my $line_test = $z + ($r * ($$ey[0]-$$ey[1]) + ($$ez[1]*$$ey[1]-$$ez[0]*$$ey[0])) / ($$ez[0] - $$ez[1]);
    my @points;
    for my $i (0..1) {
      $points[$i] = pdl ($$ez[!$i], $$ey[$i],$efield{$$ez[!$i]}{$$ey[$i]}[$dir]);
    }
    $points[2] = pdl ($$ez[0], $$ey[0],$efield{$$ez[0]}{$$ey[0]}[$dir]);
    if (sprintf "%.14f", $line_test > 0) {
      $points[2] = pdl ($$ez[1], $$ey[1],$efield{$$ez[1]}{$$ey[1]}[$dir]);
    }
    my $plane = crossp($points[0]-$points[2],$points[1]-$points[2]);
    $plane /= $plane->index(2);
    push @Ef, inner($plane,$points[2]) - inner($plane,pdl ($z,$r,0));
  }
  return pdl @Ef
}

sub arr_grab {
  my $self = shift;
  my ($pos, $array) = @_;
  my $step = $$array[1] - $$array[0];
  my @return = grep $_ >= $pos - $step, grep $_ <= $pos + $step, @$array;
  push @return, $return[0] + $step if @return == 1;
  shift @return if @return == 3;
  die "I grabbed " . @return . " elemtents when I only meant to grab 2: $!" if @return != 2;
  return \@return;
}

sub zy_pot {
  return ij_pot(zy_ij(@_));
}

sub ij_pot {
  my ($i,$j) = @_;
  my @row = split ', ', $datafile[$i];
  return $row[$j];
}
sub zy_ij {
  my ($z, $y) = @_;
  return (int ($z-$zi)/$zs, int ($y-$yi)/$ys);
}



1;
