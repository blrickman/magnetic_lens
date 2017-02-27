package Physics::ElectronProp::Mesh_Lens;
use warnings;
use strict;
use Physics::ElectronProp::Auxiliary ':constants';
use parent 'Physics::ElectronProp::EM_lens';
use Tie::File;
use Fcntl 'O_RDONLY', 'O_RDWR';
use PDL;


sub _init {
  my $self = shift;
  my @field_files = keys $self->{mesh_file};
  for (@field_files) {
    my $file = $self->mesh_file($_);
    tie my @data, 'Tie::File', $file, mode => O_RDONLY or die "Can't open " . $self->mesh_file($_) . ": $!";
    $self->{mesh_comments}{$_} = `grep -c '#' $file`; #skip commented lines
    $self->{mesh_array}{$_} = \@data;
  }
  my $fn = (split '/', $self->mesh_file($field_files[0]))[-1];
  my (undef,$zi,$zf,$zs,$ri,$rf,$rs,$zlpos) = split '_', $fn;
  $zlpos =~ s/\.\w+$//;
  $self->{mesh_length} = $zf - $zi;
  $self->{mesh_radius} = $rf - $ri;
  $self->{mesh_step  } = {'z' => $zs, 'r' => $rs};
  $self->{mesh_start } = {'z' => $zi, 'r' => $ri};
  $self->{phase} = {'B',0,'E',0} unless defined $self->{phase};
  $self->{omega} = 0 unless defined $self->{omega};
  $self->{mesh_ref_pt} = $zlpos - $zi if defined $zlpos;
  unless (ref $self->{mesh_constant} eq ref {}) {
    my $const = $self->{mesh_constant};
    $self->{mesh_constant} = { 'E' => $const, 'B' => $const,};
    warn "Extending mesh constant to E and B\n";
  }
  die unless defined $self->mesh_ref_pt;
  if ($self->mesh_ref_pt > $zf-$zi || $self->mesh_ref_pt < 0) {
    die "Check the reference point within the mesh :$!";
  }
  $self->{front_pos} = $self->front_pos - $self->mesh_ref_pt;
  $self->{mesh_pos} = $self->front_pos + $self->mesh_ref_pt;
}

## Mesh Lens Set Parameters ##

sub phase 	{ $_[0]->{phase 	}{$_[1]} }
sub omega	{ $_[0]->{omega		} }
sub mesh_constant { $_[0]->{mesh_constant}{$_[1]} }
sub mesh_file   { $_[0]->{mesh_file	}{$_[1]} }	

## Mesh Lens Auto Parameters ##

sub lens_type	{ 'mesh' 	}
sub mesh_ref_pt { $_[0]->{mesh_ref_pt } }		# Mesh Space
sub mesh_step 	{ $_[0]->{mesh_step 	}{$_[1]} }	# Mesh (and Sim) Space
sub mesh_start 	{ $_[0]->{mesh_start 	}{$_[1]} }	# Mesh Space
sub mesh_radius	{ $_[0]->{mesh_radius	} }		# Mesh Space
sub mesh_length	{ $_[0]->{mesh_length	} }		# Mesh Space
sub mesh_pos	{ $_[0]->{mesh_pos	} }		# Sim Space
sub mesh_array	{ $_[0]->{mesh_array	}{$_[1]}[$_[2] + $_[0]->{mesh_comments}{$_[1]} ] } # skip commented lines


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
  my $const = $self->mesh_constant($field);
  my ($MR,$MZ,$MZ0) = ($self->mesh_radius, $self->mesh_length,$self->front_pos);
  if ($r <= $MR && $z <= $MZ + $MZ0 && $z >= $MZ0) {
    $z = $z - $self->mesh_pos;
    return $self->get_fields($field,($r,$z)) * cos($omega*$t - $phase) * $const;
  }
  return zeros(3)
}


# Field Retrieval from mesh

sub get_fields {
  my $self = shift;
  my ($field,$r,$z) = @_;
  my @field;
  for ( qw/ r t z /) {
    my $f_val = 0;
    if (defined($self->mesh_file($field . $_))) {
      my $pos = pdl (1,$z,$r,$z*$r);
      my $box = $self->box_zr($field . $_,$z,$r); 
      $f_val = inner($pos,$self->fit_4pt($box));
    }
    push @field, $f_val;
  }
  #print pdl( \@field) . "\n" if $field eq 'B';
  return pdl \@field;
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
  return \@box; #
}

sub fit_4pt {
  my ($self,$pts) = @_;
  # $pts is a box_zr, it should have 4 elements.
  # Each element should also have 3 subelements.
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
  die "This val doesn't exist: (" . $self->mesh_array($fn_key,140) . ")" unless $self->mesh_array($fn_key,$i);
  my @row = split ', ', $self->mesh_array($fn_key,$i);
  return $row[$j];
}

sub zr_ij {
  my $self = shift;
  my ($z, $r) = @_;
  return (
    int (($z - $self->mesh_start('z'))/$self->mesh_step('z')), 
    int (($r - $self->mesh_start('r'))/$self->mesh_step('r'))
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


