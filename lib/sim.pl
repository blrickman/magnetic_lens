#!/usr/bin/perl

use warnings;
use strict;
use Physics::ElectronProp::Simulation_Cart;
use Physics::ElectronProp::Auxiliary ':constants';
use PDL::Util 'export2d';
use PDL;
#use PDL::Graphics::Prima::Simple [1000,500];
use File::chdir;
use Getopt::Long;
use File::Path 'make_path';
use File::Spec;
use File::Copy 'copy';

GetOptions(
  'export!'    		=> \(my $export = 1),
  'directory=s'  	=> \my $dir,
  'force'       	=> \my $force,
  'help'		=> \my $help,
);

my $sim_file = shift;
if ($help || !$dir || !$sim_file) {
print <<END;
Usage: $0 [options] file

Options:
  -e, --export		- should export data (default: true)
    --noexport
  -f, --force        	- should overwrite existing data (default: false)
  -d, --directory	- directory in which to import data
  -h, --help		- Shows this message

END
exit 1;
}

if ($export) { 
  unless ($dir) {
    my $clean_file = $sim_file;
    $clean_file =~ s/\./_/g;
    $dir = File::Spec->catdir( 'dir', $clean_file );
  }
 
  if (-d $dir and !$force) {
    die "directory $dir exists, stopping\n";
  }
 
  make_path $dir;
  print "Writing data files to:\n$dir\n";
 
  if ( -e File::Spec->catfile( $dir, $sim_file ) and !$force ) {
    die "file $sim_file already exists in $dir, stopping\n";
  }
 
  copy $sim_file, $dir;
}

my $sim = do $sim_file or die "Error reading file $sim_file: $!";
$sim->{dir} = $dir;
$sim->run();

## Single electron file history

if (@{$sim->electrons} == 1) {
  my $steps = $sim->steps / 1000;
  my $focus;
  print "focus: ";
  print $focus = `gnuplot -e "dir='$dir'; scale=$steps" .plot_single-ray.gp`;
  my (@rfcE, @solI);
    for (@{$sim->lens}) {
      push @rfcE, $_->E_0 if $_->lens_type eq 'rfcavity';
      push @solI, $_->current if $_->lens_type eq 'solenoid';
    }
  my $rfcE = join('; ', @rfcE) eq '' ? "na\t" : join('; ', @rfcE);
  my $solI = join('; ', @solI) eq '' ? "na\t" : join('; ', @solI);
  my $exists = -e "$dir/fit_history.dat";
  open my $FIT, '>>' . "$dir/fit_history.dat";
  print $FIT "Current,\t E-field,\t Focus\n" unless $exists;
  print $FIT "$solI,\t $rfcE,\t $focus";
}

## Run Gnuplot

#system( "gnuplot -e 'file1 = \"..\/$dir\"' plot_single-ray.gp -");

__END__
## Export Lens Shapes
for my $lens (@{ $sim->lens }) {
  my $fn_sol = join('_',
  $lens->{sol_name},
  $lens->{front_radius} . 'm',
  $lens->{back_radius} . 'm',
  $lens->{sol_length} . 'm',
  $lens->{front_pos} . 'm',
  );
  export(transpose(cat($lens->plot_lens)),$fn_sol . '.lns', 'data/lens_shape',0);
}

my ($r,undef, $z) = dog(transpose($position));
plot(
  -height1	=> ds::Pair($z*100,$r*1000,
    color 	=> cl::Red,
    plotType 	=> ppair::Lines,
  ),
  -lens	=> ds::Pair($lensx*100,$lensy*1000,
    color 	=> cl::Black,
    plotType 	=> ppair::Lines,
  ),
  x		=> {label => 'z (cm)' },# , min => $zi, max => $zf},
  y		=> {label => 'r (mm)' , min => -.0001, max => .26},
);
