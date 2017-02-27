#!/usr/bin/perl

use warnings;
use strict;
use File::chdir;
use File::Spec;
use File::Path 'make_path';
use Getopt::Long;
use File::Copy qw(copy);
use PDL;
use PDL::Fit::Polynomial;
use PDL::Fit::LM;

GetOptions(
  'files' 	=> \my $files,
  'axes=s'  	=> \(my $plot_axes = 'r(z)'),
  'scale=s'    	=> \my $scale,
  'output=s'	=> \my $output,
  'directory=s'	=> \my $subdir,
  'help'	=> \my $help,
  'minimum:i'	=> \my $find_min,
  'yrange=f{2}' => \my @yrange,
  'xrange=f{2}' => \my @xrange,
  'emfields'	=> \my $emfields,
  'gpexport:s'	=> \my $export,
  'gaussianfit' => \my $gaussian,
  'largeplot'	=> \my $largeplot,
  'linearfit:10'	=> \my $linearfit,
);

my $dir = shift;

my $HELP_MSG = <<END;
Usage: $0 [options] directory

Options:
  -f, --files		- triggers querry to select files in directory
  -a, --axes        	- select components to plot from {r,x,y,z,t,vr,vx,vy,vz}.  	
			  Ex: -a 'r(z)' (the default)
  -s, --scale		- scale the data to the ratio SCALE_X x SCALE_Y
			  Ex: -s -3x-3 (the default, scales by 10**-3 on each axis)
  -m, --minimum		- triggers the minimum find action
  -l, --linearfit	- triggers a linear fit of the back end
  -x, --xrange 		- set the x-range
			  Ex: -x (0,1)
  -y, --yrange 		- set the origin of the y-range
			  Ex: -y 0
  -e, --emfields	- include field files in electron file querry
  -g, --gpexport	- write out the plot file to the directory where the electron files are
			  Ex: -g "foo.gp"  OR  -g (default filename is plot.gp)
  -o, --output		- set output PNG filename
			  Ex: -o filename.png
  -d, --directory 	- Choose a subdirectory in pictures/
  -h, --help		- Shows this message

END

my $check_dir = opendir(my $DIR, $dir) or die "Check directory path: $dir\n\n$HELP_MSG $!";  # Makes sure $dir exists
my $HOME_DIR = $CWD; #Save home dir for later use

# Print Help Message
if ($help || !$check_dir) {
print $HELP_MSG;
exit 1;
}

my ($keyin,$plotx) = qw/inside 750/;
($keyin,$plotx) = qw/outside 1500/ if $largeplot;
my $file_query;

# Grab all electron files from DIR
my @dir_files = readdir($DIR);
my @files = grep {/_/} grep {/^e.*\.dat/} grep {$_ !~ /.+fields.dat/} @dir_files;
if ( 0 < @files ) {
  @files = sort {($a =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0] <=> ($b =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0]} grep($_=~ /e0.+\d\.dat/ , @dir_files);
  @files = (@files, sort {($a =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0] <=> ($b =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0]} grep($_=~ /e\+.+\d\.dat/ , @dir_files) );
  @files = (@files, sort {($a =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0] <=> ($b =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0]} grep($_=~ /e\-.+\d\.dat/ , @dir_files) );
} else {
  @files = sort {($a =~ /e(.+)\.dat/)[0] <=> ($b =~ /e(.+)\.dat/)[0]} grep($_=~ /e.+\d\.dat/ , @dir_files);
}
@files = sort grep($_=~ /.+fields.dat/ , @dir_files) if $emfields;
closedir($DIR);
my @plot_files = @files;

# If --emfields is triggered, --files should default to FALSE

$files = $emfields ? 1 : $files;

# If --files is flagged, keep relavent electron files
if ($files) {
  if ($files eq 'TODO') {
    @plot_files = $files[0]; # possibly add short cut to querry to only plot e10?
  } else {
    # Querries for a list of numbers (by assigned number) or shortcut 
    # and plot all using 'a', which is the default behavior.
    print "Enter numbers of files to be plotted separated by spaces or enter 'a' to plot all files:\n";
    $file_query = list_files();
    @plot_files = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
  }
}

# If --minimum is flagged, find the minimum of one electron path

if (defined $find_min) {
  if ($find_min eq 'TODO') {
    my @min_files = $files[0]; # possibly add short cut to querry to only plot e10?
  } elsif (@plot_files == 1 && !$emfields) {
    get_minimum($plot_files[0],$find_min); #fits the same path selected by the fit;
  } else {
    # Select an electron path
    print "Enter the number of the file to be fitted:\n";
    $file_query = list_files();
    my @min_files = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
    get_minimum($_,$find_min) for @min_files;
  }
}

# If --linearfit is flagged, fit linear
my @fitlabels;
if ($linearfit) {
  my @lin_fits;
  if (0) {#$linearfit eq 'TODO') {
    @lin_fits = $files[0]; # possibly add short cut to querry to only plot e10?
  } elsif (@plot_files == 1 && !$emfields) {
    @lin_fits = $plot_files[0]; #fits the same path selected by the fit;
  } else {
    # Select an electron path
    #print "Enter the number of the file to be fitted:\n";
    #list_files();
    $file_query = list_files();
    @lin_fits = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
  }
  push @fitlabels, linear_fit($_,$linearfit) for @lin_fits;
  my @elements = (0..@fitlabels-1);
  @fitlabels = map "set label \"" . $fitlabels[$_] . "\" at screen 0.2," . (.8 - $_/20) . " noenhanced", @elements;
  $fitlabels[0] = join "\n", @fitlabels;
}

# If --gaussianfit is flagged, fit gaussian w(z)
die "Specify a range in z to fit over.\n" if $gaussian && ! @xrange;
my $gauss_plot;
if ($gaussian) {
  my @gauss_fits;
  if (0) {#$linearfit eq 'TODO') {
    @gauss_fits = $files[0]; # possibly add short cut to querry to only plot e10?
  } elsif (@plot_files == 1 && !$emfields) {
    @gauss_fits = $plot_files[0]; #fits the same path selected by the fit;
  } else {
    # Select an electron path
    print "Enter the number of the file to be fitted:\n";
    $file_query = list_files();
    @gauss_fits = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
    @gauss_fits = @plot_files;
  }
  my (@w0,@z0,@zR);
  for my $i (0..@gauss_fits-1) {
    ($w0[$i],$z0[$i],$zR[$i]) = guassian_fit($gauss_fits[$i]);
  }
  my @r0 = map {($_ =~ /e[0+-]_(\d+\.+\d+)\.dat/)[0]} @gauss_fits;
  my $focus = $z0[0];
  if (@r0 > 1) {
    my $r0 = pdl @r0;
    my $z0 = pdl @z0;
    my @fit = list((fitpoly1d($r0, $z0, 3))[1]);
    $focus = $fit[0];
    printf "\nAxial focus occurs at %.3fmm\n", $focus;
  }
  for my $i (0..@z0-1) {
    my ($w0,$z0,$zR,$r0) = ($w0[$i],$z0[$i],$zR[$i],$r0[$i]);
    my $delta_r = $w0*(sqrt(1 + (($focus-$z0)/$zR)**2));#    -1);
    $gauss_plot .= sprintf ",\\\n$w0*sqrt( 1 + ((x-$z0)/$zR)**2 ) w p lt %1\$d title 'fit %1\$d'", $i+1;
    #$gauss_plot .= ",\\\n$w0*(x-$z0)/($zR**2)/sqrt( 1 + ((x-$z0)/$zR)**2 ) w l lt $fn axes x1y2 title 'fit $fn'";
    printf "%.2fmm: Min = (%.3fmm,%6.3fμm)  Δr = %6.3fμm  θ = %.3fmrad  Cs = % 9.3fm\n" , 
      $r0,
      $z0,$w0*10**3,
      $delta_r *10**3,
      abs($w0/$zR)*10**3, 
      $delta_r/abs($w0/$zR)**3 *10**-3;
    # AG model comparison
    #printf "y(21.5in): %.8fm\n", $w0 * sqrt( 1 + ((546.1-$z0)/$zR)**2 )*10**(-3);
  }
}

# Set up hash to store relavent info for each axis
my $axis = {
  'r'  => { qw% scale 10**-3  unit mm dim m  marker  (sqrt(($1)**2+($2)**2))% },
  'x'  => { qw% scale 10**-3  unit mm dim m  marker  ($1)%	},
  'y'  => { qw% scale 10**-3  unit mm dim m  marker  ($2)%	},
  'z'  => { qw% scale 10**-3  unit mm dim m  marker  ($3)%	},
  't'  => { qw% scale 10**-9  unit ns dim t  marker  ($10)%	},
  'vr' => { qw% scale 3*10**8 unit c  dim s  marker  (sqrt(($4)**2+($5)**2))% },
  'vx' => { qw% scale 3*10**8 unit c  dim m/s marker ($4)%	},
  'vy' => { qw% scale 3*10**8 unit c  dim m/s marker ($5)%	},
  'vz' => { qw% scale 3*10**8 unit c  dim m/s marker ($6)%	},
  'ax' => { qw% scale 1 unit m/s/s  dim m/s/s marker ($7)%	},
  'ay' => { qw% scale 1 unit m/s/s  dim m/s/s marker ($8)%	},
  'az' => { qw% scale 1 unit m/s/s  dim m/s/s marker ($9)%	},
  'ar' => { qw% scale 1 unit m/s/s  dim m/s/s marker (sqrt(($7)**2+($8)**2))% },
  'th' => { ( qw% scale 1 unit deg    dim deg   marker%, '(atan2(($2),($1)))' ) },
  'vth'=> { ( qw% scale 1 unit deg    dim deg   marker%, '(atan2(($5),($4)))' ) },
  'ath'=> { ( qw% scale 1 unit deg    dim deg   marker%, '(atan2(($8),($7)))' ) },
  'br' => { qw% scale 1 unit T   dim T   marker ($7)%	},
  'bth'=> { qw% scale 1 unit T   dim T   marker ($8)%	},
  'bz' => { qw% scale 1 unit T   dim T   marker ($9)%	},
  'er' => { qw% scale 1 unit V/m dim V/m marker ($4)%	},
  'eth'=> { qw% scale 1 unit V/m dim V/m marker ($5)%	},
  'ez' => { qw% scale 1 unit V/m dim V/m marker ($6)%	},
};

# Pick out the axes chosen in the --axes option
$plot_axes =~ s/\)//;
my ($y,$x) = split '\(', $plot_axes;
for ($y,$x) {
  die "$_ is not a valid axis. See --help for valid axes.\n" unless exists $axis->{$_};
}

# Set up output filename
my @filename = split "/", $dir;
shift @filename;
my $title = join ' ', @filename;
$title =~ s/_//g;
push  @filename,"${x}_$y.png";
$output = defined $output ? $output : join '-', @filename;
$output = defined $subdir ? File::Spec->catfile($subdir,$output) : $output;

# Set up scales on axes
if (defined $scale) {
  my @scale = split 'x', $scale;
  die "Invalid scale entry." unless @scale == 2;
  for ($x,$y) {
    my $scale = shift @scale;
    $axis->{$_}{scale} = 10**$scale;
    $axis->{$_}{unit}  = "10^{$scale}" . $axis->{$_}{dim};
  }
}

# Set up axes for gplot script
@yrange = ('*','*') unless @yrange;
@xrange = ('*','*') unless @xrange;
$fitlabels[0] = '' unless @fitlabels;
my $axes = join ":", map ($axis->{$_}{marker} . "/(" . $axis->{$_}{scale} .")", ($x,$y));
my $plot_files = join(",\\\n",map("'$dir${_}' u $axes w l lw 2 title '$_'",@plot_files));
$plot_files = $plot_files . $gauss_plot if defined $gauss_plot;

# Set up labels for gplot script
my @label = map {"$_ (" . $axis->{$_}{unit} . ")"} ($x,$y);

# Set up gplot script
my $GNUPLOT_SCRIPT = <<END;
set terminal pngcairo enhanced solid font "times,18" size $plotx, 500 
set output "pictures/$output"
set title "$title"
set key noenhanced $keyin
set yrange [${yrange[0]}:${yrange[1]}]
set xrange [${xrange[0]}:${xrange[1]}]
set xlabel "${label[0]}"
set ylabel "${label[1]}"
${fitlabels[0]}
plot $plot_files
END

# Store script in a temporary file
open my $TMP, "> .plot.tmp" or die $!;
print $TMP $GNUPLOT_SCRIPT;
close $TMP;
if (defined $export) {
  $export = $export eq '' ? "plot.gp" : $export;
  copy(".plot.tmp","$dir/$export") or die "Copy failed: $!";
}

# Execute script
`gnuplot .plot.tmp`;

# Open output file
#`eog pictures/$output`;

# Subroutines

sub list_files {
  return $file_query if defined $file_query;
  my $half = @files / 2 + .5*(@files % 2);
  for (1..$half) {
    if ($_ < (@files+1)/2) {
      print "$_) ${files[$_-1]} \t" . ($_ + $half) . ") ${files[$_-1+$half]}\n" 
    } else {
      print "$_) ${files[$_-1]}\n"
    }
  }
  return <>;
}

sub get_minimum {
  my ($fn,$line) = @_;
  open my $FH, "< $dir" ."$fn" or die $!;
  my @ele = split ' ', <$FH>;
  my $min = sqrt($ele[0]**2 + $ele[1]**2);
  while (<$FH>) {
    next if $. < $line && $line;
    my @ele = split ' ', $_;
    if ($min < sqrt($ele[0]**2 + $ele[1]**2)) {
      print "Minimum Found for $fn:\nr (m),\t\tz (mm),\t\tvr (m/s),\tvt (m/s)\n";
      print ((join "\t", map sprintf("%.5e", $_), ($min,$ele[2]*1000,sqrt($ele[3]**2 + $ele[4]**2),$ele[5])) . "\n");
      last
    } else {
      $min = sqrt($ele[0]**2 + $ele[1]**2);
    }  
  }
}

sub linear_fit {
  my ($fn,$lines) = @_;
  my @lines = split "\n", `tail -$lines $dir/$fn`;
  my (@x,@y);my $i = 0;
  for (@lines) {
    my @ele = split ' ', $_;
    push @y, sqrt($ele[0]**2 + $ele[1]**2);
    push @x, $ele[2];
  }
  my $y = pdl @y;
  my $x = pdl @x;
  my ($b,$m) = list((fitpoly1d($x,$y, 2))[1]);
  printf "file: $fn\ny-int: %.4fmm\nslope: %.4fmm/m\n", 1000*$b, 1000*$m;
  printf "x-int: %.8fmm\n", -$b/$m*1000 if $m;
  # AG model comparison
  #printf "y(21.5in): %.8fm\n", $m*.5461 + $b;
  return sprintf("$fn crosses at %.4f mm", -$b/$m*1000) if $m;
}

sub guassian_fit {
  my $fn = shift;
  my ($zi,$zf) = @xrange;
  my $data;
  my $wr;
  open my $DATAFILE, "< $dir/$fn" or die "No file exists";
  while (<$DATAFILE>) {
    my ($x,$y,$z) = split " ";
    $_ *= 1000 for ($x,$y,$z);
    next if $z < $zi;
    my $r = sqrt($x**2 + $y**2);
    $wr = $z > ($zf+$zi)/2 && !$wr ? $r : $wr;
    push @{$data}, [$z,$r];
    last if $z > $zf;
  }
  my $zdata = pdl($data)->slice('0')->flat();
  my $rdata = pdl($data)->slice('1')->flat();;
# set initial prameters in a pdl (order in accord with fit function below)
  my $initp = pdl [10*$wr,($zf+$zi)/2,($zf-$zi)/2];

# Weight all y data equally (else specify different weights in a pdl)
  my $wt = 1;
# Use lmfit. Fourth input argument is reference to user-defined 
# subroutine ( here \&linefit ) detailed below.
  my ($yf,$pf,$cf,$if) = lmfit $zdata, $rdata, $wt, \&beamfit, $initp;
  return list $pf;
}

sub beamfit {
  my ($x,$par,$ym,$dyda) = @_;
  my ($w0,$z0,$zR) = map { $par->slice("($_)") } (0..2);

  $ym .= $w0 * sqrt( 1 + (($x-$z0)/$zR)**2 );
  my (@dy) = map {$dyda -> slice(",($_)") } (0..2);

  $dy[0] .= sqrt( 1 + (($x-$z0)/$zR)**2 );
  $dy[1] .= -$w0 * ($x-$z0)    / ($zR**2) / sqrt( 1 + (($x-$z0)/$zR)**2 );
  $dy[2] .= -$w0 * ($x-$z0)**2 / ($zR**3) / sqrt( 1 + (($x-$z0)/$zR)**2 );
}
