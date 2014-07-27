#!/usr/bin/perl

use warnings;
use strict;
use File::chdir;
use Getopt::Long;
use File::Copy qw(copy);

GetOptions(
  'files' 	=> \my $files,
  'axes=s'  	=> \(my $plot_axes = 'r(z)'),
  'scale=s'    	=> \my $scale,
  'output=s'	=> \my $output,
  'help'	=> \my $help,
  'minimum'	=> \my $find_min,
  'yrange=s'	=> \(my $yrange = ''),
  'emfields'	=> \my $emfields,
  'gpexport:s'	=> \my $export,
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
  -y, --yrange 		- set the origin of the y-range
			  Ex: -y 0
  -e, --emfields	- include field files in electron file querry
  -g, --gpexport	- write out the plot file to the directory where the electron files are
			  Ex: -g "foo.gp"  OR  -g (default filename is plot.gp)
  -o, --output		- set output PNG filename
			  Ex: -o filename.png
  -h, --help		- Shows this message

END

my $check_dir = opendir(my $DIR, $dir) or die "Check directory path: $dir\n\n$HELP_MSG $!";  # Makes sure $dir exists
my $HOME_DIR = $CWD; #Save home dir for later use

# Print Help Message
if ($help || !$check_dir) {
print $HELP_MSG;
exit 1;
}

# Grab all electron files from DIR
my @files = sort grep($_=~ /e-?\d{1,3}.dat/ , readdir($DIR));
push @files, 'e105_fields.dat' if $emfields;
closedir($DIR);
my @plot_files = @files;

# If --emfields is triggered, --files should default to TRUE

$files = $emfields ? 1 : $files;

# If --files is flagged, keep relavent electron files
if ($files) {
  if ($files eq 'TODO') {
    @plot_files = $files[0]; # possibly add short cut to querry to only plot e10?
  } else {
    # Querries for a list of numbers (by assigned number) or shortcut 
    # and plot all using 'a', which is the default behavior.
    print "Enter numbers of files to be plotted separated by spaces or enter 'a' to plot all files:\n";
    list_files();
    my $file_query = <>;
    @plot_files = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
  }
}

# If --minimum is flagged, find the minimum of one electron path

if ($find_min) {
  if ($find_min eq 'TODO') {
    my @min_files = $files[0]; # possibly add short cut to querry to only plot e10?
  } elsif (@plot_files == 1) {
    get_minimum($plot_files[0]); #fits the same path selected by the fit;
  } else {
    # Select an electron path
    print "Enter the number of the file to be fitted:\n";
    list_files();
    my $file_query = <>;
    my @min_files = $file_query =~ 'a' ? @files : map ($files[$_-1] , (split ' ', $file_query));
    get_minimum($_) for @min_files;
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
my $axes = join ":", map ($axis->{$_}{marker} . "/(" . $axis->{$_}{scale} .")", ($x,$y));
my $plot_files = join(",\\\n",map("'$dir${_}' u $axes w l lw 2 title '$_'",@plot_files));

# Set up labels for gplot script
my @label = map {"$_ (" . $axis->{$_}{unit} . ")"} ($x,$y);

# Set up gplot script
my $GNUPLOT_SCRIPT = <<END;
set terminal pngcairo enhanced solid font "times,18" size 750, 500 
set output "pictures/$output"
set title "$title"
set yrange [$yrange:]
set xlabel "${label[0]}"
set ylabel "${label[1]}"
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
`eog pictures/$output`;

# Subroutines

sub list_files {
  my $half = @files / 2 + .5*(@files % 2);
  for (1..$half) {
    if ($_ < (@files+1)/2) {
      print "$_) ${files[$_-1]} \t" . ($_ + $half) . ") ${files[$_-1+$half]}\n" 
    } else {
      print "$_) ${files[$_-1]}\n"
    }
  }
}

sub get_minimum {
  my $fn = $_[0];
  open my $FH, "< $dir" ."@_" or die $!;
  my @ele = split ' ', <$FH>;
  my $min = sqrt($ele[0]**2 + $ele[1]**2);

  while (<$FH>) {
    my @ele = split ' ', $_;
    if ($min < sqrt($ele[0]**2 + $ele[1]**2)) {
      print "Minimum Found for $fn:\nr (m),\t\tz (mm),\t\tvr (m/s),\tvt (m/s)\n";
      print ((join "\t", map sprintf("%.4e", $_), (${min},$ele[2]*1000,sqrt($ele[3]**2 + $ele[4]**2),$ele[5])) . "\n");
      last
    } else {
      $min = sqrt($ele[0]**2 + $ele[1]**2);
    }  
  }
}

