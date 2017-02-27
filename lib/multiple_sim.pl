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
  'message=s'		=> \my $message,
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
  -m, --message		- message to write to log file describing sim
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

if ($message) {
  local $CWD = $dir;
  my $sim_dir = pop @CWD;
  open my $LOG, ">> simulation.log";
  print $LOG "$sim_dir:\n\t$message\n\n";
}

my $batch = `grep "my \$batch" $sim_file`;
{
  $batch = join "", $batch =~ m/"(.*)"/;
  local $/ = undef;
  open my $BATCH, "< $batch" or die "Can't open file $batch $!";
  $batch = @{ eval <$BATCH> } or die "Bad data in batch file $batch $!";
  $batch -= 1;
}

for our $bt_num (0..$batch) {
  my $sim = do $sim_file or die "Error reading file $sim_file: $@";
  $sim->{dir} = File::Spec->catdir( $dir, $sim->sub_dir() );
  make_path( $sim->dir() );
  print "Completed " . $sim->bt_print();
  $sim->run();
}
