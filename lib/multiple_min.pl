#!/usr/bin/perl

use warnings;
use strict;
use File::chdir;

my $dir = shift or die "Choose directory!";
opendir my $DH, $dir or die "Cannot open directory $dir: $!";
open my $OUTPUT, "> $dir/multiple_min.dat";

for my $subdir (sort readdir $DH) {
  next unless -d "$dir/$subdir";
  next if $subdir eq '.' or $subdir eq '..';
  print "$dir/$subdir\n\n";
  system("perl plot.pl -d TE011 -l $dir/$subdir/phase_0.00pi/");
  print "\n\n";
}
