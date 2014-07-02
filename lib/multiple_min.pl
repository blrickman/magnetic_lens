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
  open my $FH, "<" . "$dir/$subdir/e11.dat" or die $!;
  my @ele = split ' ', <$FH>;
  my $min = sqrt($ele[0]**2 + $ele[1]**2);

  if ($subdir =~ /I_(\d{2,3})A/) {
    $subdir = $1;
  }

  while (<$FH>) {
    my @ele = split ' ', $_;
    if ($min < sqrt($ele[0]**2 + $ele[1]**2)) {
      print $OUTPUT join(", ", ($subdir, 1000*$ele[2])) . "\n";
      last;
    } else {
      $min = sqrt($ele[0]**2 + $ele[1]**2);
    }
  }
}
