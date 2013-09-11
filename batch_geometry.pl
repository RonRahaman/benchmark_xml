#!/usr/bin/perl

use strict;
use POSIX;
$\ = "\n";

my @command = (
  'fox/main geometry',
  'v0/main geometry',
  'v1/main geometry',
  'v3/main geometry',
  'v4/main geometry',
);

unlink("benchmark_xml.txt");

# Range for the number of lines in each xml file
for (my $i = 2; $i <= 1200000 ; $i = 2*$i) {
  system("./generate_geometry.pl $i");
  # Number of trials for this xml file
  for (my $j = 1; $j <= 1; $j++) {
    foreach (@command) {
      system($_);
      unless (WIFEXITED(${^CHILD_ERROR_NATIVE})) {
        die "\nexecution of ./batch.pl has been terminated ".
          "(because ./main exited with status ${^CHILD_ERROR_NATIVE})\n";
      }
    }
  }
}


