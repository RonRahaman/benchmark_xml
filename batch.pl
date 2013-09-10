#!/usr/bin/perl

use strict;
use POSIX;
$\ = "\n";

my @command = (
  'fox/main',
  'v0/main',
  'v1/main',
  'v3/main',
  'v4/main',
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


