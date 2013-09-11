#!/usr/bin/perl

use strict;
use POSIX;
$\ = "\n";

my @command = (
  'v0/main tallies',
  'v1/main tallies',
  'v3/main tallies',
  'v4/main tallies',
);

unlink("tallies_benchmark.txt");

# Range for the number of lines in each xml file
for (my $i = 2; $i <= 20000; $i = 2*$i) {
  system("./generate_tallies.pl $i");
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
