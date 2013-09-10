#!/usr/bin/perl

use strict;
use POSIX;
$\ = "\n";

my @command = (
  'v0_v1_fox/main ./ fox_dom', 
  'v0_v1_fox/main ./ xml_fortran_old', 
  'v0_v1_fox/main ./ xml_fortran_new',
  'v3/main',
  'v0/main',
  'v1/main',
  'fox/main',
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


