#!/usr/bin/perl

use strict;
use POSIX;
$\ = "\n";

my @modes = ('fox_dom', 'xml_fortran_old', 'xml_fortran_new');

unlink("benchmark_xml.txt");

# Range for the number of lines in each xml file
for (my $i = 4; $i <= 1200000 ; $i = 2*$i) {
  system("./generate_geometry.pl $i");
  # Number of trials for this xml file
  for (my $j = 1; $j <= 3; $j++) {
    foreach (@modes) {
      system("v0_v1_fox/main ./ $_"); 
      unless (WIFEXITED(${^CHILD_ERROR_NATIVE})) {
        die "\nexecution of ./batch.pl has been terminated ".
          "(because ./main exited with status ${^CHILD_ERROR_NATIVE})\n";
      }
    }
  }
}


