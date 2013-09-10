#!/usr/bin/perl

use strict;

my $head = 'tallies_head.xml';
my $tail = 'tallies_tail.xml';
my $outfile = 'tallies.xml';
my $size = $ARGV[0];

open OUTFILE, '>', $outfile or die;

open HEAD, '<', $head or die;
while (<HEAD>) {
  print OUTFILE;
}
close HEAD;

for (my $i=100; $i < 100+$size; $i++) {
  print OUTFILE $i.' ';
  if (($i+1)%5 == 0) {
    print OUTFILE "\n";
  }
}

open TAIL, '<', $tail or die;
while(<TAIL>) {
  print OUTFILE;
}
close TAIL;

close HEAD;



