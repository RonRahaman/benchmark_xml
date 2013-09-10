#!/usr/bin/perl

use strict;
$\ = "\n";
$, = ', ';

my $base = "geometry_base.xml";
my $target = "geometry.xml";
my $size = $ARGV[0];
my @tagstack;

open INFILE, "<", $base;
open OUTFILE, ">", $target;

print OUTFILE '<?xml version="1.0"?>';

if (not defined $size) {
  die "Need to pass a size to generate_geometry.pl";
}

my $i = 1;
while (<INFILE>) {
  s/^\s+|\s+$//g;
  if (/^<!--.*-->/) {
    next;
  }
  elsif (/<.*\/>/) {
    print OUTFILE;
  }
  elsif (/<(.*)>/) {
    my $word;
    ($word) = split(/>|\s/, $1);
    if ($word !~ /^\//) {
      push @tagstack, $word;
      #print "--> pushed $word onto tagstack";

    }
    if (/<\//) {
      pop @tagstack;
      #print "<-- popped $word off of tagstack";
    }
    print OUTFILE;
  }


  if ($i++ > $size) {
    last;
  }
}

for (my $i = scalar(@tagstack)-1 ; $i >= 0; $i--) {
  print OUTFILE "</$tagstack[$i]>";
}

close INFILE;
close OUTFILE;

