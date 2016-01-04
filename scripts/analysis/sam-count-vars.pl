#!/usr/bin/env perl

# Takes input from sam-count-bases e.g.:
# NC_009648.1 21 A 66 0 0 66 0

use strict;
use warnings;
use List::Util qw( reduce );

my %dna = ('A'=>0, 'C'=>1, 'G'=>2, 'T'=>3);
my $ncalls = 0;
my $ngood_calls = 0;
my $min_cov_frac = 0.9;

for($ncalls = 0; defined(my $line = <>); $ncalls++) {
  chomp($line);
  my @cols = split('\s', $line);
  my $ref_base = $cols[2];
  my $tot_cov = $cols[3];
  my %cov = ('A'=>$cols[4+0], 'C'=>$cols[4+1], 'G'=>$cols[4+2], 'T'=>$cols[4+3]);
  my $max_base = reduce { $cov{$a} > $cov{$b} ? $a : $b } keys %cov;
  if($max_base ne $ref_base && $cov{$max_base} >= $min_cov_frac*$tot_cov) {
    $ngood_calls++;
    print "".join("\t", "GOOD", @cols, $max_base)."\n";
  } else {
    print "".join("\t", "BAD", @cols, $max_base)."\n";
  }
}

print "$ngood_calls / $ncalls (" .
      sprintf("%.2f", $ncalls ? (100*$ngood_calls)/$ncalls : 0) .
      "%)\n";
