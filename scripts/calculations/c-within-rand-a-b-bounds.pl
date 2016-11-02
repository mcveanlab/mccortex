#!/usr/bin/env perl

use strict;
use warnings;

#
# Sample three random numbers, return how often
#  c >= min(a,b) and c <= max(a,b)
#

use List::Util qw(min max);

# Human genome
my $genome_size = 3100000000;
my $num_loops = 100000; # 100M

my $num_within = 0;

for(my $i = 0; $i < $num_loops; $i++) {
  my ($a,$b,$c) = map {int(rand($genome_size))} 0..2;
  if($c >= min($a,$b) && $c <= max($a,$b)) {
    $num_within++;
  }
}


my $percent = 100 * $num_within / $num_loops;
print "$num_within / $num_loops (".sprintf("%.2f", $percent)."%)\n";
