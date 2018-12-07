#!/usr/bin/env perl

use strict;
use warnings;

# Walking through 100 copies of Alu, assume ~300 kmers each
# ~1 million copies Alu in the human genome

my $alu_length = 300;
my $n_copies = 100;

# Power in perl is **
my $k = 3;
my $m = 2**22 * 8; # 4MB, 33.6 million bits
my $success = 1;
my $num_rep_kmers = $n_copies*$alu_length;
my $false_pos_rate;

for(my $i = 0; $i < $num_rep_kmers; $i++) {
  $false_pos_rate = bloom_false_pos($k,$m,$i);
  $success *= (1-$false_pos_rate);
}

# Complete => traversing $num_rep_kmers without a false positive
print "$k hash functions; $m bits;\n";
print "$num_rep_kmers false positive rate: $false_pos_rate\n";
print "complete success rate: $success\n";
print "complete failure rate: ".(1-$success)."\n";

# k is the number of hash functions
# m is the number of total bits
# n is the number of bits set
sub bloom_false_pos
{
  my ($k,$m,$n) = @_;
  return (1 - exp(1)**(-$k * $n / $m))**$k;
}

