#!/usr/bin/env perl

use strict;
use warnings;

my $i = 0;
my ($kmersize,$total_kmers) = (undef,undef);

print "".join("\t", ('K', 'n_graph_kmers', 'n_link_kmers', 'n_links', 'link_junction_mem'))."\n";

while(my $line = <>) {
  chomp($line);
  if($line =~ /\[graph\] kmer-size: ([0-9]+)/i) {
    $kmersize = $1;
  }
  elsif($line =~ /\[GReader\] ([0-9,]+) kmers, .* filesize/i) {
    $total_kmers = $1;
    $total_kmers =~ s/,//g;
  }
  elsif($line =~ /kmers-with-paths: ([0-9,]+), num paths: ([0-9,]+), path-bytes: (.*)/gi) {
    my ($nkmers,$nlinks,$linkmem) = ($1,$2,$3);
    $nkmers =~ s/,//g;
    $nlinks =~ s/,//g;
    $linkmem =~ s/,//g;
    print "$kmersize\t$total_kmers\t$nkmers\t$nlinks\t$linkmem\n";
    ($kmersize,$total_kmers) = (undef,undef)
  }
  $i++
}

print STDERR "[$0] read $i lines\n";
