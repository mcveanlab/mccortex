#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "Usage: ./clean_bubbles.pl <kmer_size> <in.bubbles>\n";
  print STDERR "  Remove bubbles lacking valid 3p flank. Prints to STDOUT\n";
  exit(-1);
}

if(@ARGV != 2) { print_usage(); }
my $kmer_size = shift;
my $in_file = shift;
my $fh;

if($kmer_size !~ /^\d+$/ || $kmer_size % 2 == 0 || $kmer_size == 0) {
  print_usage("Invalid kmer_size: $kmer_size");
}

if(defined($in_file) && $in_file ne "-") {
  open($fh, $in_file) or print_usage("Cannot open file '$in_file'\n");
} elsif(-p STDIN) {
  open($fh, "<&=STDIN") or print_usage("Cannot read pipe");
} else {
  print_usage("Must specify or pipe in a file");
}

my $line;
my @buf;

while(1)
{
  while(defined($line = <$fh>) && $line =~ /^$/) {}
  if(!defined($line)) { last; }
  chomp($line);
  @buf = ($line);
  for(my $i = 1; $i < 9 && defined($line = <$fh>); $i++) {
    chomp($line);
    push(@buf, $line);
  }
  if(@buf != 9) { print STDERR join("\n", @buf)."\n"; die("Premature file end"); }

  while(defined($line = <$fh>)) {
    chomp($line);
    if($buf[$#buf] eq "" && $line eq "") { last; }
    push(@buf, $line);
  }

  my ($br0, $br1, $fl3) = ($buf[3], $buf[5], $buf[7]);
  $br0 .= $fl3;
  $br1 .= $fl3;
  my $minlen = min(length($br0), length($br1));
  if($minlen > $kmer_size && substr($br0, -$kmer_size) eq substr($br1, -$kmer_size)) {
    print join("\n", @buf)."\n\n";
  }
}

