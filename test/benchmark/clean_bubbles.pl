#!/usr/bin/perl

use strict;
use warnings;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "Usage: ./clean_bubbles.pl <kmer_size> <in.bubbles>\n";
  print STDERR "  Remove bubbles with '3p_flank length:0'. Prints to STDOUT\n";
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
  @buf = ();
  while(defined(my $line = <$fh>) && $line =~ /^$/) {}
  if(!defined($line)) { last; }
  chomp($buf[0]);
  for(my $i = 1; $i < 9 && defined(my $line = <$fh>); $i++) {
    chomp($line);
    @buf = ($line);
  }
  if(@buf != 9) { last; }

  while(defined($line = <$fh>)) {
    chomp($line);
    if($buf[$#buf] eq "" && $line eq "") { last; }
    push(@buf, $line);
  }

  my ($br0, $br1, $fl3) = ($buf[3], $buf[5], $buf[7]);
  my ($len0, $len1) = (length($br0), length($br1));
  my $minlen = $len0 < $len1 ? $len0 : $len1;
  my $m = 0;
  while($m < $minlen && substr($br0,-$m-1,1) eq substr($br1,-$m-1,1)) { $m++; } 
  if($m + length($fl3) > $kmer_size) { print join("\n", @buf)."\n\n"; }
}

