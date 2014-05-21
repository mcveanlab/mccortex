#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexBubbles;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./cortex_bubbles.pl in.txt\n";

  exit(-1);
}

if(@ARGV > 1) { print_usage(); }
my $file = $ARGV[0];
if(!defined($file)) { $file = "-"; }
my $fh;
open($fh, $file) or die("Cannot read file $file");

my $cb = new CortexBubbles($fh);
my ($seq5p, $seq3p, $branches, $flank5p_nkmers, $flank3p_nkmers, $branchlens, $callid);

while(1)
{
  ($seq5p, $seq3p, $branches,
   $flank5p_nkmers, $flank3p_nkmers, $branchlens, $callid) = $cb->next();
  if(!defined($seq5p)) { last; }

  print "BUBBLE $callid\n";
  print ">flank5p $flank5p_nkmers nkmers=$flank5p_nkmers\n$seq5p\n";
  print ">flank3p $flank5p_nkmers nkmers=$flank3p_nkmers\n$seq3p\n";
  print "". join('', map {">branch$_ nkmers=$branchlens->[$_]\n$branches->[$_]\n"} 0..(@$branches-1));
  print "\n";
}

close($fh);
