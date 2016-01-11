#!/usr/bin/env perl

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
"Usage: $0 <bub.txt.gz>\n" .
"  Print a contig for each bubble branch. Prints to STDOUT.\n" .
"  Contigs are named: >BUBBLENAME.branchB\n" .
"  where B is the branch number\n";

  exit(-1);
}

if(@ARGV != 1) { print_usage(); }
my ($file) = (@ARGV,"-");
my $fh;
open($fh, "gzip -fcd $file |") or die("Cannot read file $file: $!");

my $cb = new CortexBubbles($fh);
my ($seq5p, $seq3p, $branches, $flank5p_nkmers, $flank3p_nkmers, $branchlens, $callid);

while(1)
{
  ($seq5p, $seq3p, $branches,
   $flank5p_nkmers, $flank3p_nkmers, $branchlens, $callid) = $cb->next();
  if(!defined($seq5p)) { last; }

  my ($len5p,$len3p) = (length($seq5p), length($seq3p));

  for(my $i = 0; $i < @$branches; $i++) {
    print ">$callid.branch$i:$len5p:$len3p\n";
    print $seq5p.$branches->[$i].$seq3p."\n";
  }
}

close($fh);
