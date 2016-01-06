#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexBreakpoints;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: $0 <brk.gz>\n";

  exit(-1);
}

if(@ARGV > 1) { print_usage(); }
my ($file) = (@ARGV, "-");
my $fh;
open($fh, "gzip -fcd $file |") or die("Cannot read file $file: $!");

my $cb = new CortexBreakpoints($fh,$file);
my ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid);

while(1)
{
  ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid) = $cb->next();
  if(!defined($seq5p)) { last; }

  my @strs5p = map {$_->{'chrom'}.":".$_->{'start'}.'-'.$_->{'end'}} @$flank5p_refs;
  my @strs3p = map {$_->{'chrom'}.":".$_->{'start'}.'-'.$_->{'end'}} @$flank3p_refs;

  print "$callid\n";
  print ">flank5p chrs=".join(',', @strs5p)."\n$seq5p\n";
  print ">flank3p chrs=".join(',', @strs3p)."\n$seq3p\n";
  print ">path cols=".join(',', @$cols)."\n$pathseq\n";
  print "\n";
}

close($fh);
