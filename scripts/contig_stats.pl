#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexScripts;
use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "" .
"Usage: ./contig_stats.pl <in.fa|fq>
  Length distribution stats.
";

  exit(-1);
}

if(@ARGV != 1) { print_usage(); }

my $path = shift(@ARGV);
my @lengths = ();

my $fastn = open_fastn_file($path);
my ($title,$seq);
while((($title,$seq) = $fastn->read_next()) && defined($title)) {
  push(@lengths, length($seq));
}
close_fastn_file($fastn);

# min max median mode n50
@lengths = sort{$a <=> $b} @lengths;
my $ncontigs = scalar(@lengths);
my $sum = sum(@lengths);

if($ncontigs == 0) { print STDERR "[contig_stats.pl] No sequences\n"; exit -1; }

my $median = find_median(@lengths);
my $mode = find_mode(@lengths);
my $n50 = find_N50($sum,@lengths);
my $linewidth = 30;

# Some lines $linewidth+2 for 1 decimal place
print_cols("contigs:", num2str($ncontigs), $linewidth);
print_cols(" length:", num2str($sum), $linewidth);
print_cols("    min:", num2str($lengths[0]), $linewidth);
print_cols("    max:", num2str($lengths[$ncontigs-1]), $linewidth);
print_cols("   mean:", num2str($sum/$ncontigs,',',1,1), $linewidth+2);
print_cols(" median:", num2str($median,',',1,1), $linewidth+2);
print_cols("   mode:", num2str($mode), $linewidth);
print_cols("    N50:", num2str($n50), $linewidth);

sub print_cols
{
  my ($ltext,$rtext,$width) = @_;
  my $llen = length($ltext);
  my $rlen = length($rtext);
  my ($padding, $gap);
  if($llen + $rlen < $width) { $padding = $width - $llen - $rlen; }
  if($padding > 2) { $gap = " ".("."x($padding-2))." "; }
  else { $gap = " "x$padding; }
  print "[contig_stats.pl] ".$ltext.$gap.$rtext."\n";
}

# @_ must be sorted
# If there is a tie for the mode, return the first value
# e.g. 1,2,2,3,4,4 => 2
sub find_mode
{
  if(@_ == 0) { return undef; }
  my ($maxidx,$maxrun,$run) = (0,0,0);
  for(my $i = 1; $i < @_; $i++) {
    if($_[$i] == $_[$i-1]) {
      $run++;
      if($run > $maxrun) { $maxidx = $i; $maxrun = $run; }
    }
    else { $run = 1; }
  }
  return $_[$maxidx];
}

# @_ must be sorted
sub find_median
{
  if(@_ == 0) { return undef; }
  if(@_ % 2 == 0) {
    return ($_[@_ / 2 - 1] + $_[@_ / 2]) / 2;
  } else {
    return $_[int(@_ / 2)];
  }
}

# usage: find_N50($total,@array)
# @array must be sorted
sub find_N50
{
  my $total = shift;
  my ($i, $sum, $half) = (0, 0, $total/2);

  if(@_ == 0 || $half == 0) { return undef; }

  for($i = @_; $i > 0 && $sum < $half; $i--) {
    $sum += $_[$i-1];
  }

  return $_[$i];
}
