#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min max sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use FASTNFile;
use GeneticsModule;
use UsefulModule;
use CortexBreakpoints;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./check-breakpoints.pl <truth.txt> <breakpoints.txt>\n" .
"  Count how many breakpoints were correctly called. \n" .
"  truth.txt: (1-based)\n" .
"    <chr>:<pos>:<strand>\t<chr>:<pos>:<strand>\t<seq>\n" .
"  Each break should only be represented once.\n";

  exit(-1);
}

if(@ARGV != 2) { print_usage(); }

my ($truth_path, $brkpnt_path) = @ARGV;

# Load breakpoints
my @brkpnts = ();
open(TXT, $truth_path) or die("Cannot open $truth_path");
my $line;
while(defined($line = <TXT>)) {
  chomp($line);
  if($line !~ /^#/ && length($line) > 0) {
    if($line =~ /^([^:]*):(\d+):([\+\-])\t([^:]*):(\d+):([\+\-])\t?([ACGT]*)$/) {
      my $id = @brkpnts;
      push(@brkpnts, {'lchr' => $1, 'lpos' => $2, 'lstrand' => $3,
                      'rchr' => $4, 'rpos' => $5, 'rstrand' => $6,
                      'seq' => $7, 'id' => $id,
                      'counts' => [0, 0]});
    } else {
      die("Bad line: '$line'");
    }
  }
}
close(TXT);

my $num_fp = 0;
my $total_calls = 0;
my $num_with_seq = 0;

my $brkpnt_fh;
open($brkpnt_fh, $brkpnt_path) or die("Cannot open $brkpnt_path");
my $brkpnt = new CortexBreakpoints($brkpnt_fh, $brkpnt_path);

while(my ($seq5p, $seq3p, $gapseq, $hits5p, $hits3p, $cols, $callid) = $brkpnt->next())
{
  if(!defined($seq5p)) { last; }
  $total_calls++;

  if(length($gapseq) > 0) { $num_with_seq++; }

  my ($hitid, $dir, $fl5plen, $fl3plen) = find_breakpoint($seq5p, $hits5p,
                                                          $seq3p, $hits3p,
                                                          $gapseq);

  if($hitid >= 0) { $brkpnts[$hitid]->{'counts'}->[$dir]++; }
  else { $num_fp++; }
}

close($brkpnt_fh);

my $num_found = 0;
my $total_hits = 0;
my $found_fw = 0;
my $found_rv = 0;

for my $brkpnt (@brkpnts) {
  $num_found  += ($brkpnt->{'counts'}->[0] || $brkpnt->{'counts'}->[1] ? 1 : 0);
  $total_hits += $brkpnt->{'counts'}->[0] + $brkpnt->{'counts'}->[1];
  $found_fw += $brkpnt->{'counts'}->[0] > 0 ? 1 : 0;
  $found_rv += $brkpnt->{'counts'}->[1] > 0 ? 1 : 0;
  if(!$brkpnt->{'counts'}->[0] && !$brkpnt->{'counts'}->[1]) {
    print "$brkpnt->{'lchr'}:$brkpnt->{'lpos'}:$brkpnt->{'lstrand'}\t".
          "$brkpnt->{'rchr'}:$brkpnt->{'rpos'}:$brkpnt->{'rstrand'}\t".
          $brkpnt->{'seq'}." ($brkpnt->{'counts'}->[0] $brkpnt->{'counts'}->[1])\n";
  }
}

print "found_fw: $found_fw found_rv: $found_rv\n";

my $nbrkpnts = @brkpnts;
print "Expected $nbrkpnts breakpoints, had $total_calls calls\n";
print "Found ".pretty_fraction($num_found, $nbrkpnts)." breakpoints with " .
      "$total_hits calls ".
      "(" . sprintf('%.2f', $num_found ? $total_hits/$num_found : 0) . " per break)\n";
print "$num_fp false positives\n";
print "$num_with_seq calls had sequence between flanks\n";

# sub mkbrkpnt
# {
#   my ($start0,$end0,$rev0,$start1,$end1,$rev1) = @_;
#   if($start0 > $end0 || $start1 > $end1) { die("Bad blocks $start0 $end0 $start1 $end1"); }
#   return {'lpos' => $rev0 ? $start0 : $end0,   'lstrand' => $rev0 ? '-' : '+',
#           'rpos' => $rev1 ? $end1   : $start1, 'rstrand' => $rev1 ? '-' : '+'};
# }

# Find the largest match to the ref
sub get_largest_match
{
  my $is5p = shift;
  my $flank5plen = shift;
  my ($maxlen,$maxidx) = (0,0);
  for(my $i = 0; $i < @_; $i++) {
    my $h = $_[$i];
    my $len = abs($h->{'end'} - $h->{'start'}) + 1;
    # offset is 1-based
    if($len > $maxlen) {
      if(( $is5p && $h->{'offset'}-1+$len == $flank5plen) ||
         (!$is5p && $h->{'offset'} == 1))
      {
        $maxlen = $len;
        $maxidx = $i;
      }
    }
  }
  my $h = $_[$maxidx];
  my $len = abs($h->{'end'} - $h->{'start'}) + 1;
  if(( $is5p && $h->{'offset'}-1+$len != $flank5plen) ||
     (!$is5p && $h->{'offset'} != 1)) { die("Bad: $is5p $h->{'offset'} $len $flank5plen"); }
  return $_[$maxidx];
}

# returns (hitid, 'fw'/'rv')
sub find_breakpoint
{
  my ($seq5p, $hits5p, $seq3p, $hits3p, $gapseq) = @_;

  my $len5p = length($seq5p);
  my $fl5p = get_largest_match(1, $len5p, @$hits5p);
  my $fl3p = get_largest_match(0, $len5p, @$hits3p);

  # Match lengths
  my $fl5plen = abs($fl5p->{'end'}-$fl5p->{'start'})+1;
  my $fl3plen = abs($fl3p->{'end'}-$fl3p->{'start'})+1;

  # Brute force search to find matching real break
  for(my $i = 0; $i < @brkpnts; $i++) {
    my $b = $brkpnts[$i];
    if($b->{'lchr'} eq $fl5p->{'chrom'} &&
       $b->{'lpos'} == $fl5p->{'end'} && $b->{'lstrand'} eq $fl5p->{'strand'} &&
       $b->{'rchr'} eq $fl3p->{'chrom'} &&
       $b->{'rpos'} == $fl3p->{'start'} && $b->{'rstrand'} eq $fl3p->{'strand'} &&
       $gapseq eq $b->{'seq'}) {
      return ($i, 0, $fl5plen, $fl3plen);
    }
    if($b->{'lchr'} eq $fl3p->{'chrom'} &&
       $b->{'lpos'} == $fl3p->{'start'} && $b->{'lstrand'} ne $fl3p->{'strand'} &&
       $b->{'rchr'} eq $fl5p->{'chrom'} &&
       $b->{'rpos'} == $fl5p->{'end'} && $b->{'rstrand'} ne $fl5p->{'strand'} &&
       $gapseq eq rev_comp($b->{'seq'})) {
      return ($i, 1, $fl5plen, $fl3plen);
    }
  }

  return (-1, 'fw', $fl5plen, $fl3plen);
}
