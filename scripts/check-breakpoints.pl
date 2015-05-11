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
"    <chr0>:<pos0>:<strand0>\t<chr1>:<pos1>:<strand1>\t<seq>\t<chr1>:<pos1>:<strand1>\t<chr0>:<pos0>:<strand0>\t<seq>\n" .
"  Each break should only be represented once, by a single line. Each line\n" .
"  describes forward and backward breaks.\n";

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
    if($line =~ /^([^:]*):(\d+):([\+\-])\t([^:]*):(\d+):([\+\-])\t([ACGT]*)\t([^:]*):(\d+):([\+\-])\t([^:]*):(\d+):([\+\-])\t?([ACGT]*)$/) {
      my $id = @brkpnts;
      push(@brkpnts, {'lchr0' => $1, 'lpos0' => $2, 'lstrand0' => $3,
                      'rchr0' => $4, 'rpos0' => $5, 'rstrand0' => $6,
                      'seq0' => $7,
                      'lchr1' => $8, 'lpos1' => $9, 'lstrand1' => $10,
                      'rchr1' => $11, 'rpos1' => $12, 'rstrand1' => $13,
                      'seq1' => $14,
                      'id' => $id, 'counts' => [0, 0]});
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
    print "F: $brkpnt->{'lchr0'}:$brkpnt->{'lpos0'}:$brkpnt->{'lstrand0'}\t".
             "$brkpnt->{'rchr0'}:$brkpnt->{'rpos0'}:$brkpnt->{'rstrand0'}\t".
             "$brkpnt->{'seq0'}\n";
    print "R: $brkpnt->{'lchr1'}:$brkpnt->{'lpos1'}:$brkpnt->{'lstrand1'}\t".
             "$brkpnt->{'rchr1'}:$brkpnt->{'rpos1'}:$brkpnt->{'rstrand1'}\t".
             "$brkpnt->{'seq1'}\n";
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

# Find the largest match to the ref
# get_largest_match($is5p, $flank5plen, @alignments)
# returns largest alignment
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
    if($b->{'lchr0'} eq $fl5p->{'chrom'} &&
       $b->{'lpos0'} == $fl5p->{'end'} && $b->{'lstrand0'} eq $fl5p->{'strand'} &&
       $b->{'rchr0'} eq $fl3p->{'chrom'} &&
       $b->{'rpos0'} == $fl3p->{'start'} && $b->{'rstrand0'} eq $fl3p->{'strand'} &&
       $gapseq eq $b->{'seq0'}) {
      return ($i, 0, $fl5plen, $fl3plen);
    }
    if($b->{'lchr1'} eq $fl5p->{'chrom'} &&
       $b->{'lpos1'} == $fl5p->{'end'} && $b->{'lstrand1'} eq $fl5p->{'strand'} &&
       $b->{'rchr1'} eq $fl3p->{'chrom'} &&
       $b->{'rpos1'} == $fl3p->{'start'} && $b->{'rstrand1'} eq $fl3p->{'strand'} &&
       $gapseq eq $b->{'seq1'}) {
      return ($i, 1, $fl5plen, $fl3plen);
    }
  }

  return (-1, 'fw', $fl5plen, $fl3plen);
}
