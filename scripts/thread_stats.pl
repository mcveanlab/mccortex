#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use CortexScripts;

sub print_usage
{
  if(@_ > 0) { print STDERR map {"Error: $_\n"} @_; }

  print STDERR "" .
"Usage: ./thread_stats.pl [in.txt]
  Convert ctx-thread output to csv, print to STDOUT\n";

  exit(-1);
}

if(@ARGV > 1) { print_usage(); }
my $file = shift;

if(!defined($file)) { $file = "-"; }

my $fh = open_file($file);
my $line;
my $contig_data = "";
my @contigs = ();
my $store = 0;

while($line = <$fh>)
{
  if($line =~ /\[cmd\].*ctx\d+ thread/) {
    if(length($contig_data) > 0) { push(@contigs, $contig_data); $contig_data = ""; }
    $contig_data .= $line;
    $store = 1;
  }
  elsif($store) {
    $contig_data .= $line;
    if($line =~ /\[time\]/ || @contigs > 1000) {
      push(@contigs, $contig_data);
      $contig_data = "";
      $store = 0;
    }
  }
}

if(length($contig_data) > 0) { push(@contigs, $contig_data); $contig_data = ""; }

close($fh);

print STDERR "Found ".scalar(@contigs)." contig runs\n";

my @cols = qw(graphmem pathsmem totalmem kmersize numkmers
              path_type path_num path_bytes path_kmers coloured_paths
              se_reads pe_pairs
              traversal_attempts traversal_successes
              traversal_fail_paths traversal_too_short
              time);

print join(",", @cols)."\n";

for my $data (@contigs)
{
  my %stats;
  ($stats{'graphmem'}) = ($data =~ /graph: (.*?B)/i);
  ($stats{'pathsmem'}) = ($data =~ /paths: (.*?B)/i);
  ($stats{'totalmem'}) = ($data =~ /total: (.*?B)/i);
  ($stats{'kmersize'}) = ($data =~ /kmer-size: (\d+)/i);
  ($stats{'numkmers'}) = ($data =~ /Loaded ([0-9,]+).*? of kmers parsed/i);

  ($stats{'path_type'}) = map {uc($_)} ($data =~ /Paths written to: .*?((?:[sp]e)+).*?\.ctp/i);
  ($stats{'path_num'}) = ($data =~ /([0-9,]+) paths, .*? path-bytes, [0-9,]+ kmers/i);
  ($stats{'path_bytes'}) = ($data =~ /[0-9,]+ paths, (.*?) path-bytes, [0-9,]+ kmers/i);
  ($stats{'path_kmers'}) = ($data =~ /[0-9,]+ paths, .*? path-bytes, ([0-9,]+) kmers/i);
  ($stats{'coloured_paths'}) = ($data =~ /[0-9,]+ paths, .*? path-bytes, [0-9,]+ kmers, coloured paths: ([0-9,]+)/i);

  ($stats{'se_reads'}) = map {uc($_)} ($data =~ /\[stats\] single reads: ([0-9,]+); read pairs: [0-9,]+/i);
  ($stats{'pe_pairs'}) = map {uc($_)} ($data =~ /\[stats\] single reads: [0-9,]+; read pairs: ([0-9,]+)/i);

  ($stats{'traversal_attempts'}) = map {uc($_)} ($data =~ /\[gaps\] traversals succeeded: [0-9,]+ \/ ([0-9,]+)/i);
  ($stats{'traversal_successes'}) = map {uc($_)} ($data =~ /\[gaps\] traversals succeeded: ([0-9,]+) \/ [0-9,]+/i);
  ($stats{'traversal_fail_paths'}) = map {uc($_)} ($data =~ /\[gaps\] failed path check: ([0-9,]+) \/ [0-9,]+/i);
  ($stats{'traversal_too_short'}) = map {uc($_)} ($data =~ /\[gaps\] too short: ([0-9,]+) \/ [0-9,]+/i);

  ($stats{'time'}) = ($data =~ /\[time\] (.*)/i);

  # Strip commas from numbers
  map {$stats{$_} = defined($stats{$_}) ? $stats{$_} : "NA"} @cols;

  for my $col (@cols) {
    if($stats{$col} =~ /^[0-9,]+$/) { $stats{$col} =~ s/,//g; }
    $stats{$col} =~ s/^([0-9]+)([KMG]B?)$/$1.0$2/g;
  }

  print join(",", map {$stats{$_}} @cols)."\n";
}
