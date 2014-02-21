#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum min max);

sub print_usage
{
  if(@_ > 0) { print STDERR map {"Error: $_\n"} @_; }
  
  print STDERR "" .
"Usage: ./contig_stats.pl [in.txt]
  Convert ctx-contig output to csv, print to STDOUT\n";

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
  if($line =~ /\[cmd\].*ctx\d+ contigs/) {
    if(length($contig_data) > 0) { push(@contigs, $contig_data); $contig_data = ""; }
    $contig_data .= $line;
    $store = 1;
  }
  elsif($store) {
    $contig_data .= $line;
    if($line =~ /\[time\]/) {
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
              path_type path_num path_bytes path_kmers
              meancontig mediancontig N50contig mincontig maxcontig totalcontig
              resolve_straight resolve_colour resolve_path
              halt_covg halt_colcovg halt_nopaths halt_pathssplit halt_missingpaths
              paths_resolved_juncs time);

print join(",", @cols)."\n";

for my $data (@contigs)
{
  my %stats;
  ($stats{'graphmem'}) = ($data =~ /graph: (.*B)/i);
  ($stats{'pathsmem'}) = ($data =~ /paths: (.*B)/i);
  ($stats{'totalmem'}) = ($data =~ /total: (.*B)/i);
  ($stats{'kmersize'}) = ($data =~ /kmer-size: (\d+)/i);
  ($stats{'numkmers'}) = ($data =~ /Loaded (\d+).*? of kmers parsed/i);

# Loading file

  ($stats{'path_type'}) = ($data =~ /Loading file .*?(se|pe|sepe).*?\.ctp /i);
  ($stats{'path_num'}) = ($data =~ /([0-9,]+) paths, .*? path-bytes, [0-9,]+ kmers/i);
  ($stats{'path_bytes'}) = ($data =~ /[0-9,]+ paths, (.*?) path-bytes, [0-9,]+ kmers/i);
  ($stats{'path_kmers'}) = ($data =~ /[0-9,]+ paths, .*? path-bytes, ([0-9,]+) kmers/i);

  ($stats{'meancontig'}) = ($data =~ /Lengths: mean: ([0-9\.]+)/i);
  ($stats{'mediancontig'}) = ($data =~ /Lengths: .*?median: ([0-9\.]+)/i);
  ($stats{'N50contig'}) = ($data =~ /Lengths: .*?N50: ([0-9\.]+)/i);
  ($stats{'mincontig'}) = ($data =~ /Lengths: .*?min: ([0-9\.]+)/i);
  ($stats{'maxcontig'}) = ($data =~ /Lengths: .*?max: ([0-9\.]+)/i);
  ($stats{'totalcontig'}) = ($data =~ /Lengths: .*?total: ([0-9\.]+)/i);

  ($stats{'resolve_straight'}) = ($data =~ /Go straight.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'resolve_colour'}) = ($data =~ /Go colour.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'resolve_path'}) = ($data =~ /Go path.*?\[.*?([0-9\.]+%).*?\]/i);

  ($stats{'halt_covg'}) = ($data =~ /No coverage.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'halt_colcovg'}) = ($data =~ /No colour covg.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'halt_nopaths'}) = ($data =~ /No paths.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'halt_pathssplit'}) = ($data =~ /Paths split.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'halt_missingpaths'}) = ($data =~ /Paths split.*?\[.*?([0-9\.]+%).*?\]/i);
  ($stats{'paths_resolved_juncs'}) = ($data =~ /Paths resolved.*?\[.*?([0-9\.]+%).*?\]/i);

  ($stats{'time'}) = ($data =~ /\[time\] (.*)/i);

  # Strip commas from numbers
  map {$stats{$_} = defined($stats{$_}) ? $stats{$_} : "NA"} @cols;
  for my $col (@cols) { if($stats{$col} =~ /^[0-9,]+$/) { $stats{$col} =~ s/,//g; } }
  print join(",", map {$stats{$_}} @cols)."\n";
}


sub open_file
{
  my ($file) = @_;
  if(!defined($file)) { die("No file specified to open"); }
  my $handle;
  # -p checks if connected to a pipe
  if($file ne "-") {
    open($handle, $file) or die("Cannot open file: $file\n");
  }
  elsif(-p STDIN) {
    open($handle, "<&=STDIN") or die("Cannot read pipe");
  }
  else { die("Must specify or pipe in a file"); }
  return $handle;
}
