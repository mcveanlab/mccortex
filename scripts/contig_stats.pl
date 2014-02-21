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
  ($stats{'numkmers'}) = ($data =~ /Loaded (\d+).* of kmers parsed/i);

# Loading file

  ($stats{'path_type'}) = ($data =~ /Loading file .*(se|pe|sepe).*\.ctp /i);
  ($stats{'path_num'}) = ($data =~ /([\d,]+) paths, .* path-bytes, [\d,]+ kmers/i);
  ($stats{'path_bytes'}) = ($data =~ /[\d,]+ paths, (.*) path-bytes, [\d,]+ kmers/i);
  ($stats{'path_kmers'}) = ($data =~ /[\d,]+ paths, .* path-bytes, ([\d,]+) kmers/i);

  ($stats{'meancontig'}) = ($data =~ /Lengths: mean: ([\d\.]+)/i);
  ($stats{'mediancontig'}) = ($data =~ /Lengths: .*?median: ([\d\.]+)/i);
  ($stats{'N50contig'}) = ($data =~ /Lengths: .*?N50: ([\d\.]+)/i);
  ($stats{'mincontig'}) = ($data =~ /Lengths: .*?min: ([\d\.]+)/i);
  ($stats{'maxcontig'}) = ($data =~ /Lengths: .*?max: ([\d\.]+)/i);
  ($stats{'totalcontig'}) = ($data =~ /Lengths: .*?total: ([\d\.]+)/i);

  ($stats{'resolve_straight'}) = ($data =~ /Go straight.*\[.*([\d\.]+%).*\]/i);
  ($stats{'resolve_colour'}) = ($data =~ /Go colour.*\[.*([\d\.]+%).*\]/i);
  ($stats{'resolve_path'}) = ($data =~ /Go path.*\[.*([\d\.]+%).*\]/i);

  ($stats{'halt_covg'}) = ($data =~ /No coverage.*\[.*([\d\.]+%).*\]/i);
  ($stats{'halt_colcovg'}) = ($data =~ /No colour covg.*\[.*([\d\.]+%).*\]/i);
  ($stats{'halt_nopaths'}) = ($data =~ /No paths.*\[.*([\d\.]+%).*\]/i);
  ($stats{'halt_pathssplit'}) = ($data =~ /Paths split.*\[.*([\d\.]+%).*\]/i);
  ($stats{'halt_missingpaths'}) = ($data =~ /Paths split.*\[.*([\d\.]+%).*\]/i);
  ($stats{'paths_resolved_juncs'}) = ($data =~ /Paths resolved.*\[.*([\d\.]+%).*\]/i);

  ($stats{'time'}) = ($data =~ /\[time\] (.*)/i);

  # Strip commas from numbers
  map {$stats{$_} = defined($stats{$_}) ? $stats{$_} : "NA"} @cols;
  for my $col (@cols) { if($stats{$col} =~ /^[\d,]+$/) { $stats{$col} =~ s/,//g; } }
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
