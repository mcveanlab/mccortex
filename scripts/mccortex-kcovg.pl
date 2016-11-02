#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin . "/perl/";
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use McCortexScripts;
use UsefulModule;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "" .
"Usage: $0 <kmer> <genome_size> <in.ctx>
  Graph sequence coverage and read length from cortex header. Gives same result
  with raw and cleaned graph.

  Example: $0 31 3G sam.clean.ctx > sam.clean.kmercov
\n";

  exit(-1);
}

if(@ARGV != 3) { print_usage(); }

my $maxk = shift(@ARGV);
my $genome_size = shift(@ARGV);
my $graph_file = shift(@ARGV);

if($maxk !~ /^\d+$/ || !($maxk & 1)) { die("Invalid maxk value: $maxk"); }

$maxk = mccortex_maxk($maxk);
$genome_size =~ s/,//g;
$genome_size = str2num($genome_size);

my $cmd = dirname(__FILE__)."/../bin/mccortex$maxk";
my $cmdline = "$cmd view -q -i $graph_file";

if(!(-e $cmd)) {
  die("executable bin/mccortex$maxk doesn't exist -- did you compile for MAXK=$maxk?\n");
}
elsif(!(-x $cmd)) {
  die("bin/mccortex$maxk doesn't appear to be executable\n");
}

my ($ksize, $readlen, $total_seq);
# grep these lines:
#
# kmer size: (\d+)
#
# mean input contig length:\s*([0-9,]+)
# total sequence loaded:\s*([0-9,]+)
#

my $in;
open($in, '-|', $cmdline) or die $!;

while(defined(my $line = <$in>)) {
  chomp($line);
  if($line =~ /kmer size:\s*(\d+)/i) {
    if(defined($ksize)) { die("Duplicate kmer size line: $line"); }
    $ksize = $1;
  }
  if($line =~ /mean input contig length:\s*([0-9,]+)/i) {
    if(defined($readlen)) { die("Duplicate read length line: $line"); }
    $readlen = $1;
    $readlen =~ s/,//g;
  }
  if($line =~ /total sequence loaded:\s*([0-9,]+)/i) {
    if(defined($total_seq)) { die("Duplicate total seq. loaded line: $line"); }
    $total_seq = $1;
    $total_seq =~ s/,//g;
  }
}

# Number of reads * kmers per read
my $nreads = ($total_seq / $readlen);
my $kmers_per_read = ($readlen-$ksize+1);
my $nkmers_read = $nreads * $kmers_per_read;
my $kmercov = sprintf("%.2f", $nkmers_read / $genome_size);

print STDERR "[$0] ksize: $ksize\n";
print STDERR "[$0] total_seq: $total_seq\n";
print STDERR "[$0] readlen: $readlen\n";
print STDERR "[$0] nreads: $nreads\n";
print STDERR "[$0] kmers_per_read: $kmers_per_read\n";
print STDERR "[$0] nkmers_read: $nkmers_read\n";
print STDERR "[$0] kmercov: $kmercov\n";

print int($kmercov+0.5)."\n";
exit(0);
