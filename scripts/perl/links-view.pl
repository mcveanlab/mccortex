#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename; # dirname()

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use McCortexScripts;
use McCortexLinks;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: $0 [--kmer <k>] <graph.ctx> <in.ctp>\n" .
"  Interweave link file with graph file in human readable text format.\n";
  exit(-1);
}

my $k = 31;

while(@ARGV > 1 && $ARGV[0] =~ /^-./) {
  if($ARGV[0] =~ /^-?-k(mer)?$/i) {
    my $arg = shift;
    $k = shift;
    if(!defined($k) || $k !~ /^\d+$/) {
      print_usage("$arg <k> requires an argument");
    }
  }
  else { print_usage("Unknown option '$ARGV[0]'"); }
}

my $mccortex = dirname(__FILE__)."/../bin/mccortex";
if(!(-x $mccortex)) { die("Have you compiled McCortex with `make`?"); }

if(@ARGV != 2) { print_usage(); }

my ($ctx_path, $ctp_path) = @ARGV;

my $ctp_fh = open_file($ctp_path);
my $ctp_file = new McCortexLinks($ctp_fh, $ctp_path);

# Pipe graph through McCortex
# graph file reader command
my $cmdline = "$mccortex $k view --quiet --kmers $ctx_path";
my $in;
open($in, '-|', $cmdline) or die $!;

my %kmerlinks = (); # kmer->path string

# Read paths
while(1)
{
  my ($kmer, @links) = $ctp_file->next();
  if(!defined($kmer)) { last; }
  if(defined($kmerlinks{$kmer})) { die("Duplicate kmer: $kmer"); }
  $kmerlinks{$kmer} = ctp_link_to_str(@links);
}

close($ctp_fh);

# Read graph file
while(defined(my $line = <$in>))
{
  if($line =~ /^([ACGT]+)/) {
    my $kmer = $1;
    print $line;
    if(defined($kmerlinks{$kmer})) {
      print $kmerlinks{$kmer};
    }
  } else {
    die("Bad line: '$line'");
  }
}

close($in);
