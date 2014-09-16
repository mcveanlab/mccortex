#!/usr/bin/perl

use strict;
use warnings;

use File::Basename; # dirname()

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexScripts;
use CortexPaths;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./cortex_mix_paths.pl [--kmer <k>] <graph.ctx> <in.ctp>\n" .
"  Interweave path file with graph file in human readable text format.\n";
  exit(-1);
}

my $maxk = 31;

while(@ARGV > 1 && $ARGV[0] =~ /^-./) {
  if($ARGV[0] =~ /^-?-k(mer)?$/i) {
    my $arg = shift;
    $maxk = shift;
    if(!defined($maxk) || $maxk !~ /^\d+$/) {
      print_usage("$arg <k> requires an argument");
    }
  }
  else { print_usage("Unknown option '$ARGV[0]'"); }
}

# Round kmer-size up to max kmer size supported by an executable
$maxk = int(($maxk+31)/32)*32-1;

if(@ARGV != 2) { print_usage(); }

my ($ctx_path, $ctp_path) = @ARGV;

my $ctp_fh = open_file($ctp_path);
my $ctp_file = new CortexPaths($ctp_fh, $ctp_path);

# Pipe cortex graph through cortex
# graph file reader command
my $cmd = dirname(__FILE__)."/../bin/ctx$maxk";
my $cmdline = "$cmd view --quiet --kmers $ctx_path";
my $in;
open($in, '-|', $cmdline) or die $!;

my %paths = (); # kmer->path string

# Read paths
while(1)
{
  my ($kmer, @paths) = $ctp_file->next();
  if(!defined($kmer)) { last; }
  if(defined($paths{$kmer})) { die("Duplicate kmer: $kmer"); }
  $paths{$kmer} = ctp_path_to_str(@paths);
}

close($ctp_fh);

# Read graph file
while(defined(my $line = <$in>))
{
  if($line =~ /^([ACGT]+)/) {
    my $kmer = $1;
    print $line;
    if(defined($paths{$kmer})) {
      print $paths{$kmer};
    }
  } else {
    die("Bad line: '$line'");
  }
}

close($in);
