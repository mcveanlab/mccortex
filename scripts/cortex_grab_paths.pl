#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexScripts;
use CortexPaths;
use FASTNFile;
use GeneticsModule;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "" .
"Usage: ./cortex_graph_paths.pl [--print-header] <in.ctp> <kmers.fa>\n";

  exit(-1);
}

my $print_header = 0;

while(@ARGV > 2) {
  my $arg = shift;
  if($arg eq "--print-header") { $print_header = 1; }
  else { print_usage("Bad option: $arg"); }
}

if(@ARGV != 2) { print_usage(); }

my ($ctp_path,$kmer_path) = @ARGV;

my $kmer_fh = open_file($kmer_path);
my $kmer_file = new FASTNFile($kmer_fh,$kmer_path);

my $ctp_fh = open_file($ctp_path);
my $ctp_file = new CortexPaths($ctp_fh,$ctp_path);

my ($title,$seq);
my @tgt_seqs  = (); # array of sequences requested
my @tgt_kmers = (); # array of kmer-keys requested
my %paths = (); # kmer->path arr

# Load kmers requested
while((($title,$seq) = $kmer_file->next()) && defined($title)) {
  # dna_rev_comp_group() is kmer-key
  my $kmer = dna_rev_comp_group($seq);
  if(defined($paths{$kmer})) { die("Duplicate kmer in input: $kmer [$kmer_path]"); }
  $paths{$kmer} = [];
  push(@tgt_seqs, $seq);
  push(@tgt_kmers, $kmer);
}

# Read paths
while(1)
{
  my ($kmer, @paths) = $ctp_file->next();
  if(!defined($kmer)) { last; }

  # Store paths if requested
  if(defined($paths{$kmer})) {
    if(@{$paths{$kmer}} > 0) { die("Duplicate kmer: $kmer"); }
    $paths{$kmer} = \@paths;
  }
}

if($print_header) {
  print $ctp_file->{'_header'};
}

# Print paths found
for(my $i = 0; $i < @tgt_seqs; $i++) {
  my ($kmer,$seq) = ($tgt_kmers[$i],$tgt_seqs[$i]);
  my $npaths = @{$paths{$kmer}};
  if($npaths > 0) {
    if($kmer ne $seq) { print "$kmer $npaths ($seq)\n"}
    else { print "$kmer $npaths\n"; }
    for my $p (@{$paths{$kmer}}) {
      ctp_print_path($p);
    }
  }
}

close($ctp_fh);
close($kmer_fh);
