#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use CortexScripts;
use CortexLinks;
use FASTNFile;
use GeneticsModule;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "" .
"Usage: ./cortex_graph_links.pl [--print-header] <in.ctp> <kmers.fa>\n";

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
my $ctp_file = new CortexLinks($ctp_fh,$ctp_path);

my ($title,$seq);
my @tgt_seqs  = (); # array of sequences requested
my @tgt_kmers = (); # array of kmer-keys requested
my %kmerlinks = (); # kmer->link arr

# Load kmers requested
while((($title,$seq) = $kmer_file->next()) && defined($title)) {
  # dna_rev_comp_group() is kmer-key
  my $kmer = dna_rev_comp_group($seq);
  if(defined($kmerlinks{$kmer})) { die("Duplicate kmer in input: $kmer [$kmer_path]"); }
  $kmerlinks{$kmer} = [];
  push(@tgt_seqs, $seq);
  push(@tgt_kmers, $kmer);
}

# Read paths
while(1)
{
  my ($kmer, @links) = $ctp_file->next();
  if(!defined($kmer)) { last; }

  # Store paths if requested
  if(defined($kmerlinks{$kmer})) {
    if(@{$kmerlinks{$kmer}} > 0) { die("Duplicate kmer: $kmer"); }
    $kmerlinks{$kmer} = \@links;
  }
}

if($print_header) {
  print $ctp_file->{'_header'};
}

# Print links found
for(my $i = 0; $i < @tgt_seqs; $i++) {
  my ($kmer,$seq) = ($tgt_kmers[$i],$tgt_seqs[$i]);
  my $nlinks = @{$kmerlinks{$kmer}};
  if($nlinks > 0) {
    if($kmer ne $seq) { print "$kmer $nlinks ($seq)\n"}
    else { print "$kmer $nlinks\n"; }
    for my $p (@{$kmerlinks{$kmer}}) {
      ctp_print_link($p);
    }
  }
}

close($ctp_fh);
close($kmer_fh);
