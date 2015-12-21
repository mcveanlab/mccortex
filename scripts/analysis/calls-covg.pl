#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min max sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../../libs/bioinf-perl/lib';

use UsefulModule;
use FASTNFile;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./$0 <coverage.txt>\n" .
"  Print stats on VCF contig coverage output.\n";

  exit(-1);
}

if(@ARGV != 1) { print_usage(); }

my ($path) = @ARGV;

my $fastn = open_fastn_file($path);
my ($title,$seq);
my $ncols = 2;
my $entries_per_call = (1+$ncols*3)*2; # bi allelic sites

my $nread = 0;
my $ref_missing = 0;
my $ref_bubbles = 0;

my $alt_raw_missing_covg = 0;
my $alt_raw_missing_edges = 0;
my $alt_raw_simple = 0; # -----
my $alt_raw_with_links = 0; # ---><---
my $alt_raw_messy = 0; # ---<->--

#my $cln_missing = 0;
#my $cln_discoverable = 0;

my $resolve_links = 0;
my $resolve_plain = 0;

while(1)
{
  my ($names,$seqs) = read_multiple_fastn($fastn, $entries_per_call);
  if(!defined($names)) { last; }
  # Analysis here
  $nread++;

  my $refseq = $seqs->[0];
  my $refgraph_ref_edges = $seqs->[1];
  my $refgraph_ref_deg   = $seqs->[2];
  my @refgraph_ref_covg  = split(/\s+/, trim($seqs->[3]));
  my $samgraph_ref_edges = $seqs->[4];
  my $samgraph_ref_deg   = $seqs->[5];
  my @samgraph_ref_covg  = split(/\s+/, trim($seqs->[6]));

  my $altseq = $seqs->[7];
  my $refgraph_alt_edges = $seqs->[8];
  my $refgraph_alt_deg   = $seqs->[9];
  my @refgraph_alt_covg  = split(/\s+/, trim($seqs->[10]));
  my $samgraph_alt_edges = $seqs->[11];
  my $samgraph_alt_deg   = $seqs->[12];
  my @samgraph_alt_covg  = split(/\s+/, trim($seqs->[13]));

  my $c_ref_missing = (min(@refgraph_ref_covg) == 0);
  my $c_ref_bubbles = (min(@refgraph_ref_covg) > 0 && min(@refgraph_alt_covg) > 0);
  my $c_alt_raw_missing_covg = (min(@samgraph_alt_covg) == 0);
  my $c_alt_raw_missing_edges = ($samgraph_alt_deg =~ /.\../);
  my $c_alt_raw_simple = ($samgraph_alt_deg =~ /^.\-+.$/);

  my ($c_alt_raw_messy, $c_alt_raw_with_links) = (0,0);
  if($samgraph_alt_deg =~ /^.\-*[\[\{X]+.*[\]\}X]+.*.$/) { $c_alt_raw_messy = 1; }
  elsif($samgraph_alt_deg =~ /^.\-*([\]\}X]+\-*)+([\[\{X]+\-*)+.$/) { $c_alt_raw_with_links = 1; }

  $ref_missing += $c_ref_missing;
  $ref_bubbles += $c_ref_bubbles;
  $alt_raw_missing_covg += $c_alt_raw_missing_covg;
  $alt_raw_missing_edges += $c_alt_raw_missing_edges;
  $alt_raw_simple += $c_alt_raw_simple;
  $alt_raw_with_links += $c_alt_raw_with_links;
  $alt_raw_messy += $c_alt_raw_messy;

  if(!$c_ref_missing && !$c_ref_bubbles &&
     !$c_alt_raw_missing_covg && !$c_alt_raw_missing_edges &&
     !$c_alt_raw_messy)
  {
    if($c_alt_raw_with_links) { $resolve_links++; }
    else { $resolve_plain++; }
  }
}

print "ncols: $ncols\n";
print "entries_per_call: $entries_per_call\n";
print "nread: $nread\n";
print "ref_missing: $ref_missing\n";
print "ref_bubbles: $ref_bubbles\n";
print "alt_raw_missing_covg: $alt_raw_missing_covg\n";
print "alt_raw_missing_edges: $alt_raw_missing_edges\n";
print "alt_raw_with_links: $alt_raw_with_links\n";
print "alt_raw_messy: $alt_raw_messy\n";
print "resolve_plain: $resolve_plain\n";
print "resolve_links: $resolve_links\n";

close_fastn_file($fastn);


# Read multiple entries from fasta/fastq file
sub read_multiple_fastn
{
  my ($fh,$nentries) = @_;
  my @titles = ();
  my @seqs = ();
  for(my $i = 0; $i < $nentries; $i++) {
    if((($title,$seq) = $fastn->read_next()) && defined($title)) {
      push(@titles, $title);
      push(@seqs, $seq);
    } elsif($i > 0) {
      die("Bad entry: $i $nentries");
    } else {
      return (undef,undef);
    }
  }
  return (\@titles, \@seqs);
}


# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT
# >._1_58210_0_c0_edges
# 85 22 18 81 28 24 21 11 48 48 42 48 14 41 41 41 88 14 24 11 41 12 48 44 18 48 22 82 41 14 84 48 58
# >._1_58210_0_c0_degree
# {-------------------------------}
# >._1_58210_0_c0_covgs
#  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3
# >._1_58210_0_c1_edges
# 84 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 48
# >._1_58210_0_c1_degree
# -...............................-
# >._1_58210_0_c1_covgs
# 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 108
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGGCTATGAATTCTGAAATGGAACTGTTCCAGGT
# >._1_58210_1_c0_edges
# 85 22 18 81 28 24 21 11 48 48 42 48 14 41 41 41 88 14 24 11 41 12 48 44 18 48 22 82 41 14 84 48 58
# >._1_58210_1_c0_degree
# {-------------------------------}
# >._1_58210_1_c0_covgs
#  3  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  3
# >._1_58210_1_c1_edges
# 84 22 18 81 28 24 21 11 48 48 42 48 14 41 41 41 88 14 24 11 41 12 48 44 18 48 22 82 41 14 84 48 48
# >._1_58210_1_c1_degree
# ---------------------------------
# >._1_58210_1_c1_covgs
# 95 96 92 89 93 94 93 92 96 100 102 102 101 106 108 110 112 112 109 112 110 110 111 114 111 108 110 112 112 111 112 110 108

