#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min max sum shuffle);

use constant {EDGES_MISSING => 1,
              EDGES_UNRESOLVABLE => 2,
              EDGES_RESOLVABLE_PLAIN_FW => 4,
              EDGES_RESOLVABLE_PLAIN_RV => 8,
              EDGES_RESOLVABLE_PLAIN_BOTH => 12,
              EDGES_RESOLVABLE_LINKS => 16};

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
"Usage: $0 <coverage.txt>\n" .
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
my $alt_missing_covg = 0;

my $alt_missing_edges = 0; # ---.---.--
my $alt_resolve_plain = 0; # --<--<--- or --->-->-- or -------
my $alt_resolve_links = 0; # --->-<---
my $alt_unresolvable = 0;  # ---<->---

my $ref_missing_edges = 0; # ---.---.--
my $ref_resolve_plain = 0; # ---------
my $ref_resolve_links = 0; # --->-<---
my $ref_unresolvable = 0;  # ---<->---

my $ncovg = 0; # number of calls with sample coverage on alt allele
my $ncovg_and_nonref = 0; # non-ref bubble and coverage on alt

my $resolve_alt_links = 0;
my $resolve_alt_plain = 0;

my $resolve_both_links = 0;
my $resolve_both_plain = 0;

my $resolve_alt_simple = 0; # ------ (i.e. no forks)
my $resolve_both_simple = 0; # both ref/alt have no forks

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
  my $c_alt_missing_covg = (min(@samgraph_alt_covg) == 0);
  my $c_alt_edges = get_edges_state($samgraph_alt_deg);
  my $c_ref_edges = get_edges_state($refgraph_ref_deg);
  my $c_samref_edges = get_edges_state($samgraph_ref_deg);

  if(!($c_ref_edges & EDGES_RESOLVABLE_PLAIN_BOTH)) {
    # Attempt to resolve ref allele in sample graph
    $c_ref_edges = $c_samref_edges;
  }

  $ref_missing += $c_ref_missing;
  $ref_bubbles += $c_ref_bubbles;
  $alt_missing_covg += $c_alt_missing_covg;

  if(!$c_alt_missing_covg) {
    $alt_missing_edges += ($c_alt_edges == EDGES_MISSING);
    $alt_resolve_plain += ($c_alt_edges & EDGES_RESOLVABLE_PLAIN_BOTH ? 1 : 0);
    $alt_resolve_links += ($c_alt_edges == EDGES_RESOLVABLE_LINKS);
    $alt_unresolvable += ($c_alt_edges == EDGES_UNRESOLVABLE);
  }

  $ref_missing_edges += ($c_ref_edges == EDGES_MISSING);
  $ref_resolve_plain += ($c_ref_edges & EDGES_RESOLVABLE_PLAIN_BOTH ? 1 : 0);
  $ref_resolve_links += ($c_ref_edges == EDGES_RESOLVABLE_LINKS);
  $ref_unresolvable += ($c_ref_edges == EDGES_UNRESOLVABLE);

  $ncovg += (!$c_alt_missing_covg && $c_alt_edges != EDGES_MISSING);
  $ncovg_and_nonref += (!$c_alt_missing_covg &&
                        $c_alt_edges != EDGES_MISSING &&
                        !$c_ref_bubbles);

  if(!$c_ref_missing && !$c_ref_bubbles &&
     !$c_alt_missing_covg &&
     edges_resovable($c_alt_edges) && edges_resovable($c_ref_edges))
  {
    # Numbers for cortex comparison
    if($c_alt_edges == EDGES_RESOLVABLE_PLAIN_BOTH) { $resolve_alt_simple++; }
    if($c_alt_edges == EDGES_RESOLVABLE_PLAIN_BOTH &&
       $c_ref_edges == EDGES_RESOLVABLE_PLAIN_BOTH) { $resolve_both_simple++; }

    if($c_alt_edges & EDGES_RESOLVABLE_PLAIN_BOTH) { $resolve_alt_plain++; }
    else { $resolve_alt_links++; }

    # Must be able to resolve both alleles in the same direction
    if(($c_ref_edges == EDGES_RESOLVABLE_LINKS) ||
       ($c_alt_edges == EDGES_RESOLVABLE_LINKS)) {
      $resolve_both_links++;
    }
    elsif(($c_ref_edges & $c_alt_edges & EDGES_RESOLVABLE_PLAIN_BOTH) ||
          ($c_samref_edges & $c_alt_edges & EDGES_RESOLVABLE_PLAIN_BOTH)) {
      $resolve_both_plain++;
    }
  }
}

sub edges_resovable {
  return ($_[0] & (EDGES_RESOLVABLE_PLAIN_BOTH | EDGES_RESOLVABLE_LINKS));
}

sub get_edges_state
{
  my ($edges) = @_;
  if(are_edges_missing($edges)) { return EDGES_MISSING; }
  elsif(my $r = are_edge_degrees_plain_resolvable($edges)) { return $r; }
  elsif(are_edge_degrees_links_resolvable($edges)) { return EDGES_RESOLVABLE_LINKS; }
  else { return EDGES_UNRESOLVABLE; }
}

sub are_edges_missing {
  # dot anywhere but first or last position
  # return ($_[0] =~ /^.+\..+$/);
  return ($_[0] =~ /\./);
}

sub are_edge_degrees_plain_resolvable {
  # starts and ends with anything other than a dot, only contains '-'
  my $ret = 0;
  if($_[0] =~ /^[^\.][\-\]\}]+[^\.]$/) { $ret += EDGES_RESOLVABLE_PLAIN_FW; }
  if($_[0] =~ /^[^\.][\-\[\{]+[^\.]$/) { $ret += EDGES_RESOLVABLE_PLAIN_RV; }
  return $ret;
}

sub are_edge_degrees_links_resolvable {
  # <-----> is resolvable
  # -<-->-- is not (first and last are forks from bubble)
  return ($_[0] =~ /^[^\.]\-*([\]\}X].*[\[\{X]|X)\-*[^\.]$/);
}

print "ncols: $ncols\n";
print "entries_per_call: $entries_per_call\n";
print "nread: $nread\n";
print "ref_missing: $ref_missing\n";
print "ref_bubbles: $ref_bubbles\n";
print "ALT missing_covg: $alt_missing_covg\n";
print "ALT missing_edges: $alt_missing_edges\n";
print "ALT resolve_plain: $alt_resolve_plain\n";
print "ALT resolve_links: $alt_resolve_links\n";
print "ALT unresolvable: $alt_unresolvable\n";
print "REF missing_edges: $ref_missing_edges\n";
print "REF resolve_plain: $ref_resolve_plain\n";
print "REF resolve_links: $ref_resolve_links\n";
print "REF unresolvable: $ref_unresolvable\n";
print "calls with ALT covg: $ncovg\n";
print "calls with ALT covg and non-ref: $ncovg_and_nonref\n";
print "resolve_alt_plain: $resolve_alt_plain\n";
print "resolve_alt_links: $resolve_alt_links\n";
print "resolve_both_plain: $resolve_both_plain\n";
print "resolve_both_links: $resolve_both_links\n";
print "resolve_alt_simple: $resolve_alt_simple\n";
print "resolve_both_simple: $resolve_both_simple\n";

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
      warn("Bad entry: $i $nentries");
      return (undef,undef);
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

