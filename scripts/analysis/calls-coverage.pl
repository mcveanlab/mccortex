#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min max sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../../libs/bioinf-perl/lib';

use UsefulModule;

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

my $fh;
if($path eq "-") { print STDERR "Reading STDIN\n"; $fh = open_stdin(); }
else { open($fh, $path) or die("Cannot open: $path"); }

my @nalleles = ();

my $ref_missing = 0;
my $nref_bubbles = 0;
my $raw_missing = 0;
my $cln_missing = 0;

my $raw_discoverable = 0;
my $cln_discoverable = 0;

#
# ref,raw,clean
#
my @d;
my $ncalls;
for($ncalls = 0; scalar(@d = read_lines()); $ncalls++)
{
  # count number of alleles (including REF)
  my $n = (scalar(@d)-2) / 27;
  $nalleles[$n]++;

  # my @ml = ($d[10], $d[19], $d[28], $d[39], $d[48], $d[57]);
  # print "--\n";
  # for my $m (@ml) {
  #   print "$m\n";
  # }

  # Coverage
  my @ref_ref = ($d[10] =~ /(\d+)/g);
  my @raw_ref = ($d[19] =~ /(\d+)/g);
  my @cln_ref = ($d[28] =~ /(\d+)/g);
  my @ref_alt = ($d[39] =~ /(\d+)/g);
  my @raw_alt = ($d[48] =~ /(\d+)/g);
  my @cln_alt = ($d[57] =~ /(\d+)/g);

  # my @all = (@ref_ref, @ref_alt, @raw_ref, @raw_alt, @cln_ref, @cln_alt);
  # for my $t (@all) {
  #   if($t !~ /^\d+$/) { die("What: $t"); }
  # }

  # print "==\n";
  # print join(',',@ref_ref)."\n";
  # print join(',',@ref_alt)."\n";

  if(min(@ref_ref) == 0) { $ref_missing++; }
  if(min(@ref_alt) > 0) { $nref_bubbles++; }
  if(min(@raw_alt) == 0) { $raw_missing++; }
  if(min(@cln_alt) == 0) { $cln_missing++; }

  if(min(@ref_alt) == 0 && min(@raw_alt) > 0) { $raw_discoverable++; }
  if(min(@ref_alt) == 0 && min(@cln_alt) > 0) { $cln_discoverable++; }

  # if(min(@ref_ref) == 0) {
  #   my @ml = ($d[10], $d[19], $d[28], $d[39], $d[48], $d[57]);
  #   print "--\n";
  #   print "$d[0]\n";
  #   print "$d[1]\n";
  #   for my $m (@ml) {
  #     print "$m\n";
  #   }
  # }
}

close($fh);

print "ref_missing: ".pretty_fraction($ref_missing, $ncalls),"\n";
print "nref_bubbles: ".pretty_fraction($nref_bubbles, $ncalls),"\n";
print "raw_missing: ".pretty_fraction($raw_missing, $ncalls),"\n";
print "cln_missing: ".pretty_fraction($cln_missing, $ncalls),"\n";
print "discoverable in raw: ".pretty_fraction($raw_discoverable, $ncalls),"\n";
print "discoverable in clean: ".pretty_fraction($cln_discoverable, $ncalls),"\n";

print "Number of alleles:\n";
for(my $i = 1; $i < @nalleles; $i++) {
  print "  $i: ".(defined($nalleles[$i]) ? $nalleles[$i]-1 : 0)."\n";
}

#
# Functions
#

# Buffer for peeking at next line
my $nl = undef;

sub mypeekline
{
  if(!defined($nl)) { $nl = <$fh>; }
  return $nl;
}

sub myreadline
{
  my $line;
  if(defined($nl)) { $line = $nl; $nl = undef; }
  else { $line = <$fh>; }
  if(defined($line)) { chomp($line); }
  return $line;
}

sub read_lines
{
  my @lines = ();
  my $line;
  my $i;
  while(1) {
    for($i = 0; $i < 29; $i++) {
      $line = myreadline();
      if(!defined($line)) {
        if($i == 0) { return @lines; }
        else { die("Missing lines"); }
      }
      push(@lines, $line);
    }
    my $next = mypeekline();
    if(varcode($lines[0]) ne varcode($next)) { last; }
  }
  return @lines;
}

sub varcode
{
  my ($line) = @_;
  if(!defined($line)) { return ""; }
  my ($code) = ($line =~ /^>(.*)_\d+$/i);
  return $code;
}


# 3*3*3+2=29 lines per ref/alt
# So 2*29 per biallelic

# ref allele example:
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 85 22 18 81 28 24 21 11 48 48 42 48 14 41 41 41 88 14 24 11 41 12 48 44 18 48 22 82 41 14 84 48 58
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# <===============================>
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
#  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 84 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 48
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# =!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
# 103  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 114
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 84 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 48
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# =!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=
# >._1_58210_0
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
# 103  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 114
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 85 22 18 81 28 24 21 11 48 48 42 48 14 41 41 41 88 14 24 11 41 12 48 44 18 48 22 82 41 14 84 48 58
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# <===============================>
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
#  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 84 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 48
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# =!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
# 103  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 114
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_edges
# 84 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 48
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_degree
# =!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=
# >._1_58210_1
# CATCCCAGGGGAGGGTACAGAGGAGCTGATGACTATGAATTCTGAAATGGAACTGTTCCAGGT_covgs
# 103  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 114