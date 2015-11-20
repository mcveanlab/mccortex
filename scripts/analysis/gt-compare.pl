#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min max sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../';

use LineReader;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./$0 <a.vcf> <b.vcf>\n" .
"  Generate table of GT concordance\n";

  exit(-1);
}

if(@ARGV != 2) { print_usage(); }

my ($path1, $path2) = @ARGV;
my ($fh1, $fh2);
my ($line1, $line2);
my ($gt1, $gt2);

open($fh1, $path1) or die("Cannot open: $path1");
open($fh2, $path2) or die("Cannot open: $path2");

my $file1 = new LineReader($fh1);
my $file2 = new LineReader($fh2);

my $hdr1 = read_vcf_hdr($file1);
my $hdr2 = read_vcf_hdr($file2);

# print $hdr1;
# print $hdr2;

my %gts = ();
my %gtypes = (); # genotypes seen
my @gstrs = qw(0/0 0/1 0/. 1/0 1/1 1/. ./0 ./1 ./.);
for my $a (@gstrs) { $gts{$a} = {}; }

while(1)
{
  $line1 = $file1->read_line();
  $line2 = $file2->read_line();

  if(!defined($line1) && !defined($line2)) { last; }
  if(!defined($line1)) { die("File $path1 ran out of lines"); }
  if(!defined($line2)) { die("File $path2 ran out of lines"); }
  my @cols1 = split('\t', $line1);
  my @cols2 = split('\t', $line2);
  if(@cols1 < 10) { die("File $path1 missing sample column"); }
  if(@cols2 < 10) { die("File $path2 missing sample column"); }
  for my $c (0, 1, 3, 4) {
    if($cols1[$c] ne $cols2[$c]) { die("Lines don't match:\n$line1\n$line2\n"); }
  }

  if($cols1[9] =~ /^([01\.](?:[\/|][01\.])*)/) {
    $gt1 = $1;
    $gt1 =~ s/\|/\//g;
  } else { die("Bad line [$path1]: $line1"); }
  if($cols2[9] =~ /^([0-9\.]+(?:[\/|][0-9\.]+)*)/) {
    $gt2 = $1;
    $gt2 =~ s/\|/\//g;
  } else { die("Bad line [$path2]: $line2"); }

  $gtypes{$gt1} = 1;
  $gtypes{$gt2} = 1;
  $gts{$gt1}->{$gt2}++;
}

close($fh1);
close($fh2);

# print table
my @keys = keys %gtypes;
print "\t".join("\t", @keys)."\n";
for my $gt1 (@keys) {
  print "$gt1";
  for my $gt2 (@keys) {
    my $v = defined($gts{$gt1}) && defined($gts{$gt1}->{$gt2}) ? $gts{$gt1}->{$gt2} : 0;
    print "\t$v";
  }
  print "\n";
}

sub read_vcf_hdr
{
  my ($lr) = @_;
  my $line;
  my $hdr = "";
  while(defined($line = $lr->read_line())) {
    if($line =~ /^#/) { $hdr .= $line; }
    elsif($line !~ /^\s*$/) {
      $lr->unread_line($line);
      last;
    }
  }
  return $hdr;
}

