#!/usr/bin/env perl

#
# Stand alone script print maximal substrings in O(N^3) time
# An O(N^2) time dynamic programming solution written in C can be found in:
#   libs/maximal_substrs
#

use strict;
use warnings;

if(@ARGV != 1) {
  print STDERR "usage ./find_maximal_substrings.pl <file>\n";
  print STDERR "  Find all maximal substrings in O(n^3) time and O(n^2) space\n";
  exit(-1);
}

my $file = shift(@ARGV);

# noccur(str,s)
#   Return number of (potentially overlapping) occurrance of s in str
#   see: http://stackoverflow.com/a/2114234/431087
sub noccur
{
  my ($str, $s) = @_;
  my $adj = length($s) - 1;
  die "Search string cannot be empty!" if $adj < 0;

  my $count = 0;
  while ( $str =~ /\Q$s/g ) {
    pos $str -= $adj;
    $count++;
  }
  return $count;
}

# Return substring starting at i up to and including j
sub getsubstr
{
  my ($str,$i,$j) = @_;
  return substr($str,$i,$j-$i+1);
}

# find all occurances of a given substring
#  str  string to search
#  l    length of string str
#  i    start location (inclusive)
#  i    end location (inclusive)
sub findall
{
  my ($str,$l,$i,$j) = @_;
  my $n = noccur($str,getsubstr($str,$i,$j));
  if($n == 1) { return undef; }
  if($n <  1) { die("Wat: $n"); }
  while($i > 0    && noccur($str,getsubstr($str,$i-1,$j  )) == $n) { $i--; }
  while($j+1 < $l && noccur($str,getsubstr($str,$i  ,$j+1)) == $n) { $j++; }
  return (getsubstr($str, $i, $j), $n);
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

# Read all lines into a string
my $str = "";
my $fh = open_file($file);
while(defined(my $line = <$fh>)) { $str .= $line; }
close($fh);

my ($i,$j);
my $l = length($str);
my %h = ();
for($i=0; $i < $l; $i++) {
  for($j=$i; $j < $l; $j++) {
    my ($s,$n) = findall($str,$l,$i,$j);
    if(defined($s) && !defined($h{$s})) { $h{$s} = $n; }
  }
}

# map {$h{$_} = noccur($str,$_)} keys(%h);

# Print unique repeats
# sort by counts, then by substring alphabetically
for my $s (sort { ($h{$a} <=> $h{$b}) || ($a cmp $b) } keys(%h)) {
  print "$s\t$h{$s}\n";
}
