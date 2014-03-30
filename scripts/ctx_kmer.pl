#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use CortexScripts;

# usage: ./ctx_kmer.pl [-a] <in.ctx>

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./ctx_kmer.pl [-a] <in.ctx>
  Print kmer size of a cortex graph file. If -a given, also print whole file.\n";

  exit(-1);
}

my $print_file = 0;

if(@ARGV == 2 && $ARGV[0] =~ /^-?-a(ll)?$/) {
  shift;
  $print_file = 1;
}

if(@ARGV != 1) { print_usage(); }

#
# Pull out the kmer-size of a graph file and print it
# Then proceed to print the rest of the file
# Used in cortex_to_graphviz.pl
#

# file format:
# uint8_t  |    6   | the string "CORTEX" (Note: not null-terminated)
# uint32_t |    1   | version number
# uint32_t |    1   | kmer size (<kmer_size>)

my $path = shift;
if(!defined($path)) { $path = "-"; }

my $fh = open_file($path);

my ($data0,$data1);
if(!read($fh, $data0, 10) || !read($fh, $data1, 4)) { die("Cannot read: $path"); }

print unpack("I", $data1)."\n";

if($print_file) {
  print $data0;
  print $data1;

  # Remaining data
  while(read($fh, $data0, 1)) { print $data0; }
}

close($fh);


# Other comments
# ./scripts/ctx_kmer.pl ./test/subgraph/subgraph1.ctx | { read k; k=$[(($k+31)/32)*32-1]; cat | ./bin/ctx$k view --kmers -; }

# cat ./test/subgraph/subgraph1.ctx | { k=$(./scripts/ctx_kmer.pl); k=$[(($k+31)/32)*32-1]; cat | ./bin/ctx$k view --kmers; }

## attempts to do the same with straight bash and named pipes
## read binary files
# cat test/subgraph/subgraph1.ctx | od -An -t u1 | head

## Names pipes
# mkfifo my_pipe
# cat my_pipe | ./scripts/ctx_kmer.pl &
# cat test/subgraph/subgraph1.ctx | tee pipe | 

# cat input |
#  { .. read n bytes to tmp file; k= use tmp file; cat tmp_file - | ... using $k }
