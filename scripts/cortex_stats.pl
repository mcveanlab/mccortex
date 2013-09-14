#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use IPC::Open2;
use List::Util qw(sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

sub print_usage
{
  if(@_ > 0) { print STDERR map {"Error: $_\n"} @_; }
  
  print STDERR "" .
"Usage: ./cortex_stats.pl <in.ctx>
  Prints cortex graph stats\n";

  exit(-1);
}

if(@ARGV != 1) { print_usage(); }
my $file = shift;

if(!(-r $file)) { print_usage("Cannot read file: $file\n"); }

my $cmd = dirname(__FILE__)."/../bin/ctx31";

if(!(-e $cmd)) {
  print_usage("Executable bin/ctx31 doesn't exist -- did you compile?");
} elsif(!(-x $cmd)) {
  print_usage("bin/ctx31 doesn't appear to be executable");
}

my $cmdline = "$cmd view --print_kmers $file";
my ($pid, $in, $out);

my ($num_nodes) = (0);
my @node_degree = (0)x5; # 0..4

$pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

while(defined(my $line = <$in>))
{
  my ($kmer, $covgs, $edges, $shades, $num_cols) = parse_ctx_line($line);
  if(defined($kmer))
  {
    $num_nodes++;
    # Merge edges
    my @edges = (0)x8;
    for(my $col = 0; $col < $num_cols; $col++) {
      for(my $i = 0; $i < 8; $i++) {
        $edges[$i] |= (substr($edges->[$col],$i,1) ne '.');
      }
    }
    my $indegree = sum(@edges[0..3]);
    my $outdegree = sum(@edges[4..7]);
    $node_degree[$indegree]++;
    $node_degree[$outdegree]++;
  }
}

close($in);
close($out);

waitpid($pid, 1);

print "Number of nodes: $num_nodes\n";
map {print "  out-degree $_: $node_degree[$_]\n"} 0..4;

sub parse_ctx_line
{
  my ($line) = @_;
  chomp($line);
  my @columns = split(' ', $line);

  if($line =~ /Error/i)
  {
    print STDERR "$line\n";
    return undef;
  }

  my ($kmer, $num_cols);
  my @covgs;
  my @edges;
  my @shades;

  if($line =~ /^([acgt]+) ((?: ?\d+)+)( [a-z\.]+)+$/i)
  {
    my $covgtxt = $2;
    @covgs = split(' ', $covgtxt);
    $num_cols = scalar(@covgs);
  }
  else
  {
    print STDERR "Cannot parse line:\n";
    print STDERR "  $line\n";
    exit(-1);
  }

# <kmer> <covg0> <covg1> <edges0> <edges1> <shades0> <shades1>

  $kmer = $columns[0];
  @edges = @columns[$num_cols+1..2*$num_cols];

  if(2*$num_cols+1 == @columns) {
    @shades = () x $num_cols;
  } else {
    @shades = @columns[2*$num_cols+1..$#columns];
  }

  # print STDERR "kmer:$kmer covgs:".join(';', @covgs)."; ".
  #              "edges:".join(';', @edges)."; ".
  #              "shades:".join(';',@shades).";\n";

  return ($kmer, \@covgs, \@edges, \@shades, $num_cols);
}
