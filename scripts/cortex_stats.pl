#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use IPC::Open2;
use List::Util qw(sum min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use CortexScripts;

sub print_usage
{
  if(@_ > 0) { print STDERR map {"Error: $_\n"} @_; }
  
  print STDERR "" .
"Usage: ./cortex_stats.pl [--covg <maxcovg>] <in.ctx>
  Prints cortex graph stats\n";

  exit(-1);
}

my $maxcovg = 10;
if(@ARGV == 3 && $ARGV[0] =~ /--covg/) {
  shift(@ARGV);
  $maxcovg = shift(@ARGV);
  if($maxcovg !~ /^\d+$/ || $maxcovg == 0) { print_usage("Invalid --covg"); }
}

if(@ARGV != 1) { print_usage(); }
my $file = shift;

my ($path) = ($file =~ /^([^:]+):?/);
if(!(-r $path)) { print_usage("Cannot read file: $path [$file]\n"); }

my $cmd = dirname(__FILE__)."/../bin/mccortex31";

if(!(-e $cmd)) {
  print_usage("Executable bin/mccortex31 doesn't exist -- did you compile?");
} elsif(!(-x $cmd)) {
  print_usage("bin/mccortex31 doesn't appear to be executable");
}

my $cmdline = "$cmd view --kmers $file";
my ($pid, $in, $out);

my ($num_nodes) = (0);
my @node_degree = (0)x5; # 0..4

$pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");
my $num_cols = 0;
my @kmers_per_col = ();

my @covg_distrib = (0) x ($maxcovg+1);

while(defined(my $line = <$in>))
{
  my ($kmer, $covgs, $edges, $shades, $ncols) = parse_ctx_line($line);
  if(defined($kmer))
  {
    $num_nodes++;
    $num_cols = $ncols;

    # Merge edges
    my @medges = (0)x8;
    for(my $col = 0; $col < $num_cols; $col++) {
      for(my $i = 0; $i < 8; $i++) {
        $medges[$i] |= (substr($edges->[$col],$i,1) ne '.');
      }
    }

    my $indegree = sum(@medges[0..3]);
    my $outdegree = sum(@medges[4..7]);
    $node_degree[$indegree]++;
    $node_degree[$outdegree]++;

    # Update kmers per colour
    for(my $i = 0; $i < $num_cols; $i++) {
      $kmers_per_col[$i] += ($covgs->[$i] > 0 || $edges->[$i] !~ /^\.+$/);
    }

    # Sum covg
    my $sumcovg = sum(map {$covgs->[$_]} 0..($num_cols-1));
    $covg_distrib[min($maxcovg,$sumcovg)]++;
  }
}

close($in);
close($out);

waitpid($pid, 1);

print "Number of nodes: ".num2str($num_nodes)."\n";
print "Number of colours: ".num2str($num_cols)."\n";
print "Node Out-degree:\n";
map {print " $_: ".num2str($node_degree[$_])."\n"} 0..4;
# Print num of kmers per colour
print "Kmers per colour:\n";
for(my $i = 0; $i < $num_cols; $i++) {
  print " Colour $i: ".num2str($kmers_per_col[$i])."\n";
}

print "Coverage:\n";
for(my $i = 1, my $s = 0; $i <= $maxcovg && $s < $num_nodes; $i++) {
  print " $i: $covg_distrib[$i]\n";
  $s += $covg_distrib[$i];
}

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
