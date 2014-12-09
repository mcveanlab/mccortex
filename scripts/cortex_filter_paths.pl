#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use List::Util qw(sum);

use GeneticsModule;

use CortexScripts;
use CortexPaths;
use CortexLinkCleaning;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./cortex_filter_paths.pl [plot|filter] [options] <in.ctp>
  Plot or filter graph path trees. Only works on colour 0 atm.

  Options:
    --keep-subsets    Don't through out redundant substrings
    --genome <GSize>  Genome size of G bases [required for filter]
";
  exit(-1);
}

if(@ARGV < 2) { print_usage(); }
my $cmd = shift;
if($cmd ne "plot" && $cmd ne "filter") { print_usage("Bad command: $cmd"); }

# Options
my $print_leaves_only = 1;
my $genome_size;

while(@ARGV > 1) {
  my $arg = shift;
  if($arg eq "--keep-subsets") { $print_leaves_only = 0; }
  elsif($arg eq "--genome") {
    $genome_size = shift;
    if(!defined($genome_size) || $genome_size !~ /^\d+$/) {
      print_usage("--genome <G> requires positive number: '$genome_size'");
    }
  }
  else { print_usage("Bad arg: $arg"); }
}

if(@ARGV != 1) { print_usage(); }
my ($ctp_path) = @ARGV;

my $filter_paths = ($cmd eq "filter");
my $ctp_fh = open_file($ctp_path);
my $ctp_file = new CortexPaths($ctp_fh, $ctp_path);

my ($kmer, @paths) = $ctp_file->next();

if(!defined($kmer)) {
  print STDERR "No links in .ctp file\n";
  exit;
}

my $colour = 0;
my $kmer_size = length($kmer);
my $max_count = 150;
my $threshold = 0.01;
my $expect_tbl;

if($filter_paths)
{
  my $hdr_txt = $ctp_file->ctp_get_header();
  print STDERR "Loading JSON header...\n";
  my ($contigs) = json_hdr_load_contig_hist($hdr_txt,$colour);
  print STDERR "Calculating likelihoods...\n";
  print STDERR "  (Generating ".scalar(@$contigs)." x $max_count table)\n";
  ($expect_tbl) = calc_exp_runlens($genome_size, $kmer_size, $max_count, @$contigs);
  # Find length cutoff (keep links <$cutoff)
  my $cutoff = find_link_cutoff($expect_tbl,$kmer_size,$threshold);
  print STDERR "Cutoff is <$cutoff\n";
  print STDERR "Filtering paths...\n";
}

# Format is:
#   [kmer] [num_paths] ...(ignored)
#   [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] ...(ignored)
# CACTGATGA 1
# R 5 1 0,0,0,1 G
# CAGTGGCCG 1
# R 5 1 0,0,0,1 A

# Each node has: {'A' => undef, 'C' => undef, 'G' => undef, 'T' => undef,
#                 'dist' => 12, 'count' => 12, 'id' = "node0"}
while(defined($kmer))
{
  my %trees = ();
  $trees{'F'} = {'A'=>undef, 'C'=>undef, 'G'=>undef, 'T'=>undef,
                 'dist' => -1, 'count' => 0, 'id' => 0,
                 'label' => $kmer, 'seq' => ''};
  $trees{'R'} = {'A'=>undef, 'C'=>undef, 'G'=>undef, 'T'=>undef,
                 'dist' => -1, 'count' => 0, 'id' => 0,
                 'label' => rev_comp($kmer), 'seq' => ''};

  my $nodeid = 1;
  # print STDERR "$kmer ".scalar(@paths)."\n";
  # print STDERR "".ctp_path_to_str(@paths)."\n";

  for my $path (@paths)
  {
    my ($seqstr) = ($path->{'other'} =~ /seq=([ACGT]+)/);
    my ($juncstr) = ($path->{'other'} =~ /juncpos=([^ ]+)/);
    if(!defined($seqstr) || !defined($juncstr)) { die("Cannot find seq= juncpos="); }
    my @juncs = split(',', $juncstr);
    if(length($path->{'juncs'}) != @juncs) { die("Mismatch in lengths"); }
    $trees{$path->{'dir'}}->{'seq'} = substr($seqstr,length($kmer),$juncs[0]);
    $nodeid = add_link_to_tree($trees{$path->{'dir'}}, $path->{'juncs'},
                               \@juncs, $seqstr,
                               $path->{'counts'}->[$colour],
                               $nodeid);
  }

  if($filter_paths)
  {
    # Threshold paths
    # TODO

    # Now print .ctp to STDOUT
    print_links_in_ctp_format($trees{'F'},"F",$print_leaves_only,1);
    print_links_in_ctp_format($trees{'R'},"R",$print_leaves_only,1);
  }
  else {
    # Print graphviz .dot format to STDOUT
    print_tree_in_dot_format($trees{'F'});
    print_tree_in_dot_format($trees{'R'});
    # Now print .ctp to STDERR
    print_links_in_ctp_format($trees{'F'},"F",$print_leaves_only,0);
    print_links_in_ctp_format($trees{'R'},"R",$print_leaves_only,0);
  }

  ($kmer, @paths) = $ctp_file->next();
}

close($ctp_fh);

# If leaf_only, only print links that end at leaf nodes.
#   Otherwise print all original links.
# Adds links to $linksarr
sub ctp_emit_links
{
  my ($node,$dir,$juncs,$seq,$juncstr,$leaf_only,$linksarr) = @_;

  $juncstr .= ($juncstr eq "" ? "" : ",").$node->{'dist'};
  $juncs .= $node->{'label'};
  $seq .= $node->{'label'}.$node->{'seq'};

  my $nxt_counts = sum(map {defined($node->{$_}) ? $node->{$_}->{'count'} : 0} qw(A C G T));

  if($nxt_counts == 0 || (!$leaf_only && $nxt_counts < $node->{'count'}))
  {
    # Link ends here -> output it
    my $link = ctp_create_path($dir, $node->{'dist'}+2, length($juncs),
                               [$node->{'count'}], $juncs,
                               "seq=$seq juncpos=$juncstr");
    push(@$linksarr, $link);
  }

  for my $base (qw(A C G T)) {
    if(defined($node->{$base})) {
      ctp_emit_links($node->{$base}, $dir, $juncs, $seq, $juncstr,$leaf_only,$linksarr);
    }
  }
}

sub print_links_in_ctp_format
{
  my ($tree,$dir,$leaf_only,$printstdout) = @_;
  my @links = ();
  for my $base (qw{A C G T}) {
    if(defined($tree->{$base})) {
      ctp_emit_links($tree->{$base}, $dir, "", $tree->{'label'}.$tree->{'seq'}, "",
                     $leaf_only, \@links);
    }
  }
  # Print
  if($printstdout) { map {print STDOUT ctp_path_to_str($_)} @links; }
  else             { map {print STDERR ctp_path_to_str($_)} @links; }
}

sub dot_print_node
{
  my ($node) = @_;
  print "  ".$node->{'id'}." [label=\"$node->{'label'} $node->{'seq'} (dist: ".($node->{'dist'}+1)." count: $node->{'count'})\"]\n";
  for my $base (qw{A C G T}) {
    if(defined($node->{$base})) { dot_print_node($node->{$base}); }
  }
}

sub dot_print_node_edges
{
  my ($node) = @_;
  for my $base (qw{A C G T}) {
    if(defined($node->{$base})) {
      print "  node".$node->{'id'}." -> node".$node->{$base}->{'id'}."\n";
      dot_print_node_edges($node->{$base});
    }
  }
}

sub print_tree_in_dot_format
{
  my ($tree) = @_;
  print "digraph G {\n";
  print "  node [shape=none fontname=\"Courier New\" fontsize=9]\n";
  dot_print_node($tree);
  dot_print_node_edges($tree);
  print "}\n";
}

sub create_tree_node
{
  my ($juncs,$dists,$seq,$count,$nodeid) = @_;

  # Trim off prev seq
  # $seq = substr($seq,$dists->[0]);

  my $subseq = "";
  if(@$dists > 1) {
    $subseq = substr($seq, length($kmer)+$dists->[0]+1, $dists->[1]-$dists->[0]-1);
  }

  my $node = {'A'=>undef, 'C'=>undef, 'G'=>undef, 'T'=>undef,
              'dist' => $dists->[0], 'count' => $count, 'id' => $nodeid,
              'label' => substr($juncs,0,1), 'seq' => $subseq};
  $nodeid++;

  if(length($juncs) > 1)
  {
    # Trim off first junction choice
    my @tmp = @$dists;
    @tmp = @tmp[1..$#tmp];
    ($node->{substr($juncs,1,1)},$nodeid) = create_tree_node(substr($juncs,1),\@tmp,$seq,$count,$nodeid);
  }

  return ($node,$nodeid);
}


sub add_link_to_tree
{
  my ($tree,$juncs,$dists,$seq,$count,$nodeid) = @_;

  # my @thing = @$dists;
  # print "$tree $juncs @thing $count\n";
  $tree->{'count'} += $count;

  # Look up current juncs
  my $node = $tree;
  my $offset = 0;
  for(my $i = 0; $i < length($juncs); $offset += $dists->[$i], $i++) {
    my $base = substr($juncs,$i,1);
    if(defined($node->{$base})) {
      $node = $node->{$base};
      $node->{'count'} += $count;
    }
    else {
      my @tmp = @$dists;
      @tmp = @tmp[$i..$#tmp];
      ($node->{$base},$nodeid) = create_tree_node(substr($juncs,$i),\@tmp,$seq,$count,$nodeid);
      last;
    }
  }

  return $nodeid;
}
