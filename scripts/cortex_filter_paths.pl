#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use List::Util qw(min max sum);
use Fcntl qw(SEEK_SET);
use Cwd qw(getcwd abs_path);
use Sys::Hostname qw(hostname);
use Config;
use JSON;
use POSIX;

# bioinf-perl modules
use GeneticsModule;
use UsefulModule;

# McCortex modules
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

# Save args for later
my @CMDARGS = ($0, @ARGV);

if(@ARGV < 2) { print_usage(); }
my $cmd = shift(@ARGV);
if($cmd ne "plot" && $cmd ne "filter") { print_usage("Bad command: $cmd"); }

# Options
my $print_leaves_only = 1;
my $genome_size;

while(@ARGV > 1) {
  my $arg = shift(@ARGV);
  if($arg eq "--keep-subsets") { $print_leaves_only = 0; }
  elsif($arg eq "--genome") {
    $genome_size = shift(@ARGV);
    if(!defined($genome_size)) {
      print_usage("--genome <G> requires positive number: '$genome_size'");
    }
    $genome_size = str2num($genome_size);
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
my $max_count = 100;
my $threshold = 0.001;
my $expect_tbl;
my $cutoff_dist;

my $hdr_txt = $ctp_file->ctp_get_header();
my $hdr_json;
my $tmp_file_path;
my $tmp_fh;

if($filter_paths)
{
  print STDERR "Loading JSON header...\n";
  $hdr_json = decode_json($hdr_txt);
  my ($contigs) = json_hdr_load_contig_hist($hdr_json,$colour);

  # Open a temporary file for read/write
  for(my $i = 0; $i < 20; $i++) {
    $tmp_file_path = "tmp.filter.".int(rand(10000)).".ctp";
    if(!(-e $tmp_file_path)) { last; }
  }
  if(-e $tmp_file_path) { die("Cannot create random temp file."); }

  # Create and immediately remove a temporary file
  # File will exists only to us until we close
  print STDERR "Opening temporary file: $tmp_file_path\n";
  open($tmp_fh, "+>$tmp_file_path") or die("Cannot open tmp file: $tmp_file_path");
  unlink($tmp_file_path) or warn("Cannot delete temporary file $tmp_file_path");

  print STDERR "Calculating likelihoods...\n";
  print STDERR "  (Generating ".scalar(@$contigs)." x $max_count table)\n";
  ($expect_tbl,undef) = calc_exp_runlens($genome_size, $kmer_size, $max_count, @$contigs);
  # Find length cutoff (keep links <$cutoff)
  $cutoff_dist = find_link_cutoff($expect_tbl,$kmer_size,$threshold);
  print STDERR "Cutoff is length <$cutoff_dist bp\n";
  print STDERR "Filtering paths...\n";
  # print "$hdr_txt\n\n";
  # print_expect_table($expect_tbl, $kmer_size);
}

# Format is:
#   [kmer] [num_paths] ...(ignored)
#   [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] ...(ignored)
# CACTGATGA 1
# R 5 1 0,0,0,1 G
# CAGTGGCCG 1
# R 5 1 0,0,0,1 A

my $num_kmers_with_paths = 0;
my $num_paths = 0;
my $num_path_bytes = 0;

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

  for my $path (@paths)
  {
    my ($seqstr) = ($path->{'other'} =~ /seq=([ACGT]+)/);
    my ($juncstr) = ($path->{'other'} =~ /juncpos=([^ ]+)/);
    if(!defined($seqstr) || !defined($juncstr)) { die("Cannot find seq= juncpos="); }
    my @juncs = split(',', $juncstr);
    if(length($path->{'juncs'}) != @juncs) { die("Mismatch in lengths"); }
    $nodeid = add_link_to_tree($trees{$path->{'dir'}}, $path->{'juncs'},
                               \@juncs, $seqstr,
                               $path->{'counts'}->[$colour],
                               $nodeid);
  }

  if($filter_paths)
  {
    # print STDERR "-- PREFILTER --\n";
    # print_links_in_ctp_format(\*STDERR, $kmer, $trees{'F'}, $trees{'R'}, $print_leaves_only);

    # Threshold paths
    threshold_tree($trees{'F'});
    threshold_tree($trees{'R'});
  }

  if(tree_has_branches($trees{'F'}) || tree_has_branches($trees{'R'}))
  {
    if($filter_paths)
    {
      # Now print .ctp to STDOUT
      print_links_in_ctp_format($tmp_fh, $kmer, $trees{'F'}, $trees{'R'}, $print_leaves_only);
    }
    else {
      # Print graphviz .dot format to STDOUT
      print_tree_in_dot_format($trees{'F'});
      print_tree_in_dot_format($trees{'R'});

      print_links_in_ctp_format(\*STDERR, $kmer, $trees{'F'}, $trees{'R'}, $print_leaves_only);
    }
  }
  # else {  print STDERR "  bare.\n"; }

  ($kmer, @paths) = $ctp_file->next();
}

close($ctp_fh);


if($filter_paths)
{
  # Need to update header and merge with temporary file
  print STDERR "Writing filtered paths with new header...\n";

  # Reheader file
  my $hex = "0123456789abcdef";
  my $file_key = gen_rand_key(16);
  my $cmd_key  = gen_rand_key(8);
  $hdr_json->{'file_key'} = $file_key;

  my $datetime = strftime("%Y-%m-%d %H:%M:%S", localtime(time));
  my $user = getlogin() || getpwuid($<) || "Unknowable";
  my $cwd = abs_path(getcwd());
  my $host = hostname();

  # my $os = $Config{'osname'};
  # my $osversion = $Config{'osvers'};
  # my $osrelease = $^O;
  # my $hardware = $Config{'archname'};
  my $os = `uname -s` or warn("Cannot call uname");
  my $osrelease = `uname -r` or warn("Cannot call uname");
  my $osversion = `uname -v` or warn("Cannot call uname");
  my $hardware = `uname -m` or warn("Cannot call uname");

  chomp($os);
  chomp($osrelease);
  chomp($osversion);
  chomp($hardware);

  # Add new command
  my $new_cmd = {
    'key' => $cmd_key,
    'cmd' => \@CMDARGS,
    'out_key' => $file_key,
    'out_path' => "-",
    'date' => $datetime,
    'os' => $os,
    'osversion' => $osversion,
    'osrelease' => $osrelease,
    'hardware' => $hardware,
    'cwd' => $cwd,
    'user' => $user,
    'host' => $host,
    'prev' => [$hdr_json->{'commands'}->[0]->{'key'}]
    };
  unshift(@{$hdr_json->{'commands'}}, $new_cmd);

  # Update path info
  $hdr_json->{'paths'}->{'num_kmers_with_paths'} = $num_kmers_with_paths;
  $hdr_json->{'paths'}->{'num_paths'}            = $num_paths;
  $hdr_json->{'paths'}->{'path_bytes'}           = $num_path_bytes;

  # Print updated JSON header
  my $new_hdr_txt = to_json($hdr_json, {utf8 => 1, pretty => 1});
  print "$new_hdr_txt\n\n";

  # Print lines from temporary file
  my $tmp_line;
  seek($tmp_fh, 0, SEEK_SET);
  while(defined($tmp_line = <$tmp_fh>)) {
    print $tmp_line;
  }

  # Close temporary file (should disappear on close since already deleted)
  close($tmp_fh) or warn("Cannot close temporary file $tmp_file_path");
  if(-e $tmp_file_path) { warn("Temporary file still exists: $tmp_file_path"); }

  print STDERR "Done.\n";
}

exit;


sub gen_rand_key
{
  my ($length) = @_;
  my $key = "";
  my $hex = "0123456789abcdef";
  map {$key .= substr($hex, rand(16), 1)} 1..$length;
  return $key;
}

sub tree_has_branches
{
  my ($node) = @_;
  return (defined($node->{'A'}) || defined($node->{'C'}) ||
          defined($node->{'G'}) || defined($node->{'T'}));
}

sub threshold_tree
{
  my ($node) = @_;
  for my $base (qw(A C G T)) {
    if(defined($node->{$base})) {
      my $dist = $kmer_size + $node->{$base}->{'dist'} + 1;
      my $count = $node->{$base}->{'count'};
      if(!should_keep_link($expect_tbl, $max_count, $threshold, $cutoff_dist,
                           $dist, $count))
      {
        my $exp = $expect_tbl->[$dist]->[min($max_count,$count)];
        # print STDERR "Drop dist:$dist count:$count exp:$exp [cutoff_dist: $cutoff_dist]\n";
        $node->{$base} = undef;
      }
      else {
        threshold_tree($node->{$base});
      }
    }
  }
}

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

sub ctp_emit_links_root
{
  my ($tree,$dir,$leaf_only,$linksarr) = @_;
  for my $base (qw{A C G T}) {
    if(defined($tree->{$base})) {
      ctp_emit_links($tree->{$base}, $dir, "", $tree->{'label'}.$tree->{'seq'}, "",
                     $leaf_only, $linksarr);
    }
  }
}

sub print_links_in_ctp_format
{
  my ($fh,$kmer,$tree_fw,$tree_rv,$leaf_only,$printstdout) = @_;
  my @links = ();

  # Collect all links
  ctp_emit_links_root($tree_fw,'F',$leaf_only,\@links);
  ctp_emit_links_root($tree_rv,'R',$leaf_only,\@links);

  # Print
  if(@links > 0) {
    print $fh "$kmer ".scalar(@links)."\n";
    for my $link (@links) { print $fh ctp_path_to_str($link); }

    # Update stats
    $num_kmers_with_paths++;
    $num_paths += @links;
    # 4 bases per byte, round up
    $num_path_bytes += sum(map {int(($_->{'num_juncs'} + 3) / 4)} @links);
  }
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
      if($i > 0) {
        # Add seq to new node
        $node->{'seq'} = substr($seq, length($kmer)+$dists->[$i-1]+1,
                                      $dists->[$i]-$dists->[$i-1]-1);
      }
      last;
    }
  }

  return $nodeid;
}
