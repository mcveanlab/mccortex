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

use constant {CMD_LIST => 1, CMD_PLOT => 2, CMD_CLEAN => 3};
use constant {EMIT_NODES_ALL => 0, EMIT_NODES_LEAVES => 2, EMIT_NODES_INPUT => 3};

# bioinf-perl modules
use GeneticsModule;
use UsefulModule;

# McCortex modules
use CortexScripts;
use CortexLinks;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./cortex_links.pl [cmd] [args]
  Plot, list, clean graph path trees. Only works on colour 0 atm.

  Commands:
    plot   [options] <file> <min_coverage>

      Produce a Graphviz (DOT) plot of the links
        --keep-subsets    Don't throw out redundant substrings

    list   [options] <file> <eff_covg.csv> <links.csv>

      Print all branches as CSV entries

    clean [options] <file> <min_coverage>

      Clean errors from the link file. input.csv is of form: <len>,<cutoff>
        --keep-subsets    Don't throw out redundant substrings

  General Options:
    --limit <N>       Stop after print <N> kmers [default: 0 (off)]
";
  exit(-1);
}

# Save args for later
my @CMDARGS = ($0, @ARGV);

if(@ARGV < 3) { print_usage(); }
my $action = 0;
my $reqargs = 0;
my $cmd = shift(@ARGV);

if($cmd eq "plot")      { $reqargs = 2; $action = CMD_PLOT;  }
elsif($cmd eq "list")   { $reqargs = 3; $action = CMD_LIST;  }
elsif($cmd eq "clean")  { $reqargs = 2; $action = CMD_CLEAN; }
else { print_usage("Bad command: $cmd"); }

# Options
my $emit_nodes_mode = ($action == CMD_LIST ? EMIT_NODES_ALL : EMIT_NODES_LEAVES);
my $use_err_model = 0;
my $limit = 0;

while(@ARGV > $reqargs) {
  my $arg = shift(@ARGV);
  if($arg eq "--keep-subsets") {
    if($action == CMD_LIST) { print_usage("$arg not for use with `list`"); }
    $action = EMIT_NODES_INPUT;
  }
  elsif($arg eq "--use-err-model") { $use_err_model    = 1; }
  elsif($arg eq "--limit") {
    $limit = shift(@ARGV);
    if($limit !~ /^\d+$/) { print_usage("Invalid --limit <N> arg: $limit"); }
  }
  else { print_usage("Bad arg: $arg"); }
}

if(@ARGV != $reqargs) { print_usage(); }

my ($ctp_path,$min_coverage);
my ($eff_covg_path,$link_csv_path);
my ($eff_covg_fh,$link_csv_fh);

if($action == CMD_PLOT) {
  ($ctp_path,$min_coverage) = @ARGV;
}
elsif($action == CMD_LIST) {
  ($ctp_path,$eff_covg_path,$link_csv_path) = @ARGV;
  open($eff_covg_fh,">$eff_covg_path") or die("Cannot open $eff_covg_path");
  open($link_csv_fh,">$link_csv_path") or die("Cannot open $link_csv_path");
}
elsif($action == CMD_CLEAN) {
  ($ctp_path,$min_coverage) = @ARGV;
}

if(defined($min_coverage)) { print STDERR "Min coverage $min_coverage\n"; }

my $ctp_fh = open_file($ctp_path);
my $ctp_file = new CortexLinks($ctp_fh, $ctp_path);

my ($kmer, @links_read) = $ctp_file->next();

if(!defined($kmer)) {
  print STDERR "No links in .ctp file\n";
}

my $colour = 0;

# Temp file used by filter command
my $tmp_file_path;
my $tmp_fh;

print STDERR "Loading JSON header...\n";
my $hdr_json = $ctp_file->ctp_get_hdr_json();
my ($contigs_hist) = $ctp_file->ctp_get_hdr_contig_hist($colour);

# Get kmer size
my $kmer_size = $hdr_json->{'graph'}->{'kmer_size'};
if(defined($kmer) && length($kmer) != $kmer_size) {
  die("Kmer size mismatch with header (".length($kmer)." != ".$kmer_size.")");
}

if($action == CMD_CLEAN)
{
  # Open a temporary file for read/write
  for(my $i = 0; $i < 20; $i++) {
    $tmp_file_path = "tmp.clean.".int(rand(10000)).".ctp";
    if(!(-e $tmp_file_path)) { last; }
  }
  if(-e $tmp_file_path) { die("Cannot create random temp file."); }

  # Create and immediately remove a temporary file
  # File will exists only to us until we close
  print STDERR "Opening temporary file: $tmp_file_path\n";
  open($tmp_fh, "+>$tmp_file_path") or die("Cannot open tmp file: $tmp_file_path");
  unlink($tmp_file_path) or warn("Cannot delete temporary file $tmp_file_path");
}

if($limit > 0) { print STDERR "Printing only the first $limit kmers\n"; }

if($action == CMD_CLEAN)   { print STDERR "Cleaning links...\n"; }
elsif($action == CMD_PLOT) { print STDERR "Plotting links...\n"; }
elsif($action == CMD_LIST)
{
  print STDERR "Calculating likelihoods...\n";

  my @eff_covg_hist = calc_eff_covg_hist(@$contigs_hist);

  print STDERR "Saving read segment histogram to $eff_covg_path...\n";
  print $eff_covg_fh "ReadLength,Counts\n";
  for(my $i = $kmer_size; $i < @eff_covg_hist; $i++) {
    print $eff_covg_fh "$i,$eff_covg_hist[$i]\n";
  }
  close($eff_covg_fh);
  print STDERR "Saving link CSV to $link_csv_path...\n";
  print $link_csv_fh "LinkLength,Count\n";

  print STDERR "Listing links...\n";
}

# Statistics on links printed
my $num_inithdr_link_kmers = $hdr_json->{'paths'}->{'num_kmers_with_paths'};
my $num_inithdr_links      = $hdr_json->{'paths'}->{'num_paths'};
my $num_inithdr_link_bytes = $hdr_json->{'paths'}->{'path_bytes'};

my $num_init_link_kmers = 0;
my $num_init_links      = 0;
my $num_init_link_bytes = 0;
my $num_init_link_nodes = 0;

my $num_output_link_kmers = 0;
my $num_output_links      = 0;
my $num_output_link_bytes = 0;
my $num_output_link_nodes = 0;

#
# dist of a node is the index of the junction in the link, where index 0 is the
# initial kmer where the link started
#

for(my $kmer_num = 0; ($limit == 0 || $kmer_num < $limit) && defined($kmer); $kmer_num++)
{
  my %trees = ();
  $trees{'F'} = {'A'=>undef, 'C'=>undef, 'G'=>undef, 'T'=>undef,
                 'dist' => -1, 'count' => 0, 'id' => 0,
                 'label' => $kmer, 'seq' => ''};
  $trees{'R'} = {'A'=>undef, 'C'=>undef, 'G'=>undef, 'T'=>undef,
                 'dist' => -1, 'count' => 0, 'id' => 0,
                 'label' => rev_comp($kmer), 'seq' => ''};

  my $nodeid = 1;

  for my $link (@links_read)
  {
    my ($seqstr) = ($link->{'other'} =~ /seq=([ACGT]+)/);
    my ($juncstr) = ($link->{'other'} =~ /juncpos=([^ ]+)/);
    if(!defined($seqstr) || !defined($juncstr)) { die("Cannot find seq= juncpos="); }
    my @juncs = split(',', $juncstr);
    if(length($link->{'juncs'}) != @juncs) { die("Mismatch in lengths"); }
    $nodeid = add_link_to_tree($trees{$link->{'dir'}}, $link->{'juncs'},
                               \@juncs, $seqstr,
                               $link->{'counts'}->[$colour],
                               $nodeid);
  }

    # Collect all links
  my @links_out = ();
  ctp_emit_links_root($trees{'F'},'F',\@links_out);
  ctp_emit_links_root($trees{'R'},'R',\@links_out);

  $num_init_link_kmers++;
  $num_init_links += scalar(@links_out);
  $num_init_link_nodes += tree_count_nodes_root($trees{'F'}) +
                          tree_count_nodes_root($trees{'R'});
  $num_init_link_bytes += sum(map {int(($_->{'num_juncs'} + 3) / 4)} @links_out);

  # print_links_in_ctp_format(\*STDERR, $kmer, $trees{'F'}, $trees{'R'});

  if($action == CMD_CLEAN)
  {
    # Threshold links
    threshold_tree($trees{'F'});
    threshold_tree($trees{'R'});
  }

  if(tree_has_branches($trees{'F'}) || tree_has_branches($trees{'R'}))
  {
    if($action == CMD_CLEAN)
    {
      # Now print .ctp to STDOUT
      print_links_in_ctp_format($tmp_fh, $kmer, $trees{'F'}, $trees{'R'});
    }
    elsif($action == CMD_PLOT) {
      # Print graphviz .dot format to STDOUT
      if(   tree_has_branches($trees{'F'})) { print_tree_in_dot_format($trees{'F'}); }
      elsif(tree_has_branches($trees{'R'})) { print_tree_in_dot_format($trees{'R'}); }

      print_links_in_ctp_format(\*STDERR, $kmer, $trees{'F'}, $trees{'R'});
      # Exit after printing one link
      last;
    }
    elsif($action == CMD_LIST) {
      if(tree_has_branches($trees{'F'})) { list_links($trees{'F'}); }
      if(tree_has_branches($trees{'R'})) { list_links($trees{'R'}); }
    }
  }

  ($kmer, @links_read) = $ctp_file->next();
}

close($ctp_fh);

print STDERR "-- Storage after duplicate removal --\n";
print STDERR "  Number of kmers with links: ".pretty_fraction($num_init_link_kmers,   $num_inithdr_link_kmers)."\n";
print STDERR "  Number of links:            ".pretty_fraction($num_init_links,        $num_inithdr_links)     ."\n";
print STDERR "  Number of link bytes:       ".pretty_fraction($num_init_link_bytes,   $num_inithdr_link_bytes)."\n";
print STDERR "  Number of link nodes:       ".num2str($num_init_link_nodes)."\n";

print STDERR "-- Output --\n";
print STDERR "  Number of kmers with links: ".pretty_fraction($num_output_link_kmers, $num_init_link_kmers)."\n";
print STDERR "  Number of links:            ".pretty_fraction($num_output_links,      $num_init_links)     ."\n";
print STDERR "  Number of link bytes:       ".pretty_fraction($num_output_link_bytes, $num_init_link_bytes)."\n";
print STDERR "  Number of link nodes:       ".pretty_fraction($num_output_link_nodes, $num_init_link_nodes)     ."\n";

print STDERR "-- Cleaning --\n";
# print STDERR "  Number of dropped cutoff:   ".pretty_fraction($num_dropped_cutoff,    $num_init_links)."\n";
# print STDERR "  Number of dropped unlikely: ".pretty_fraction($num_dropped_unlikely,  $num_init_links)."\n";

if($action == CMD_CLEAN)
{
  # Need to update header and merge with temporary file
  print STDERR "Writing cleaned links with new header...\n";

  # Reheader file
  my $hex = "0123456789abcdef";
  my $file_key = get_rand_hex_key(16);
  my $cmd_key  = get_rand_hex_key(8);
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
  $hdr_json->{'paths'}->{'num_kmers_with_paths'} = $num_output_link_kmers;
  $hdr_json->{'paths'}->{'num_paths'}            = $num_output_links;
  $hdr_json->{'paths'}->{'path_bytes'}           = $num_output_link_bytes;

  # Print updated JSON header
  my $new_hdr_txt = to_json($hdr_json, {utf8 => 1, pretty => 1});
  print "$new_hdr_txt\n
# This file was generated with McCortex
#   written by Isaac Turner <turner.isaac\@gmail.com>
#   url: https://github.com/mcveanlab/mccortex
#
# Comment lines begin with a # and are ignored, but must come after the header
# Format is:
#   [kmer] [num_paths] ...(ignored)
#   [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] ...(ignored)\n\n";

  # Print lines from temporary file
  my $tmp_line;
  seek($tmp_fh, 0, SEEK_SET);
  while(defined($tmp_line = <$tmp_fh>)) {
    print $tmp_line;
  }

  # Close temporary file (should disappear on close since already deleted)
  close($tmp_fh) or warn("Cannot close temporary file $tmp_file_path");
  if(-e $tmp_file_path) { warn("Temporary file still exists: $tmp_file_path"); }
}
elsif($action == CMD_LIST) {
  close($link_csv_fh);
}

print STDERR "Done.\n";
exit;


sub calc_eff_covg_hist
{
  my (@hist) = @_;
  my @eff_covg = (0) x scalar(@hist);
  for(my $i = 0; $i < @hist; $i++) {
    for(my $j = $i; $j < @hist; $j++) {
      $eff_covg[$i] += ($j-$i+1)*$hist[$j];
    }
  }
  return @eff_covg;
}


sub tree_has_branches
{
  my ($node) = @_;
  return (defined($node->{'A'}) || defined($node->{'C'}) ||
          defined($node->{'G'}) || defined($node->{'T'}));
}

# Count links resulting from a tree
sub tree_count_links
{
  my ($node) = @_;
  my $count = sum(map {defined($node->{$_}) ? tree_count_links($node->{$_}) : 0} qw(A C G T));
  return $count == 0 ? 1 : $count;
}

sub tree_count_links_root
{
  my ($tree) = @_;
  return tree_has_branches($tree) ? tree_count_links($tree) : 0;
}

sub tree_count_nodes
{
  my ($node) = @_;
  return 1 + sum(map {defined($node->{$_}) ? tree_count_nodes($node->{$_}) : 0} qw(A C G T));
}

sub tree_count_nodes_root
{
  my ($tree) = @_;
  my $count = tree_count_nodes($tree);
  return $count == 1 ? 0 : $count;
}

sub should_keep_branch
{
  my ($node,$base) = @_;
  # my $dist = $kmer_size + $node->{$base}->{'dist'} + 1;
  # my $count = $node->{$base}->{'count'};
  # return ($count > $min_coverage);
  return ($node->{$base}->{'count'} > $min_coverage);
}

sub threshold_tree
{
  my ($node) = @_;
  for my $base (qw(A C G T)) {
    if(defined($node->{$base})) {
      if(!should_keep_branch($node,$base)) {
        $node->{$base} = undef;
      } else {
        threshold_tree($node->{$base});
      }
    }
  }
}

sub emit_link_from_node
{
  my ($node) = @_;
  if($node->{'dist'} == 0) { return 0; }
  if($emit_nodes_mode == EMIT_NODES_ALL) { return 1; }
  my $nxt_counts = sum(map {defined($node->{$_}) ? $node->{$_}->{'count'} : 0} qw(A C G T));
  return ($nxt_counts == 0 ||
          ($emit_nodes_mode == EMIT_NODES_INPUT && $nxt_counts < $node->{'count'}));
}

sub list_links
{
  my ($node) = @_;

  for my $base (qw(A C G T)) {
    if(defined($node->{$base})) {
      my $dist = $kmer_size + $node->{$base}->{'dist'} + 1;
      my $count = $node->{$base}->{'count'};
      print $link_csv_fh "$dist,$count\n";
      list_links($node->{$base});
    }
  }
}

# If print_leaves_only, only print links that end at leaf nodes.
#   Otherwise print all original links.
# Adds links to $linksarr
sub ctp_emit_links
{
  my ($node, $dir, $juncs, $seq, $juncstr, $linksarr) = @_;

  $juncstr .= ($juncstr eq "" ? "" : ",").$node->{'dist'};
  $juncs .= $node->{'label'};
  $seq .= $node->{'label'}.$node->{'seq'};

  if(emit_link_from_node($node))
  {
    # Link ends here -> output it
    my $link = ctp_create_link($dir, $node->{'dist'}+2, length($juncs),
                               [$node->{'count'}], $juncs,
                               "seq=$seq juncpos=$juncstr");
    push(@$linksarr, $link);
  }

  for my $base (qw(A C G T)) {
    if(defined($node->{$base})) {
      ctp_emit_links($node->{$base}, $dir, $juncs, $seq, $juncstr, $linksarr);
    }
  }
}

sub ctp_emit_links_root
{
  my ($tree,$dir,$linksarr) = @_;
  for my $base (qw{A C G T}) {
    if(defined($tree->{$base})) {
      ctp_emit_links($tree->{$base}, $dir, "", $tree->{'label'}.$tree->{'seq'},
                     "", $linksarr);
    }
  }
}

# Returns number of links printed
sub print_links_in_ctp_format
{
  my ($fh,$kmer,$tree_fw,$tree_rv) = @_;

  # Collect all links
  my @links = ();
  ctp_emit_links_root($tree_fw,'F',\@links);
  ctp_emit_links_root($tree_rv,'R',\@links);

  # Print
  if(@links > 0) {
    print $fh "$kmer ".scalar(@links)."\n";
    for my $link (@links) { print $fh ctp_link_to_str($link); }

    # Update stats
    $num_output_link_kmers++;
    $num_output_links += scalar(@links);
    $num_output_link_nodes += tree_count_nodes_root($tree_fw) +
                              tree_count_nodes_root($tree_rv);
    $num_output_link_bytes += sum(map {int(($_->{'num_juncs'} + 3) / 4)} @links);
    # ^ 4 bases per byte, round up
  }

  return scalar(@links);
}

sub dot_print_node
{
  my ($node) = @_;
  my $label = $node->{'label'} . ' ' . $node->{'seq'} .
              " (dist: " . ($node->{'dist'}+1) . " count: $node->{'count'})";
  print "  node" . $node->{'id'} . " [label=\"$label\"]\n";
  for my $base (qw{A C G T}) {
    if(defined($node->{$base})) { dot_print_node($node->{$base}); }
  }
}

sub dot_print_node_edges
{
  my ($node) = @_;
  for my $base (qw{A C G T}) {
    if(defined($node->{$base})) {
      print "  node".$node->{'id'}." -> node".$node->{$base}->{'id'}.
            (should_keep_branch($node,$base) ? "" : " [color=\"red\"]")."\n";
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
    ($node->{substr($juncs,1,1)},$nodeid)
      = create_tree_node(substr($juncs,1),\@tmp,$seq,$count,$nodeid);
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
      ($node->{$base}, $nodeid) = create_tree_node(substr($juncs,$i), \@tmp,
                                                   $seq, $count, $nodeid);

      if($i > 0) {
        $node->{'seq'} = substr($seq, length($kmer)+$dists->[$i-1]+1,
                                      $dists->[$i]-$dists->[$i-1]-1);
      } else {
        $node->{'seq'} = substr($seq, length($kmer), $dists->[$i]);
      }
      last;
    }
  }

  return $nodeid;
}
