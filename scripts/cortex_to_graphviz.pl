#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use IPC::Open2;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexGraph;

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./cortex_to_graphviz.pl [--point|--shaded|--simplify] <in.ctx>
  Prints graphviz `dot' output.  Not to be used with large graphs!

  --point    Don't print kmer values, only points
  --shaded   Print shades
  --simplify Simplify supernodes

  Example: ./cortex_to_graphviz.pl small.ctx > small.dot
           dot -Tpng small.dot > small.png\n";

  exit(-1);
}

my $use_points = 0;
my $print_shades = 0;
my $simplify = 0;

while(@ARGV > 1 && $ARGV[0] =~ /^-/) {
  if($ARGV[0] =~ /^-?-p(oints?)?$/i) {
    shift;
    $use_points = 1;
  }
  elsif($ARGV[0] =~ /^-?-shade[sd]??$/i) {
    shift;
    $print_shades = 1;
  }
  elsif($ARGV[0] =~ /^-?-simplify$/i) {
    shift;
    $simplify = 1;
  }
  else { print_usage("Unknown option '$ARGV[0]'"); }
}

if(@ARGV > 1) { print_usage(); }

my $file = shift;
# DEV: if not defined or eq "-", use STDIN

if(!(-r $file))
{
  print STDERR "Error: Cannot read file: $file\n";
}

my $cmd = dirname(__FILE__)."/../bin/ctx31";

if(!(-e $cmd))
{
  print STDERR "Error: executable bin/ctx31 doesn't exist -- " .
               "did you compile?\n";
}
elsif(!(-x $cmd))
{
  print STDERR "Error: bin/ctx31 doesn't appear to be executable\n";
  exit(-1);
}

# colours for shades
my @cols = qw(red green blue orange purple pink brown black);

# Print a warning only once if shades mismatch
my $shade_mismatch = 0;

# cortex_bin_reader command
my $cmdline = "$cmd view --kmers $file";
my ($pid, $in, $out);

print "digraph G {\n";
print "  edge [dir=both arrowhead=none arrowtail=none]\n";
print "  node [".($use_points ? "shape=point label=none" : "shape=none")." ".
      "fontname=Courier fontsize=9]\n";

if($simplify)
{
  $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

  my $graph = new CortexGraph();
  my $num_cols;

  # Construct graph
  while(defined(my $line = <$in>))
  {
    my ($kmer, $covgs, $edges, $shades);
    ($kmer, $covgs, $edges, $shades, $num_cols) = parse_ctx_line($line);

    $graph->add_kmer($kmer);
    $graph->{$kmer}->{'shades'} = $shades;

    for(my $i = 0; $i < $num_cols; $i++)
    {
      my @edges_arr = split('', uc($edges->[$i]));

      for my $prev_edge (grep {$_ ne '.'} @edges_arr[0..3]) {
        $graph->add_edges_between(uc($prev_edge).substr($kmer,0,-1), $kmer);
      }

      for my $next_edge (grep {$_ ne '.'} @edges_arr[4..7]) {
        $graph->add_edges_between($kmer, substr($kmer,1).$next_edge);
      }
    }
  }

  # $graph->dump();
  # exit;

  # Get kmer size
  my $kmer_size = $graph->get_kmer_size();

  # print "kmer size: $kmer_size\n";

  # Simplify graph into supernodes
  # Hash of edge kmers -> supernodes
  my %super_graph = ();
  my @supernodes = ();

  for my $key (keys %$graph) {
    if(!defined($graph->{$key}->{'visited'})) {
      my $contig = $graph->get_supernode($key);
      $graph->mark_kmers_visited($contig);
      my $supernode = {'seq' => $contig};
      push(@supernodes, $supernode);
      my $key0 = kmer_key(substr($contig, 0, $kmer_size));
      my $key1 = kmer_key(substr($contig, -$kmer_size));
      $super_graph{$key0} = $supernode;
      $super_graph{$key1} = $supernode;
    }
  }

  # Print nodes
  if($print_shades)
  {
    for my $supernode (@supernodes)
    {
      my $seq = $supernode->{'seq'};
      my $kmer0 = substr($seq, 0, $kmer_size);
      my $shades0 = $graph->{kmer_key($kmer0)}->{'shades'};

      my $kmer_len = length($seq)+1-$kmer_size;
      my $num_shaded_nodes = $kmer_len < 3 ? $kmer_len : 3;
      my $num_of_shades = defined($shades0) ? length($shades0->[0]) : 0;
      my $num_of_columns = $num_shaded_nodes*$num_of_shades+($num_shaded_nodes-1);

      print $seq . ' [shape=none label=<<table ' .
  'border="'.(defined($shades0) && $shades0 =~ /^\-+$/ ? '1' : '0').'" '.
  'cellborder="0" cellpadding="0" cellspacing="0">
  <tr><td PORT="top" colspan="'.$num_of_columns.'" cellpadding="0" cellspacing="0" border="0">
  '.($use_points ? '.' : $seq).'</td>
  </tr>';

      for(my $i = 0; $i < $num_cols; $i++)
      {
        # print first kmer shades
        print '<tr>';
        print_kmer_shades($shades0->[$i]);

        if($num_shaded_nodes == 3) {
          # Get middle shades
          my $kmer1 = substr($seq,1,$kmer_size);
          my $shades1 = $graph->{kmer_key($kmer1)}->{'shades'};

          for(my $j = 2; $j < length($seq)-$kmer_size; $j++)
          {
            $kmer1 = substr($seq,$j,$kmer_size);
            my $tmp_shades = $graph->{kmer_key($kmer1)}->{'shades'};

            # Merge if not equal
            if($tmp_shades->[$i] ne $shades1->[$i])
            {
              die("Shades mismatch within a supernode");
              # if(!$shade_mismatch) {
              #   # warn("Shades mismatch within a supernode");
              #   $shade_mismatch = 1;
              # }

              # # Merge
              # for(my $k = 0; $k < length($shades1); $k++) {
              #   my ($a,$b) = map {substr($_,$j,1)} ($shades1,$tmp_shades);
              #   my ($uc,$lc) = (0,0);
              #   if($a ne '.') {
              #     if($a eq uc($a)) {$uc = $a;}
              #     if($a eq lc($a)) {$lc = $a;}
              #   }
              #   if($b ne '.') {
              #     if($b eq uc($b)) {$uc = $b;}
              #     if($b eq lc($b)) {$lc = $b;}
              #   }
              #   my $c = '-';
              #   if($uc && $lc) { $c = '-'; }
              #   elsif($uc) { $c = $uc; }
              #   elsif($lc) { $c = $lc; }
              #   substr($shades1,$j,1) = $c;
              # }
            }
          }

          # Print middle kmer shades
          print '<td>|</td>'."\n";
          print_kmer_shades($shades1->[$i]);
        }

        if($num_shaded_nodes > 1) {
          # Print last kmer shades
          my $kmer2 = substr($seq,-$kmer_size);
          my $shades2 = $graph->{kmer_key($kmer2)}->{'shades'};
          print '<td>|</td>'."\n";
          print_kmer_shades($shades2->[$i]);
        }

        print "</tr>\n";
      }

      print "</table>>];\n";
    }
  }
  else
  {
    for my $supernode (@supernodes) {
      print "  $supernode->{'seq'}\n";
    }
  }

  # Print edges
  for my $supernode (@supernodes)
  {
    my $kmer0 = substr($supernode->{'seq'}, 0, $kmer_size);
    my $kmer1 = substr($supernode->{'seq'}, -$kmer_size);
    my ($key0, $key1) = map {kmer_key($_)} ($kmer0, $kmer1);
    my $reverse0 = get_orientation($kmer0, $key0);
    my $reverse1 = get_orientation($kmer1, $key1);

    my @prev_edges = $graph->get_edges($key0,!$reverse0);
    my @next_edges = $graph->get_edges($key1,$reverse1);

    # print "@prev_edges:$kmer0  $kmer1:@next_edges\n";

    for my $next (@next_edges) {
      my $kmer = substr($supernode->{'seq'},-$kmer_size+1).$next;
      my $next_supernode = $super_graph{kmer_key($kmer)};
      print_supernode($supernode, $next, $next_supernode, $kmer, 1);
    }

    for my $prev (@prev_edges) {
      my $kmer = revcmp($prev).substr($supernode->{'seq'},0,$kmer_size-1);
      my $prev_supernode = $super_graph{kmer_key($kmer)};
      print_supernode($supernode, $prev, $prev_supernode, $kmer, 0);
    }
  }
}
else
{
  # Not 'simplifying' contigs
  if($print_shades)
  {
    $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

    while(defined(my $line = <$in>))
    {
      my ($kmer, $covgs, $edges, $shades, $num_cols) = parse_ctx_line($line);
      if(defined($kmer))
      {
        my $num_shades = defined($shades) ? length($shades->[0]) : 1;
        print $kmer . ' [shape=none label=<<table ' .
              'border="'.(defined($shades) && $shades =~ /^\-+$/ ? '1' : '0').'" '.
              'cellborder="0">
<tr><td PORT="top" colspan="'.$num_shades.'" cellpadding="0" cellspacing="0">
'.($use_points ? '.' : $kmer)."</td></tr>\n";

        for(my $i = 0; $i < $num_cols; $i++)
        {
          print '<tr>';
          print_kmer_shades($shades->[$i]);
          print "</tr>\n";
        }

        print "</table>>];\n";
      }
    }

    close($in);
    close($out);
  }

  $pid = open2($in, $out, $cmdline) or die("Cannot run cmd: '$cmdline'");

  while(defined(my $line = <$in>))
  {
    my ($kmer, $covgs, $edges, $shades, $num_cols) = parse_ctx_line($line);
    if(defined($kmer))
    {
      my $num_edges_printed = 0;

      for(my $i = 0; $i < $num_cols; $i++)
      {
        for(my $base = 0; $base < 4; $base++)
        {
          if((my $edge = substr($edges->[$i], $base, 1)) ne ".")
          {
            my $prev_kmer = uc($edge) . substr($kmer,0,-1);
            my $right_base = substr($kmer,-1);
            dump_edge($prev_kmer, $right_base, 0);
            $num_edges_printed++;
          }
        }

        for(my $base = 4; $base < 8; $base++)
        {
          if((my $edge = substr($edges->[$i], $base, 1)) ne ".")
          {
            dump_edge($kmer, uc($edge), 1);
            $num_edges_printed++;
          }
        }

        if($num_edges_printed == 0)
        {
          print "  ".kmer_key($kmer)."\n";
        }
      }
    }
  }
}

close($in);
close($out);

print "}\n";

waitpid($pid, 1);


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

sub print_kmer_shades
{
  my ($shades_txt) = @_;

  my @sharr = defined($shades_txt) ? split('', $shades_txt) : ('.');

  for(my $i = 0; $i < @sharr; $i++) {
    print '<td fixedsize="true" width="3" height="3" ' .
          'cellpadding="0" cellspacing="0" border="1" ';
    if($sharr[$i] ne '.') {
      if($sharr[$i] eq '-') { print 'bgcolor="black"'; }
      elsif($sharr[$i] eq lc($sharr[$i])) { print 'style="rounded"'; }
      else { print 'style="rounded" bgcolor="'.$cols[$i % @cols].'"'; }
      print ' color="'.$cols[$i % @cols].'"';
    }
    else { print 'color="white"'; }
    print '></td>'."\n";
  }
}

sub print_supernode
{
  my ($supernode,$rbase,$next_supernode,$kmer0,$going_right) = @_;

  my $kmer_size = length($kmer0);

  my $kmer1a = substr($next_supernode->{'seq'}, 0, $kmer_size);
  my $kmer1b = substr($next_supernode->{'seq'}, -$kmer_size);
  my ($key0,$key1a,$key1b) = map {kmer_key($_)} ($kmer0, $kmer1a, $kmer1b);

  my ($kmer1, $key1);
  if($key0 eq $key1a) { $kmer1 = $kmer1a; $key1 = $key1a; }
  elsif($key0 eq $key1b) { $kmer1 = $kmer1b; $key1 = $key1b; }
  else { die("Mismatch in supernode edges"); }

  my $arrive_left = $kmer0 eq ($going_right ? $kmer1a : revcmp($kmer1a));

  if(($supernode->{'seq'} lt $next_supernode->{'seq'}) ||
     ($supernode->{'seq'} le $next_supernode->{'seq'} &&
      ($arrive_left != $going_right || $arrive_left && $going_right)))
  {
    my $from = $supernode->{'seq'};
    my $to = $next_supernode->{'seq'};
    if($print_shades) { ($from,$to) = map {"$_:top"} ($from,$to); }
    print "  $from:" . ($going_right ? 'e' : 'w') . " -> " .
          "$to:"  . ($arrive_left ? 'w' : 'e') . "\n";
  }
}

sub dump_edge
{
  my ($kmer1, $rbase, $going_right) = @_;

  my $kmer2 = substr($kmer1, 1) . $rbase;

  my $key1 = kmer_key($kmer1);
  my $key2 = kmer_key($kmer2);

  my $rev1 = ($kmer1 ne $key1);
  my $rev2 = ($kmer2 ne $key2);

  # When doing right hand edges, do only those that go to left
  #                              or those to a greater key
  # When doing left hand edges, do only those that go to a greater key
  if(($going_right && $rev1 == $rev2) || ($rev1 != $rev2 && $key1 le $key2))
  {
    # Print a coloured edge for each colour that is in both nodes
    $key1 = $print_shades ? "$key1:top:" : "$key1:";
    $key2 = $print_shades ? "$key2:top:" : "$key2:";
    print "  $key1" . ($rev1 ? 'w' : 'e') . " -> " .
             $key2  . ($rev2 ? 'e' : 'w') . "\n";
  }
}
