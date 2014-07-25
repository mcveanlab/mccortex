#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules as well as bioinf-perl module
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../libs/bioinf-perl/lib';

use CortexBreakpoints;
use GeneticsModule;
use FASTNFile;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./cortex_brkpnts_check_ref.pl <in.txt> <ref.fa> [ref2.fa ...]
  Check breakpoint calls match the reference\n";

  exit(-1);
}

if(@ARGV < 2) { print_usage(); }

my $file = shift(@ARGV);
my @ref_paths = @ARGV;

if($file =~ /^-?-h(elp)?/) { print_usage(); }
my $fh;
open($fh, $file) or die("Cannot read file $file");

my $cb = new CortexBreakpoints($fh);
my ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid);

# Load references
my %refs = ();
for my $ref_path (@ref_paths)
{
  my $fastn = open_fastn_file($ref_path);
  my ($title,$seq);
  while((($title,$seq) = $fastn->read_next()) && defined($title)) {
    $title =~ s/^(\W+).*$/$1/g; # Strip everything after the first whitespace
    if(defined($refs{$title})) { die("Ref names clash: $title"); }
    $refs{$title} = $seq;
  }
  close_fastn_file($fastn);
}

my $num_checked = 0;

while(1)
{
  ($seq5p, $seq3p, $pathseq, $flank5p_refs, $flank3p_refs, $cols, $callid) = $cb->next();
  if(!defined($seq5p)) { last; }

  check_ref_regions($callid, $seq5p, @$flank5p_refs);
  check_ref_regions($callid, $seq3p, @$flank3p_refs);
  $num_checked++;
}

print "Checked $num_checked entries. All passed.\n";

close($fh);

sub check_ref_regions
{
  my $callid = shift(@_);
  my $seq = shift(@_);
  for my $reg (@_)
  {
    my $str = "call" . $callid . " => " . $reg->{'chrom'} . ":" .
              $reg->{'start'} . '-' . $reg->{'end'} . ':' .
              $reg->{'strand'} . ':' . $reg->{'offset'};

    my $fwstrand = ($reg->{'strand'} eq '+');
    if(($fwstrand  && $reg->{'start'} > $reg->{'end'}) ||
       (!$fwstrand && $reg->{'start'} < $reg->{'end'})) { die("BAD Strand: $str"); }

    my $len = abs($reg->{'start'} - $reg->{'end'}) + 1;
    if($reg->{'offset'}-1 + $len > length($seq)) { die("Seq too short: $str"); }

    # Check chrom exists
    if(!defined($refs{$reg->{'chrom'}})) { die("Missing chrom: $str"); }

    my $query = substr($seq, $reg->{'offset'}-1, $len);
    if(!$fwstrand) { $query = rev_comp($query); }

    my $start = min($reg->{'start'}, $reg->{'end'});
    my $ref = substr($refs{$reg->{'chrom'}}, $start-1, $len);

    if(uc($ref) ne uc($query)) {
      print STDERR "region: $str\nQ: $query\nR: $ref\n";
      die("Mismatch");
    }
  }
}
