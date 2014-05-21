package CortexBubbles;

use strict;
use warnings;
use Carp;

sub new
{
  my ($class,$handle,$path) = @_;

  my $next_line;
  my @header = ();
  while(defined($next_line = <$handle>)) {
    chomp($next_line);
    if($next_line =~ /^$/) {}
    elsif($next_line =~ /^#/) { push(@header, $next_line); }
    else { last; }
  }

  my $self = {
      _handle => $handle,
      _path => $path,
      _header => \@header,
      _next_line => $next_line
  };

  bless $self, $class;
  return $self;
}

sub read_line
{
  my ($self) = @_;
  my $tmp_line = $self->{_next_line};
  my $handle = $self->{_handle};
  $self->{_next_line} = <$handle>;
  return $tmp_line;
}

sub read_fasta
{
  my ($self) = @_;
  my ($header, $seq);
  # Read past empy lines
  while(defined($header = $self->read_line()) && $header =~ /^$/) {}
  # No more entries
  if(!defined($header)) { return undef; }
  chomp($header);
  if($header !~ /^>/) { die("Bad line: $header"); }
  if(!defined($seq = $self->read_line())) { die("Missing seq line, truncated?"); }
  chomp($seq);
  return ($header, $seq);
}

sub next
{
  my ($self) = @_;
  my ($callid,$flank5p_nkmers,$flank3p_nkmers);
  my @branches = ();
  my @branchlens = ();

  my ($hdr5p,$seq5p) = $self->read_fasta();
  if(!defined($hdr5p)) { return undef; }

  my ($hdr3p,$seq3p) = $self->read_fasta();
  if(!defined($hdr3p)) { die("Missing 3p flank"); }

  # >bubble.0.5pflank kmers=8
  if($hdr5p =~ /^>bubble\.(\d+)\.5pflank kmers=(\d+)/i) {
    $callid = $1;
    $flank5p_nkmers = $2;
  }
  else { die("Cannot find 5p flank: $hdr5p"); }

  if($hdr3p =~ /^>bubble\.$callid\.3pflank kmers=(\d+)/i) {
    $flank3p_nkmers = $1;
  }
  else { die("Cannot find 3p flank: $hdr3p"); }

  while(defined($self->{'_next_line'}) && $self->{'_next_line'} =~ /^>/) {
    my ($branchhdr,$branchseq) = $self->read_fasta();
    my $nkmers = 0;
    # >bubble.0.path cols=0
    if($branchhdr =~ /^>bubble\.$callid\.branch\.(\d+) kmers=(\d+)/i) { $nkmers=$2; }
    else { die("Cannot find callid : $branchhdr"); }
    push(@branches, $branchseq);
    push(@branchlens, $nkmers);
  }

  if(@branches == 0) { die("No branches"); }

  return ($seq5p, $seq3p, \@branches,
          $flank5p_nkmers, $flank3p_nkmers, \@branchlens, $callid);
}

1;