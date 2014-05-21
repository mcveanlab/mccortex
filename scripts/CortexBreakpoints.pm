package CortexBreakpoints;

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

sub chr_to_hash
{
  my @arr = ();
  for my $line (@_) {
    if($line =~ /^(.*?):(\d+)-(\d+):([+-]):(\d+)/) {
      push(@arr, {'chrom' => $1, 'start' => $2, 'end' => $3,
                  'strand' => $4, 'offset' => $5});
    } else { die("Bad chrom: $line"); }
  }
  return @arr;
}

sub next
{
  my ($self) = @_;
  my $callid;
  my @flank5p_refs = ();
  my @flank3p_refs = ();
  my @cols = ();

  my ($hdr5p,$seq5p) = $self->read_fasta();
  if(!defined($hdr5p)) { return undef; }

  my ($hdr3p,$seq3p) = $self->read_fasta();
  if(!defined($hdr3p)) { die("Missing 3p flank"); }

  my ($pathhdr,$pathseq) = $self->read_fasta();
  if(!defined($pathhdr)) { die("Missing path sequence"); }

  # >call.0.path cols=0
  if($pathhdr =~ /^>call\.(\d+)\.path cols=(\d+(?:,\d+)*)/i) {
    $callid = $1;
    @cols = split(/,/, $2);
  }
  else { die("Cannot find cols=... : $pathhdr"); }

  if($hdr5p =~ /chr=(.*(?:,.*)*)/i) { @flank5p_refs = split(/,/, $1); }
  else { die("Cannot find 5p flank chrs=... : $hdr5p"); }

  if($hdr3p =~ /chr=(.*(?:,.*)*)/i) { @flank3p_refs = split(/,/, $1); }
  else { die("Cannot find 3p flank chrs=... : $hdr3p"); }

  @flank5p_refs = chr_to_hash(@flank5p_refs);
  @flank3p_refs = chr_to_hash(@flank3p_refs);

  return ($seq5p, $seq3p, $pathseq, \@flank5p_refs, \@flank3p_refs, \@cols, $callid);
}

1;