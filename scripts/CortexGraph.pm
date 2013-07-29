package CortexGraph;

use strict;
use warnings;
use Carp;

my %complement = ('a' => 't', 'A' => 'T',
                  'c' => 'g', 'C' => 'G',
                  'g' => 'c', 'G' => 'C',
                  't' => 'a', 'T' => 'A');

use base 'Exporter';
our @EXPORT = qw(revcmp kmer_key get_orientation);

sub new
{
  my ($class) = @_;

  my $self = {};

  bless $self, $class;
  return $self;
}

sub revcmp
{
  my ($seq) = @_;
  for(my $i = 0; $i < length($seq); $i++)
  {
    my $b = substr($seq, $i, 1);
    my $x = $complement{$b};
    substr($seq, $i, 1) = defined($x) ? $x : $b;
  }
  my $rev = reverse($seq);
  return $rev;
}

sub kmer_key
{
  my ($kmer) = @_;
  my $kmer_revcmp = revcmp($kmer);
  return $kmer lt $kmer_revcmp ? $kmer : $kmer_revcmp;
}

sub get_orientation
{
  my ($kmer,$key) = @_;
  return $kmer eq $key ? 0 : 1;
}

sub add_edges_between
{
  my ($self, $kmer1, $kmer2) = @_;
  my ($base2, $base1) = ($complement{substr($kmer1,0,1)}, substr($kmer2,-1));
  my ($key1, $key2) = map{kmer_key($_)} ($kmer1, $kmer2);
  my $reverse1 = get_orientation($kmer1,$key1);
  my $reverse2 = get_orientation($kmer2,$key2);
  $self->{$key1}->{ $reverse1 ? 'prev_'.$base1 : 'next_'.$base1} = 1;
  $self->{$key2}->{!$reverse2 ? 'prev_'.$base2 : 'next_'.$base2} = 1;
}

sub add_kmer
{
  my ($self) = @_;
  for(my $i = 1; $i < @_; $i++) {
    my $key = kmer_key($_[$i]);
    if(!defined($self->{$key})) { $self->{$key} = {}; }
  }
}

sub get_edges
{
  my ($self,$key,$reverse) = @_;
  my @k = keys %{$self->{$key}};
  my $lookup = $reverse ? 'prev_' : 'next_';
  my @bases = grep {defined($self->{$key}->{$lookup.$_})} qw(A C G T);
  return @bases;
}

sub get_kmer_size
{
  my ($self) = @_;
  my ($key) = keys(%$self);
  return defined($key) ? length($key) : 0;
}


sub get_supernode_extension
{
  my ($self,$kmer) = @_;
  my $first_key = kmer_key($kmer);

  my $key = $first_key;
  my $reverse = get_orientation($kmer,$key);
  my $supernode = '';
  my (@edges_fw, @edges_rv);

  while(scalar(@edges_fw = $self->get_edges($key,$reverse)) == 1)
  {
    my $base = $edges_fw[0];
    $kmer = substr($kmer,1).$base;
    $key = kmer_key($kmer);
    if($key eq $first_key) { last; }
    $reverse = get_orientation($kmer,$key);
    if(scalar(@edges_rv = $self->get_edges($key,!$reverse)) != 1) {last;}
    $supernode .= $base;
  }

  return $supernode;
}

sub get_supernode
{
  my ($self,$kmer) = @_;

  my $left = $self->get_supernode_extension(revcmp($kmer));
  my $right = $self->get_supernode_extension($kmer);
  my $contig = revcmp($left).$kmer.$right;
  # print "$kmer $contig\n";

  return kmer_key($contig);
}

sub mark_kmers_visited
{
  my ($self,$contig) = @_;
  my $kmer_size = $self->get_kmer_size();
  my $n_kmers = length($contig) + 1 - $kmer_size;
  for(my $i = 0; $i < $n_kmers; $i++) {
    my $key = kmer_key(substr($contig,$i,$kmer_size));
    $self->{$key}->{'visited'} = 1;
  }
}

sub dump_kmer
{
  my ($self,$key) = @_;
  my $k = $self->{$key};
  my @prev_edges = map {defined($k->{'prev_'.$_}) ? $_ : '.'} qw(A C G T);
  my @next_edges = map {defined($k->{'next_'.$_}) ? $_ : '.'} qw(A C G T);
  my $prev = lc(revcmp(join('', @prev_edges)));
  my $next = join('', @next_edges);
  print "$key $prev$next\n";
}

sub dump
{
  my ($self) = @_;
  for my $key (sort keys %$self) {
    $self->dump_kmer($key);
  }
}

1;
