package CortexPaths;

use strict;
use warnings;
use Carp;
use CortexScripts; # load_json_hdr()

use base 'Exporter';
our @EXPORT = qw(ctp_print_path ctp_path_to_str);

sub new
{
  my ($class,$fh,$path) = @_;

  my $header = load_json_hdr($fh,$path);
  my $next_line = <$fh>;

  my $self = {
      _fh => $fh,
      _path => $path,
      _header => $header,
      _next_line => $next_line
  };

  bless $self, $class;
  return $self;
}

sub read_line
{
  my ($self) = @_;
  my $tmp_line = $self->{_next_line};
  my $fh = $self->{_fh};
  $self->{_next_line} = <$fh>;
  return $tmp_line;
}

sub next
{
  my ($self) = @_;

  # Read over empty lines and comment lines
  my $line;
  while(defined($line = $self->read_line()) && $line =~ /^\s*(?:#.*)?$/) {}
  if(!defined($line)) { return undef; } # empty file

  my ($kmer) = ($line =~ /^([ACGT]+)/);
  if(!defined(!$kmer)) { croak("Bad line [$self->{_path}]: '$line'"); }
  my @paths = ();

  while(defined($self->{_next_line}) && $self->{_next_line} !~ /^[ACGT]/) {
    $line = $self->read_line();

    # Read over empty lines and comment lines
    if($line !~ /^\s*(?:#.*)?$/)
    {
      my ($dir,$nkmers,$njuncs,$counts_str,$seq,$other)
        = ($line =~ /^([FR]) (\d+) (\d+) (\d+(?:,\d+)*) ([ACGT]+)(?: ?(.*))/);

      if(!defined($dir)) { die("Quak: Bad line [$self->{_path}]: $line"); }

      push(@paths, {'dir'=>$dir, 'num_kmers'=>$nkmers, 'num_juncs'=>$njuncs,
                    'counts'=>[split(',',$counts_str)], 'seq'=>$seq,
                    'other'=>$other});
    }
  }

  return ($kmer,@paths);
}

sub ctp_path_to_str
{
  my $str = "";
  for my $p (@_) {
    my $other = ""; # Other info

    if(defined($p->{'other'}) && length($p->{'other'}) > 0) {
      $other = " ".$p->{'other'};
    }

    $str .= $p->{'dir'} . " " . $p->{'num_kmers'} . " ". $p->{'num_juncs'} . " " .
            join(',', @{$p->{'counts'}}) . " " . $p->{'seq'} . $other . "\n";
  }
  return $str;
}

sub ctp_print_path
{
  my ($p,$out) = @_;
  my $txt = ctp_print_path_str($p);
  if(defined($out)) { print $out $txt; }
  else { print $txt; }
}

1;