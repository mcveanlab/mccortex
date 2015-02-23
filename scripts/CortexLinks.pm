package CortexLinks;

use strict;
use warnings;
use Carp;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use List::Util qw(min max);
use JSON;
use CortexScripts; # load_json_hdr()

use base 'Exporter';
our @EXPORT = qw(ctp_print_link ctp_link_to_str ctp_create_link
                 ctp_get_hdr_txt ctp_get_hdr_json ctp_get_hdr_contig_hist);

sub new
{
  my ($class,$fh,$path) = @_;

  my $hdr_txt = load_json_hdr($fh,$path);
  my $hdr_json = decode_json($hdr_txt);

  my $next_line = <$fh>;

  my $self = {
      _fh => $fh,
      _path => $path,
      _hdr_txt => $hdr_txt,
      _hdr_json => $hdr_json,
      _next_line => $next_line
  };

  bless $self, $class;
  return $self;
}

sub ctp_get_hdr_txt
{
  my ($self) = @_;
  return $self->{'_hdr_txt'};
}

sub ctp_get_hdr_json
{
  my ($self) = @_;
  return $self->{'_hdr_json'};
}

sub ctp_get_hdr_contig_hist
{
  my ($self,$colour) = (@_,0); # colour defaults to zero

  my $json = $self->{'_hdr_json'};
  my $length_arr = $json->{'paths'}->{'contig_hists'}->[$colour]->{'lengths'};
  my $counts_arr = $json->{'paths'}->{'contig_hists'}->[$colour]->{'counts'};

  if(!defined($length_arr)) { die("Missing lengths array [$colour]"); }
  if(!defined($counts_arr)) { die("Missing counts array [$colour]"); }

  if(@$length_arr != @$counts_arr) { die("Mismatch [$colour]: @$length_arr != @$counts_arr"); }

  my $max_len = max(@$length_arr);
  my @contigs = map {0} 0..$max_len;

  for(my $i = 0; $i < @$length_arr; $i++) {
    $contigs[$length_arr->[$i]] += $counts_arr->[$i];
  }

  return \@contigs;
}

sub read_line
{
  my ($self) = @_;
  my $tmp_line = $self->{_next_line};
  my $fh = $self->{_fh};
  $self->{_next_line} = <$fh>;
  return $tmp_line;
}

sub ctp_create_link
{
  my ($dir,$nkmers,$njuncs,$countsarr,$juncstr,$other) = @_;
  return {'dir' => $dir, 'num_kmers' => $nkmers, 'num_juncs' => $njuncs,
          'counts' => $countsarr, 'juncs' => $juncstr, 'other' => $other};
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
      my ($dir,$nkmers,$njuncs,$counts_str,$juncs,$other)
        = ($line =~ /^([FR]) (\d+) (\d+) (\d+(?:,\d+)*) ([ACGT]+)(?: ?(.*))/);

      if(!defined($dir)) { die("Quak: Bad line [$self->{_path}]: $line"); }

      push(@paths, {'dir'=>$dir, 'num_kmers'=>$nkmers, 'num_juncs'=>$njuncs,
                    'counts'=>[split(',',$counts_str)], 'juncs'=>$juncs,
                    'other'=>$other});
    }
  }

  return ($kmer,@paths);
}

sub ctp_link_to_str
{
  my $str = "";
  for my $p (@_) {
    my $other = ""; # Other info

    if(defined($p->{'other'}) && length($p->{'other'}) > 0) {
      $other = " ".$p->{'other'};
    }

    $str .= $p->{'dir'} . " " . $p->{'num_kmers'} . " ". $p->{'num_juncs'} . " " .
            join(',', @{$p->{'counts'}}) . " " . $p->{'juncs'} . $other . "\n";
  }
  return $str;
}

sub ctp_print_link
{
  my ($p,$out) = @_;
  my $txt = ctp_link_to_str($p);
  if(defined($out)) { print $out $txt; }
  else { print $txt; }
}

1;