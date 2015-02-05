package CortexScripts;

use strict;
use warnings;
use Carp;

use base 'Exporter';
our @EXPORT = qw(open_file get_rand_hex_key load_json_hdr);

sub open_file
{
  my ($file) = @_;
  if(!defined($file)) { die("No file specified to open"); }
  my $handle;
  # -p checks if connected to a pipe
  if($file ne "-") {
    open($handle, $file) or die("Cannot open file: $file\n");
  }
  elsif(-p STDIN) {
    open($handle, "<&=STDIN") or die("Cannot read pipe");
  }
  else { die("Must specify or pipe in a file"); }
  return $handle;
}

sub get_rand_hex_key
{
  my ($length) = @_;
  my $key = "";
  my $hex = "0123456789abcdef";
  map {$key .= substr($hex, rand(16), 1)} 1..$length;
  return $key;
}

sub load_json_hdr
{
  my ($fh,$path) = @_;

  # Read JSON header
  my $line = <$fh>;
  if(!defined($line)) { croak("Empty file: $path"); }
  if(length($line) == 0 || substr($line, 0, 1) ne '{') { croak("Bad JSON: $path"); }
  my $header = $line;

  my ($i, $prev_offset) = (0,0);
  my ($num_curly_open, $num_curly_close) = (0,0); # '{' and '}'
  my ($num_brkt_open, $num_brkt_close) = (0,0); # '[' and ']'
  my ($in_string, $escape_char) = (0,0); # '\'

  while(1)
  {
    my $len = length($header);

    for($i = $prev_offset; $i < $len; $i++) {
      my $c = substr($header, $i, 1);
      if($in_string) {
        if($escape_char)  { $escape_char = 0; }
        elsif($c eq '\\') { $escape_char = 1; }
        elsif($c eq '"')  { $in_string = 0;   }
      }
      elsif($c eq '"') { $in_string = 1;     }
      elsif($c eq '{') { $num_curly_open++;  }
      elsif($c eq '}') { $num_curly_close++; }
      elsif($c eq '[') { $num_brkt_open++;   }
      elsif($c eq ']') { $num_brkt_close++;  }
    }
    $prev_offset = $len;

    if($num_curly_open == $num_curly_close && $num_brkt_open == $num_brkt_close) {
      last;
    }

    # header is not finished yet
    # Do some checks
    if($num_curly_open <= $num_curly_close) { croak("'}' before '{': $path"); }
    if($num_brkt_open  <  $num_brkt_close) { croak("']' before '[': $path"); }
    if($len >= 50000) { croak("Large JSON header: $path"); }

    # Read next line
    if(!defined($line = <$fh>)) { croak("Premature end of JSON: $path"); }
    $header .= $line;
  }

  return $header;
}

1;
