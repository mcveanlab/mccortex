package CortexScripts;

use base 'Exporter';
our @EXPORT = qw(open_file);

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

1;
