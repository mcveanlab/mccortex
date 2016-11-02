package LineReader;

use strict;
use warnings;
use Carp;

sub new
{
  my ($class,$fh,$path) = @_;

  my $self = {
      _fh => $fh,
      _path => $path,
      _next => []
  };

  bless $self, $class;
  return $self;
}

sub read_line
{
  my ($self) = @_;
  my $fh = $self->{_fh};
  my $next = shift(@{$self->{_next}});
  return defined($next) ? $next : <$fh>;
}

sub unread_line
{
  my ($self,$line) = @_;
  unshift(@{$self->{_next}}, $line);
}

1;