#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use McCortexScripts; # load_json_hdr()

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage $0 <file>

  Read only the JSON header from a file.

";

  exit(-1);
}

if(@ARGV != 1) { print_usage(); }
my $path = shift(@ARGV);

use IO::Zlib;
my $gz = new IO::Zlib;
$gz->open($path, "rb") or die("Cannot open file: $path");

my $hdr = load_json_hdr($gz, $path);

print $hdr;

$gz->close();
