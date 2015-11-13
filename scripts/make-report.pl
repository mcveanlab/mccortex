#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max min sum);

# TODO:
# [x] parse args, print usage
# [ ] Crawl directory to find:
#     [ ] list of kmers used
#     [ ] list of samples
#     [ ] vcf files produced
#     [ ] plots for each sample
# [x] print latex
#

# Get kmers:
# $dir/k*

# Get samples:
# $dir/k$k/graphs/*.raw.ctx

# Get kmer coverage:
# $dir/k$k/graphs/$sample.clean.kmercov

# Include plots:
# $dir/k$k/graphs/$sample.raw.cov.pdf
# $dir/k$k/graphs/$sample.clean.cov.pdf
# $dir/k$k/graphs/$sample.raw.len.pdf
# $dir/k$k/graphs/$sample.clean.len.pdf
# $dir/k$k/graphs/$sample.se.links.pdf
# $dir/k$k/graphs/$sample.pe.links.pdf

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage $0 [options] <dir> [srcpath]

  Generate a report in latex on a McCortex pipeline run. Prints latex to STDOUT.

  Options:
    None.

  Example:
    $0 mcrun > run_report.tex
    pdflatex run_report
";

  exit(-1);
}

if(@ARGV < 1 || @ARGV > 2) { print_usage(); }

my $dir = shift(@ARGV);
my $srcpath = shift(@ARGV);

if(!defined($srcpath)) { $srcpath = '.'; }

my @samples = ();
my @kmers = ();

# Crawl
my @files;

# find kmers
@files = ls($dir);
@kmers = map {$_ =~ /^k(\d+)$/} (grep {/^k\d+$/} ls($dir));
# find samples
@files = map {lspath("$dir/k$_/graphs/", "$dir/k$_/links/")} @kmers;
@samples = rmdups(map {$_ =~ /\/([^\/]*).raw.ctx$/} (grep {/graphs\/.*?.raw.ctx$/} @files));
my %pathhash = ();
map {$pathhash{$_} = 1} @files;

# Debug: list all paths found
# my @tmp = sort keys(%pathhash);
# print STDERR "files:".join("\n", @tmp)."\n";

print STDERR "Kmers: @kmers\n";
print STDERR "Samples: @samples\n";

# Get cleaning thresholds:
my %kmer_cleaning = (); # {$sample}->{$k} => kmer cleaning threshold
my %link_cleaning = (); # {$sample}->{$k}->{"se"|"pe"} => cleaning threshold

for my $sample (@samples) {
  $kmer_cleaning{$sample} = {};
  $link_cleaning{$sample} = {};
  for my $k (@kmers) {
    $link_cleaning{$sample}->{$k} = {};
  }
}

for my $sample (@samples) {
  for my $k (@kmers) {
    my $logpath = "$dir/k$k/graphs/$sample.clean.ctx.log";
    if(defined($pathhash{$logpath})) {
      my $thresh;
      open(FH, $logpath) or die("Cannot read $logpath: $!");
      while(defined(my $line = <FH>)) {
        if($line =~ /Removing unitigs with coverage < (\d+)/i) { $thresh = $1; last; }
      }
      close(FH);
      if(defined($thresh)) { $kmer_cleaning{$sample}->{$k} = $thresh; }
    }

    for my $type (qw(se pe)) {
      $logpath = "$dir/k$k/links/$sample.$type.thresh.txt";
      if(defined($pathhash{$logpath})) {
        my $thresh;
        open(FH, $logpath) or die("Cannot read $logpath: $!");
        while(defined(my $line = <FH>)) {
          if($line =~ /suggested_cutoff=(\d+)/i) { $thresh = $1; last; }
        }
        close(FH);
        if(defined($thresh)) { $link_cleaning{$sample}->{$k}->{$type} = $thresh; }
      }
    }
  }
}

# List of files in a directory, including path
sub lspath
{
  my @files = ();
  for my $path (@_) {
    $path =~ s/\/+$//g; # remove trailing slashes
    @files = (@files, map {$path.'/'.$_} (grep {$_ ne "." && $_ ne ".."} ls($path)));
  }
  return @files;
}

# List of files in a directory
sub ls
{
  my ($path) = @_;
  my $dh;
  opendir($dh, $path) or die "Couldn't open dir '$path': $!";
  my @files = readdir $dh;
  closedir($dh);
  return @files;
}

sub rmdups
{
  my %m = ();
  for my $x (@_) { $m{$x} = 1; }
  return keys(%m);
}

# Latex is weird, dots in filenames confuse it
# "wow.much.pdf" -> "{wow.much}.pdf"
sub latexpath
{
  my ($path) = @_;
  $path =~ s/\/([^\/]*).pdf$/\/\{$1\}.pdf/g;
  return $path;
}

sub mksubfloat
{
  my ($images, $labels, $captions, $label, $caption, $width) = @_;
  print '
\begin{figure}[hp]
\centering
  \subfloat['.$captions->[0].']{
    ';
  if(defined($pathhash{$images->[0]})) {
    print '\includegraphics[width='.$width.'\textwidth]{\srcpath/'.latexpath($images->[0]).'}';
  } else { print '\hspace{'.$width.'\textwidth}'; } print '
    \label{'.$labels->[0].'}
  }
  \hfill
  \subfloat['.$captions->[1].']{
    ';
  if(defined($pathhash{$images->[1]})) {
    print '\includegraphics[width='.$width.'\textwidth]{\srcpath/'.latexpath($images->[1]).'}';
  } else { print '\hspace{'.$width.'\textwidth}'; } print '
    \label{'.$labels->[1].'}
  }
  \caption{'.$caption.'}
  \label{'.$label.'}
\end{figure}
';
}

sub mk_cov_fig
{
  my ($k,$sample) = @_;
  my @images = ("$dir/k$k/graphs/$sample.raw.cov.pdf",
                "$dir/k$k/graphs/$sample.clean.cov.pdf");
  my @labels = ("fig:$sample.k$k.raw.cov",
                "fig:$sample.k$k.clean.cov");
  my @captions = ("Raw coverage",
                  "Clean coverage");

  my $label = "fig:$sample.k$k.cov";
  my $caption = "Sample `$sample' (k=$k) coverage";

  my $thresh;
  if(defined($thresh = $kmer_cleaning{$sample}->{$k})) {
    $caption .= " (cleaned off \$<$thresh\$)";
  }

  mksubfloat(\@images, \@labels, \@captions, $label, $caption, 0.9);
}

sub mk_len_fig
{
  my ($k,$sample) = @_;
  my @images = ("$dir/k$k/graphs/$sample.raw.len.pdf",
                "$dir/k$k/graphs/$sample.clean.len.pdf");
  my @labels = ("fig:$sample.k$k.raw.len",
                "fig:$sample.k$k.clean.len");
  my @captions = ("Raw unitig",
                  "Clean unitig");

  my $label = "fig:$sample.k$k.unitig.lens";
  my $caption = "Sample `$sample' (\$k=$k\$) unitig lengths";

  my $thresh;
  if(defined($thresh = $kmer_cleaning{$sample}->{$k})) {
    $caption .= " (cleaned off coverage \$<$thresh\$)";
  }

  mksubfloat(\@images, \@labels, \@captions, $label, $caption, 0.4);
}

sub mk_links_fig
{
  my ($k,$sample) = @_;
  my @images = ("$dir/k$k/links/$sample.se.links.pdf",
                "$dir/k$k/links/$sample.pe.links.pdf");
  my @labels = ("fig:$sample.k$k.se.links",
                "fig:$sample.k$k.pe.links");
  my @captions = ("Raw SE links",
                  "Raw PE links");

  my $label = "fig:$sample.k$k.link.cov";
  my $caption = "Sample `$sample' (\$k=$k\$) link coverage";

  my ($thresh_se, $thresh_pe);
  if(defined($thresh_se = $link_cleaning{$sample}->{$k}->{'se'})) {
    $captions[0] .= " (cleaned off \$<$thresh_se\$)";
  }
  if(defined($thresh_pe = $link_cleaning{$sample}->{$k}->{'pe'})) {
    $captions[1] .= " (cleaned off \$<$thresh_pe\$)";
  }

  mksubfloat(\@images, \@labels, \@captions, $label, $caption, 0.4);
}


print '
\documentclass[a4paper]{article}

\usepackage{mathtools}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage{url}
\usepackage{amssymb}

\newcommand{\srcpath}{'.$srcpath.'}

\title{McCortex Pipeline Report: '.$dir.'}
\author{McCortex make-report.pl}
\begin{document}

\maketitle
\tableofcontents

\section{Summary}

This file was generated by McCortex: \url{https://github.com/mcveanlab/mccortex}.
\begin{itemize}
\item kmers: $'.join(', ', @kmers).'$
\item samples: '.join(', ', @samples).'
\end{itemize}

';

for my $sample (@samples) {
  print "\\section{Sample: `$sample'}\n";
  for my $k (@kmers) {
    print "\\subsection{`$sample' k=\$$k\$}\n";
    mk_cov_fig($k, $sample);
    mk_len_fig($k, $sample);
    mk_links_fig($k, $sample);
    print '\clearpage'."\n";
  }
}

print '

\end{document}

';
