#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max min sum);
use File::Path qw(make_path);
use File::Copy;
use File::Basename; # dirname()

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin."/perl/";

use McCortexScripts; # mccortex_maxk()

# TODO:
# [ ] Get raw/clean number of kmers from log files
# [ ] Get raw/clean number of kmers by running `mccortex31 view ...`
# [ ] Generate our own link coverage plots using `mccortex31 links ...`
# [ ] Plot links of length 1, 2, 3 each in their own plot
# [ ] Plot distribution of link lengths (no coverage info)

# Get kmers:
#   $srcdir/k*
# Get samples:
#   $srcdir/k$k/graphs/*.raw.ctx
# Get kmer coverage:
#   $srcdir/k$k/graphs/$sample.clean.kmercov
# Get readlen:
#   $srcdir/k$k/graphs/$sample.{raw|clean}.ctx
# Get cleaning thresholds:
#   $srcdir/k$k/graphs/$sample.clean.ctx

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage $0 [options] <srcdir> <reportdir>

  Generate a report in latex on a McCortex pipeline run.

  Example:
    $0 mcrun mcrun_report
    cd mcrun_report
    make

";

  exit(-1);
}

if(@ARGV != 2) { print_usage(); }

my $srcdir = shift(@ARGV);
my $outdir = shift(@ARGV);

my $ctxdir = dirname(__FILE__)."/..";

#
# Gather intel
#
my @samples = ();
my @kmers = ();

# Crawl
my @files;

# find kmers
@files = ls($srcdir);
@kmers = map {$_ =~ /^k(\d+)$/} (grep {/^k\d+$/} ls($srcdir));
# find samples
@files = map {lspath("$srcdir/k$_/graphs/", "$srcdir/k$_/links/")} @kmers;
@samples = rmdups(map {$_ =~ /\/([^\/]*).(?:raw|clean).ctx$/} (grep {/graphs\/.*?.(?:raw|clean).ctx$/} @files));
my %pathhash = ();
map {$pathhash{$_} = 1} @files;

# Debug: list all paths found
# my @tmp = sort keys(%pathhash);
# print "files:".join("\n", @tmp)."\n";

print "Kmers: @kmers\n";
print "Samples: @samples\n";

# Get cleaning thresholds:
my %kmer_cleaning = (); # {$sample}->{$k} => kmer cleaning threshold
my %link_cleaning = (); # {$sample}->{$k}->{"se"|"pe"} => cleaning threshold
# Get kmer coverage:
my %kmer_cov = (); # {$sample}->{$k} => kmer cov
my %readlen = (); # {$sample}->{$k} => readlen

for my $sample (@samples) {
  $kmer_cleaning{$sample} = {};
  $link_cleaning{$sample} = {};
  $kmer_cov{$sample} = {};
  $readlen{$sample} = {};
  for my $k (@kmers) {
    $link_cleaning{$sample}->{$k} = {};
  }
}

for my $sample (@samples) {
  for my $k (@kmers) {
    my $path;
    my $maxk = mccortex_maxk($k);
    my $mccortex = dirname(__FILE__)."/../bin/mccortex$maxk";

    if(!(-e $mccortex)) {
      die("Please compile with `make MAXK=$maxk` (missing: $mccortex)");
    }

    # read length
    my $path1 = "$srcdir/k$k/graphs/$sample.raw.ctx";
    my $path2 = "$srcdir/k$k/graphs/$sample.clean.ctx";
    $path = undef;
    if(defined($pathhash{$path1})) { $path = $path1; }
    if(defined($pathhash{$path2})) { $path = $path2; }
    if(defined($path)) {
      my $cmd = "$mccortex view -q --info $path";
      my $info = `$cmd` or die("Command failed: $! ($cmd)");
      my ($rlen) = ($info =~ /mean input contig length:\s*([0-9\.]+)/);
      if(!defined($rlen)) { die("Can't get readlen length: $path"); }
      $readlen{$sample}->{$k} = int($rlen+0.5);
    }

    # Kmer cleaning threshold
    $path = "$srcdir/k$k/graphs/$sample.clean.ctx.log";
    if(defined($pathhash{$path})) {
      my $thresh;
      open(FH, $path) or die("Cannot read $path: $!");
      while(defined(my $line = <FH>)) {
        if($line =~ /Removing unitigs with coverage < (\d+)/i) { $thresh = $1; last; }
      }
      close(FH);
      if(defined($thresh)) { $kmer_cleaning{$sample}->{$k} = $thresh; }
    }

    # link cleaning threshold
    for my $type (qw(se pe)) {
      $path = "$srcdir/k$k/links/$sample.$type.thresh.txt";
      if(defined($pathhash{$path})) {
        my $thresh;
        open(FH, $path) or die("Cannot read $path: $!");
        while(defined(my $line = <FH>)) {
          if($line =~ /suggested_cutoff=(\d+)/i) { $thresh = $1; last; }
        }
        close(FH);
        if(defined($thresh)) { $link_cleaning{$sample}->{$k}->{$type} = $thresh; }
      }
    }

    # kmer coverage
    $path = "$srcdir/k$k/graphs/$sample.clean.kmercov";
    if(defined($pathhash{$path})) {
      open(FH, $path) or die("Cannot read $path: $!");
      my $line = <FH>;
      if(defined($line)) {
        chomp($line);
        if($line =~ /^(\d+)$/) {
          $kmer_cov{$sample}->{$k} = $1;
        } else {
          print STDERR "Cannot process: $line\n";
        }
      }
      close(FH);
    }
  }
}



print "Creating report directory...\n";

#
# Make output directory
#
create_dir($outdir);
create_dir("$outdir/plots");
create_dir("$outdir/data");
create_dir("$outdir/scripts");


sub create_dir
{
  if(!(-d $_[0])) {
    make_path($_[0]) or die("Cannot create dir: $_[0]");
  }
}

#
# Copy data to:
#
# graph:
#   data/<sample>.k<k>.raw.cov.csv
#   data/<sample>.k<k>.clean.cov.csv
#   data/<sample>.k<k>.raw.len.csv
#   data/<sample>.k<k>.clean.len.csv
#   data/<sample>.k<k>.kmercov
#   data/<sample>.k<k>.kthresh
# links:
#   data/<sample>.k<k>.se.links.csv
#   data/<sample>.k<k>.pe.links.csv
#   data/<sample>.k<k>.se.links.thresh
#   data/<sample>.k<k>.pe.links.thresh
#

# List of plots we are able to create
my %plots = ();

sub attempt_copy
{
  my ($src, $dst, $plot) = @_;
  if(defined($pathhash{$src})) {
    copy($src, $dst) or die("Copy failed: $src -> $dst");
    if(defined($plot)) { $plots{$plot} = 1; }
  }
}

sub write_to_file
{
  my ($path, $text) = @_;
  open(FH, ">$path") or die("Cannot write to $path");
  print FH $text;
  close(FH);
}

sub copy_script
{
  my ($src, $dst) = @_;
  copy($src, $dst) or die("Copy failed: $src -> $dst");
  chmod(0777, $dst) or die("chmod failed: $dst [0777]");
}

# Copy scripts
copy_script("$ctxdir/scripts/report/make-kmer-plot.sh",
            "$outdir/scripts/make-kmer-plot.sh");

copy_script("$ctxdir/scripts/report/make-link-plot.sh",
            "$outdir/scripts/make-link-plot.sh");

copy_script("$ctxdir/scripts/R/install-deps.R",
            "$outdir/scripts/install-deps.R");

copy_script("$ctxdir/scripts/R/plot-covg-hist.R",
            "$outdir/scripts/plot-covg-hist.R");

copy_script("$ctxdir/scripts/R/plot-length-hist.R",
            "$outdir/scripts/plot-length-hist.R");

copy_script("$ctxdir/scripts/R/link-cov-heatmap.R",
            "$outdir/scripts/link-cov-heatmap.R");

# Copy data
for my $sample (@samples) {
  for my $k (@kmers) {
    attempt_copy("$srcdir/k$k/graphs/$sample.raw.cov.csv",
                 "$outdir/data/$sample.k$k.raw.cov.csv",
                 "plots/$sample.k$k.raw.cov.pdf");
    attempt_copy("$srcdir/k$k/graphs/$sample.clean.cov.csv",
                 "$outdir/data/$sample.k$k.clean.cov.csv",
                 "plots/$sample.k$k.clean.cov.pdf");
    attempt_copy("$srcdir/k$k/graphs/$sample.raw.len.csv",
                 "$outdir/data/$sample.k$k.raw.len.csv",
                 "plots/$sample.k$k.raw.len.pdf");
    attempt_copy("$srcdir/k$k/graphs/$sample.clean.len.csv",
                 "$outdir/data/$sample.k$k.clean.len.csv",
                 "plots/$sample.k$k.clean.len.pdf");
    if(defined(my $t = $kmer_cov{$sample}->{$k})) {
      write_to_file("$outdir/data/$sample.k$k.kmercov", "$t\n");
    }
    if(defined(my $t = $readlen{$sample}->{$k})) {
      write_to_file("$outdir/data/$sample.k$k.readlen", "$t\n");
    }
    if(defined(my $t = $kmer_cleaning{$sample}->{$k})) {
      write_to_file("$outdir/data/$sample.k$k.kthresh", "$t\n");
    }
    attempt_copy("$srcdir/k$k/links/$sample.se.links.csv",
                 "$outdir/data/$sample.k$k.se.links.csv",
                 "plots/$sample.k$k.se.links.pdf");
    attempt_copy("$srcdir/k$k/links/$sample.pe.links.csv",
                 "$outdir/data/$sample.k$k.pe.links.csv",
                 "plots/$sample.k$k.pe.links.pdf");
    for my $type (qw(se pe)) {
      if(defined(my $t = $link_cleaning{$sample}->{$k}->{$type})) {
        write_to_file("$outdir/data/$sample.k$k.$type.links.thresh", "$t\n");
      }
    }
  }
}

# Debug: print all plots we have
# print STDERR "plots:".join("\n", sort keys(%plots))."\n";

#
# Write Makefile
#
open(FH, ">$outdir/Makefile") or die("Cannot write Makefile");
print FH "
PDFLATEX=pdflatex
KCOV_PLOTTER=Rscript scripts/plot-covg-hist.R
UNITIG_PLOTTER=Rscript scripts/plot-length-hist.R
LINK_PLOTTER=Rscript scripts/link-cov-heatmap.R
MKLINKPLOT=scripts/make-link-plot.sh
MKKMERPLOT=scripts/make-kmer-plot.sh

CSV_FILES=\$(wildcard data/*.csv)
PLOTS=\$(patsubst data/%.csv,plots/%.pdf,\$(CSV_FILES))

all: report.pdf

# Print any variable with `make -f file.mk print-VARNAME`
print-%:
  \@echo '\$*=\$(\$*)'

# will also use data/%.kthresh if available
# will also use data/%.kmercov if available
plots/%.cov.pdf: data/%.cov.csv
\t\$(MKKMERPLOT) \"\$(KCOV_PLOTTER)\" \$< \$@

plots/%.len.pdf: data/%.len.csv
\t\$(UNITIG_PLOTTER) \$< \$@

# will also use data/%.{se|pe}.links.thresh if available
# will also use data/%.kmercov if available
plots/%.links.pdf: data/%.links.csv
\t\$(MKLINKPLOT) \"\$(LINK_PLOTTER)\" \$< \$@

report.pdf: report.tex \$(PLOTS)
\t\$(PDFLATEX) \$< \$@
\t\$(PDFLATEX) \$< \$@

clean:
\trm -rf report.pdf plots/*.pdf
";
close(FH);

#
# Write report.tex
#
my $fh;
open($fh, ">$outdir/report.tex") or die("Cannot write report.tex");
print $fh '
\documentclass[a4paper]{article}

\usepackage{mathtools}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage{url}
\usepackage{amssymb}

\title{McCortex Pipeline Report: '.$srcdir.'}
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
  print $fh "\\section{Sample: `$sample'}\n";
  for my $k (@kmers) {
    print $fh "\\subsection{`$sample' k=\$$k\$}\n";
    if(defined(my $cov = $kmer_cov{$sample}->{$k})) {
      print $fh "Kmer coverage = $cov\n";
    }
    mk_cov_fig($fh, $sample, $k);
    mk_len_fig($fh, $sample, $k);
    mk_links_fig($fh, $sample, $k);
    print $fh '\clearpage'."\n";
  }
  print $fh "\n";
}
print $fh '
\end{document}
';
close($fh);

print "Done.\n";


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
  my ($fh, $images, $labels, $captions, $label, $caption, $width) = @_;
  print $fh '
\begin{figure}[hp]
\centering
  \subfloat['.$captions->[0].']{
    ';
  if(defined($plots{$images->[0]})) {
    print $fh '\includegraphics[width='.$width.'\textwidth]{'.latexpath($images->[0]).'}';
  } else { print $fh '\hspace{'.$width.'\textwidth}'; } print $fh '
    \label{'.$labels->[0].'}
  }
  \hfill
  \subfloat['.$captions->[1].']{
    ';
  if(defined($plots{$images->[1]})) {
    print $fh '\includegraphics[width='.$width.'\textwidth]{'.latexpath($images->[1]).'}';
  } else { print $fh '\hspace{'.$width.'\textwidth}'; } print $fh '
    \label{'.$labels->[1].'}
  }
  \caption{'.$caption.'}
  \label{'.$label.'}
\end{figure}
';
}

sub mk_cov_fig
{
  my ($fh, $sample, $k) = @_;
  my @images = ("plots/$sample.k$k.raw.cov.pdf",
                "plots/$sample.k$k.clean.cov.pdf");
  my @labels = ("fig:$sample.k$k.raw.cov",
                "fig:$sample.k$k.clean.cov");
  my @captions = ("Raw coverage",
                  "Clean coverage");

  my $label = "fig:$sample.k$k.cov";
  my $caption = "Sample `$sample' (k=$k) coverage.";

  my ($thresh, $cov);
  if(defined($thresh = $kmer_cleaning{$sample}->{$k})) {
    $caption .= " Cleaned \$<$thresh\$ (solid line).";
  }
  if(defined($cov = $kmer_cov{$sample}->{$k})) {
    $caption .= " Mean kmer coverage \$=$cov\$ (dashed line).";
  }

  mksubfloat($fh, \@images, \@labels, \@captions, $label, $caption, 0.9);
}

sub mk_len_fig
{
  my ($fh, $sample, $k) = @_;
  my @images = ("plots/$sample.k$k.raw.len.pdf",
                "plots/$sample.k$k.clean.len.pdf");
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

  mksubfloat($fh, \@images, \@labels, \@captions, $label, $caption, 0.4);
}

sub mk_links_fig
{
  my ($fh, $sample, $k) = @_;
  my @sources = ("plots/$sample.k$k.se.links.csv",
                 "plots/$sample.k$k.pe.links.csv");
  my @images = ("plots/$sample.k$k.se.links.pdf",
                "plots/$sample.k$k.pe.links.pdf");
  my @labels = ("fig:$sample.k$k.se.links",
                "fig:$sample.k$k.pe.links");
  my @captions = ("Raw SE links",
                  "Raw PE links");

  my $label = "fig:$sample.k$k.link.cov";
  my $caption = "Sample `$sample' (\$k=$k\$) link coverage.";

  my $rlen;
  if(defined($rlen = $readlen{$sample}->{$k})) {
    $caption .= " Mean read length: $rlen.";
  }

  my ($thresh_se, $thresh_pe);
  if(defined($thresh_se = $link_cleaning{$sample}->{$k}->{'se'})) {
    $captions[0] .= " (cleaned off \$<$thresh_se\$)";
  }
  if(defined($thresh_pe = $link_cleaning{$sample}->{$k}->{'pe'})) {
    $captions[1] .= " (cleaned off \$<$thresh_pe\$)";
  }

  if(!defined($pathhash{$sources[0]})) {
    $captions[0] .= " (not generated).";
  }
  if(!defined($pathhash{$sources[1]})) {
    $captions[1] .= " (not generated).";
  }

  mksubfloat($fh, \@images, \@labels, \@captions, $label, $caption, 0.4);
}

