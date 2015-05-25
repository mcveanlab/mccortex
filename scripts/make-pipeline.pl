#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use POSIX qw/strftime/;
use CortexScripts;

#
# TODO:
# [ ] Merge info fields when merging VCF files (waiting on bcftools request)
# [ ] Add pooled cleaning (for low coverage samples)
# [ ] 1-by-1 bubble/breakpoint calling for lower memory
# [ ] genotyping sites
# [ ] pass genome size / fetch from ref FASTA
# [x] add option to use stampy to map
# [x] thread: take fragment length min/max length
# [x] bubbles: take max allele + flank lengths
# [x] calls2vcf: take min mapq value
# [x] pop bubbles for diploid when assembling contigs
# [x] option to call variants without using links
# [x] unitigs target to dump unitigs
#

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./make-pipeline.pl [options] <list-of-kmers> <out-dir> <samples.txt>
  Generate a Makefile to run common McCortex pipelines

  Options:
    -r,--ref <ref.fa>             Reference sequence
    -s,--stampy <path/stampy.py>  Use stampy instead of BWA to place variants
    -S,--stampy-base <B>          Stampy hashes <B>.stidx and <B>.sthash

  Example:
    ./make-pipeline.pl 31:39:2 my_proj samples.txt > job.mk
    make -f job.mk bubbles-vcf

  To list all the commands without running:
    make -f job.mk --always-make --dry-run bubbles-vcf

  <kmers> specifies which kmers are to be used. It must be a comma separated
  list e.g. '21,33', or of the form <firstK[:lastK[:stepK]]>. Examples:
    '27' => 27;  '27:31' => 27,29,31; '27:39:4' => 27,31,35,49

  <samples.txt> should space or tab separated with 2-4 columns of the format:
    # comment
    <sample_name> <se_file,...> <pefile1:pefile2,...> <interleaved_file,...>
    ...
";
  exit(-1);
}

my $cmd = "$0 @ARGV";

my $default_mem = "1G";
my $default_ctxdir = "~/mccortex/";
my $default_nthreads = 2;

my $ref_path; # path to reference FASTA if available

my $stampy;
my $stampy_base; # base to stampy hash files (.stidx, .sthash)

# Parse command line args
while(@ARGV > 3) {
  my $arg = shift;
  if($arg =~ /^(-r|--ref)$/ && !defined($ref_path)) { $ref_path = shift; }
  elsif($arg =~ /^(-s|--stampy)$/ && !defined($stampy)) { $stampy = shift; }
  elsif($arg =~ /^(-S|--stampy-base)$/ && !defined($stampy_base)) { $stampy_base = shift; }
  else { print_usage("Unknown argument: $arg"); }
}

if(@ARGV != 3) { print_usage(); }

# Check if stampy is used, set it up
# Otherwise we use BWA instead
if(defined($stampy) && !defined($ref_path)) { die("Gave --stampy <S> without --ref <R>"); }
if(defined($stampy_base) && !defined($stampy)) { die("Gave --stampy-base <B> without --stampy <S>"); }

if(defined($ref_path) && !defined($stampy_base)) {
  if($ref_path =~ /^(.*?)(.fa|.fa.gz|.fq|.fq.gz)?$/) {
    $stampy_base = $1;
  } else {
    $stampy_base = $ref_path;
  }
  if($stampy_base eq "") { die("Please pass --stampy-base <B>"); }
}

my @kmers = parse_kmer_list($ARGV[0]);
my $proj = $ARGV[1];
my $sample_path = $ARGV[2];

print STDERR "kmers: ".join(', ', @kmers)."\n";
print STDERR "outdir: $proj\n";
print STDERR "sample_file: $sample_path\n";

# Load samples file
# returns array of ({'name','se_files','pe_files','i_files'}, ...)
my @samples = load_samples_file($sample_path);
if(@samples == 0) { die("No samples given in: $sample_path"); }
print STDERR "sample_names: ".join(', ', map {$_->{'name'}} @samples)."\n";

my $union_bubble_links_vcf = "$proj/vcfs/bubbles.links.".join('.',map {"k$_"} @kmers).".vcf.gz";
my $union_brkpnt_links_vcf = "$proj/vcfs/breakpoints.links.".join('.',map {"k$_"} @kmers).".vcf.gz";
my $union_bubble_plain_vcf = "$proj/vcfs/bubbles.plain.".join('.',map {"k$_"} @kmers).".vcf.gz";
my $union_brkpnt_plain_vcf = "$proj/vcfs/breakpoints.plain.".join('.',map {"k$_"} @kmers).".vcf.gz";

print '# '.strftime("%F %T", localtime($^T)).'
#
# Generated with:
#     '.$cmd.'
#
# To use this file:
#     make -f <thisfile> [options] [target]
#
# Valid targets:
#   graphs          <- build and clean graphs
#   links           <- build and clean links
#   bubbles         <- make bubble calls
#   breakpoints     <- make breakpoint calls
#   bubbles-vcf     <- make bubble vcf
#   breakpoints-vcf <- make breakpoint vcf
#   contigs         <- assemble contigs for each sample
#   contigs-pop     <- assemble contigs after popping bubbles
#   unitigs         <- dump unitigs for each sample
#
# Make targets without using links:
#   plain-bubbles         <- make bubble calls
#   plain-breakpoints     <- make breakpoint calls
#   plain-bubbles-vcf     <- make bubble vcf
#   plain-breakpoints-vcf <- make breakpoint vcf
#
# Options:
#   --dry-run              Print commands but not run them
#   --always-make          List all commands even if dependencies exist.
#   CTXDIR=<mccortex-dir>  e.g. CTXDIR=~/bin/mccortex
#   MEM=<mem-to-use>       e.g. MEM=80G
#   NTHREADS=<nthreads>
#

SHELL=/bin/bash -eou pipefail

# General options
CTXDIR='.$default_ctxdir.'
MEM='.$default_mem.'
NTHREADS='.$default_nthreads.'
# Reference sequence (FASTA/FASTQ file) leave blank if none
REF_FILE='.(defined($ref_path) ? $ref_path : '').'
# Matepair orientation of library (FR,FF,RR,RF)
MATEPAIR=FR

# Command arguments
BUILD_ARGS=--fq-cutoff 10 --cut-hp 10 --keep-pcr --matepair $(MATEPAIR)
KMER_CLEANING_ARGS=--fallback 2
POP_BUBBLES_ARGS=--max-diff 50 --max-covg 5
THREAD_ARGS=--min-frag-len 150 --max-frag-len 1000 --fq-cutoff 5 --matepair $(MATEPAIR) --one-way --gap-diff-const 5 --gap-diff-coeff 0.1
LINK_CLEANING_ARGS=--limit 5000 --threshold 0.001
BREAKPOINTS_ARGS=--minref 20
BUBBLES_ARGS=--max-allele 3000 --max-flank 1000
CALL2VCF_ARGS=--max-align 500 --max-allele 100 --min-mapq 30
CONTIG_ARGS=--no-missing-check --confid-step 0.99
CONTIG_POP_ARGS=--confid-step 0.99

# Paths to scripts
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
VCFSORT=$(CTXDIR)/libs/biogrok/vcf-sort
VCFRENAME=$(CTXDIR)/libs/biogrok/vcf-rename

# Third party libraries packaged in McCortex
BWA=$(CTXDIR)/libs/bwa/bwa
BGZIP=$(CTXDIR)/libs/htslib/bgzip
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
STAMPY='.(defined($stampy) ? $stampy : '').'
STAMPY_BASE='.(defined($stampy_base) ? $stampy_base : '').'

# Set up memory, threads and number of kmers in the graph
CTX_ARGS=
ifdef MEM
  CTX_ARGS:=$(CTX_ARGS) -m $(MEM)
endif
ifdef NKMERS
  CTX_ARGS:=$(CTX_ARGS) -n $(NKMERS)
endif
ifdef NTHREADS
  CTX_ARGS:=$(CTX_ARGS) -t $(NTHREADS)
endif

';

for my $k (@kmers) {
  print "# Files at k=$k\n";
  print "RAW_GRAPHS_K$k=".join(' ', map {"$proj/k$k/graphs/$_->{'name'}.raw.ctx"} @samples)."\n";
  print "CLEAN_GRAPHS_K$k=\$(RAW_GRAPHS_K$k:.raw.ctx=.clean.ctx)\n";
  print "CLEAN_UNITIGS_K$k=\$(CLEAN_GRAPHS_K$k:.ctx=.unitigs.fa.gz)\n";
  print "RAW_SE_LINKS_K$k=". join(' ', map {"$proj/k$k/links/$_->{'name'}.se.raw.ctp.gz"} @samples)."\n";
  print "RAW_PE_LINKS_K$k=". join(' ', map {"$proj/k$k/links/$_->{'name'}.pe.raw.ctp.gz"} @samples)."\n";
  print "CLEAN_SE_LINKS_K$k=\$(RAW_SE_LINKS_K$k:.raw.ctp.gz=.clean.ctp.gz)\n";
  print "CLEAN_PE_LINKS_K$k=\$(RAW_PE_LINKS_K$k:.raw.ctp.gz=.clean.ctp.gz)\n";
  print "BUBBLES_LINKS_K$k=$proj/k$k/bubbles/bubbles.txt.gz\n";
  print "BUBBLES_PLAIN_K$k=$proj/k$k/bubbles_plain/bubbles.txt.gz\n";
  print "CONTIGS_K$k=".join(' ', map {"$proj/k$k/contigs/$_->{'name'}.rmdup.fa.gz"} @samples)."\n";
  print "CONTIGS_POP_K$k=".join(' ', map {"$proj/k$k/contigs/$_->{'name'}.pop.rmdup.fa.gz"} @samples)."\n";
  if(defined($ref_path)) {
    print "BREAKPOINTS_LINKS_K$k=$proj/k$k/breakpoints/breakpoints.txt.gz\n";
    print "BREAKPOINTS_PLAIN_K$k=$proj/k$k/breakpoints_plain/breakpoints.txt.gz\n";
  } else {
    print "BREAKPOINTS_LINKS_K$k=\n";
    print "BREAKPOINTS_PLAIN_K$k=\n";
  }
  print "\n";
}

print "RAW_GRAPHS=" .join(' ', map {"\$(RAW_GRAPHS_K$_)"}  @kmers)."\n";
print "CLEAN_GRAPHS=\$(RAW_GRAPHS:.raw.ctx=.clean.ctx)\n";
print "CLEAN_UNITIGS=\$(CLEAN_GRAPHS:.ctx=.unitigs.fa.gz)\n";
print "RAW_LINKS="  .join(' ', map {"\$(RAW_SE_LINKS_K$_) \$(RAW_PE_LINKS_K$_) "} @kmers)."\n";
print "CLEAN_SE_LINKS=".join(' ', map {"\$(CLEAN_SE_LINKS_K$_)"} @kmers)."\n";
print "CLEAN_PE_LINKS=".join(' ', map {"\$(CLEAN_PE_LINKS_K$_)"} @kmers)."\n";
print "FINAL_LINKS=\$(CLEAN_PE_LINKS)\n";
print "BUBBLES_LINKS="    .join(' ', map {"\$(BUBBLES_LINKS_K$_)"}        @kmers)."\n";
print "BUBBLES_PLAIN="    .join(' ', map {"\$(BUBBLES_PLAIN_K$_)"}        @kmers)."\n";
print "BREAKPOINTS_LINKS=".join(' ', map {"\$(BREAKPOINTS_LINKS_K$_)"}    @kmers)."\n";
print "BREAKPOINTS_PLAIN=".join(' ', map {"\$(BREAKPOINTS_PLAIN_K$_)"}    @kmers)."\n";
print "CONTIGS="    .join(' ', map {"\$(CONTIGS_K$_)"}        @kmers)."\n";
print "CONTIGS_POP=".join(' ', map {"\$(CONTIGS_POP_K$_)"}    @kmers)."\n";

my @dirlist = ();
for my $k (@kmers) {
  my $dirs = join(' ', "$proj/k$k/graphs/", "$proj/k$k/links/",
                       "$proj/k$k/contigs/",
                       "$proj/k$k/bubbles/", "$proj/k$k/breakpoints/",
                       "$proj/k$k/bubbles_plain/", "$proj/k$k/breakpoints_plain/",
                       "$proj/k$k/ref/");
  push(@dirlist, $dirs);
}
push(@dirlist, "$proj/vcfs/");

print 'DIRS='.join(" \\\n     ", @dirlist).'

COVG_CSV_FILES=$(RAW_GRAPHS:.raw.ctx=.raw.covg.csv)

# .csi are index files (for VCF in this case)
BUBBLES_LINKS_VCFS=$(BUBBLES_LINKS:.txt.gz=.norm.vcf.gz)
BUBBLES_LINKS_CSIS=$(BUBBLES_LINKS_VCFS:=.csi)
BREAKPOINTS_LINKS_VCFS=$(BREAKPOINTS_LINKS:.txt.gz=.norm.vcf.gz)
BREAKPOINTS_LINKS_CSIS=$(BREAKPOINTS_LINKS_VCFS:=.csi)

BUBBLES_PLAIN_VCFS=$(BUBBLES_PLAIN:.txt.gz=.norm.vcf.gz)
BUBBLES_PLAIN_CSIS=$(BUBBLES_PLAIN_VCFS:=.csi)
BREAKPOINTS_PLAIN_VCFS=$(BREAKPOINTS_PLAIN:.txt.gz=.norm.vcf.gz)
BREAKPOINTS_PLAIN_CSIS=$(BREAKPOINTS_PLAIN_VCFS:=.csi)

CALL_FILES=$(BUBBLES_LINKS) $(BUBBLES_PLAIN) $(BREAKPOINTS_LINKS) $(BREAKPOINTS_PLAIN)
RAW_VCFS=$(CALL_FILES:.txt.gz=.raw.vcf)

# CALL_VCFS=$(CALL_FILES:.txt.gz=.norm.vcf.gz)
# CALL_CSIS=$(BUBBLES_LINKS_CSIS) $(BREAKPOINTS_LINKS_CSIS)
# VCF_TMP_FILES=$(BUBBLES:.txt.gz=.flanks.fa.gz) $(BUBBLES:.txt.gz=.flanks.sam) \
#               $(CALL_FILES:.txt.gz=.sort.vcf) $(CALL_FILES:.txt.gz=.norm.vcf)

# Referece Graphs
';

if(defined($ref_path)) {
  for my $k (@kmers) { print "REF_GRAPH_K$k=$proj/k$k/ref/ref.ctx\n"; }
} else {
  for my $k (@kmers) { print "REF_GRAPH_K$k=\n"; }
}

print '# REF_GRAPHS='.join(' ', map {'$(REF_GRAPH_K'.$_.')'} @kmers).'

# HAVE_LOGS=$(RAW_GRAPHS) $(CLEAN_GRAPHS) $(CLEAN_UNITIGS) $(REF_GRAPHS) $(RAW_LINKS) $(CLEAN_SE_LINKS) $(CLEAN_PE_LINKS) $(LINK_TMP_FILES) $(CALL_FILES) $(RAW_VCFS)
# LOG_FILES=$(HAVE_LOGS:=.log)

# Mark all dependencies as secondary
# It means don\'t re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:

# Delete files if their recipe fails
.DELETE_ON_ERROR:

# Remove in-built rules for certain file suffixes
.SUFFIXES:

all: ' .(defined($ref_path) ? 'bubbles-vcf breakpoints-vcf' : 'bubbles').' unitigs | checks

graphs: $(CLEAN_GRAPHS) | checks

unitigs: $(CLEAN_UNITIGS) | checks

links: $(FINAL_LINKS) | checks

bubbles: links-bubbles
bubbles-vcf: links-bubbles-vcf
breakpoints: links-breakpoints
breakpoints-vcf: links-breakpoints-vcf

links-bubbles: $(BUBBLES_LINKS) | checks
plain-bubbles: $(BUBBLES_PLAIN) | checks

contigs: $(CONTIGS) | checks
contigs-pop: $(CONTIGS_POP) | checks

checks:'."\n";
my @ctx_maxks = get_maxk_values(@kmers);
for my $maxk (@ctx_maxks) {
  print "\t@[ -x \$(CTXDIR)/bin/mccortex$maxk ] || ( echo 'Error: Please compile cortex with `make MAXK=$maxk` or pass CTXDIR=<path/to/mccortex/>' 1>&2 && false )\n";
}
print "\n";

# Can only create VCFs if we have a reference
if(defined($ref_path)) {
  print "links-breakpoints: \$(BREAKPOINTS_LINKS)\n\n";
  print "plain-breakpoints: \$(BREAKPOINTS_PLAIN)\n\n";
  print "links-bubbles-vcf: $union_bubble_links_vcf $union_bubble_links_vcf.csi\n\n";
  print "links-breakpoints-vcf: $union_brkpnt_links_vcf $union_brkpnt_links_vcf.csi\n\n";
  print "plain-bubbles-vcf: $union_bubble_plain_vcf $union_bubble_plain_vcf.csi\n\n";
  print "plain-breakpoints-vcf: $union_brkpnt_plain_vcf $union_brkpnt_plain_vcf.csi\n\n";
} else {
  my @tgts = qw(links-breakpoints     plain-breakpoints
                links-breakpoints-vcf plain-breakpoints-vcf
                links-bubbles-vcf     plain-bubbles-vcf);
  for my $tgt (@tgts) {
    print "$tgt:\n\t\@echo 'Need to give make-pipeline.pl --ref <r.fa> to run $tgt 2>1 && false\n\n";
  }
}

print "\$(DIRS):
\tmkdir -p \$@

clean:
\t\@echo To delete: rm -rf $proj

.PHONY: all clean checks graphs links unitigs contigs contigs-pop
.PHONY: bubbles breakpoints bubbles-vcf breakpoints-vcf
.PHONY: links-bubbles links-breakpoints links-bubbles-vcf links-breakpoints-vcf
.PHONY: plain-bubbles plain-breakpoints plain-bubbles-vcf plain-breakpoints-vcf

";

# Create and clean graph files
print "#\n# Build graph files\n#\n";
for my $k (@kmers) {
  my $ctx = get_mccortex($k);

  # Build reference
  if(defined($ref_path)) {
    print "# reference at k=$k\n";
    print "$proj/k$k/ref/ref.ctx: $ref_path | \$(DIRS)\n";
    print "\t$ctx build \$(CTX_ARGS) \$(BUILD_ARGS) -k $k --sample ref --seq \$< \$@ >& \$@.log\n\n";
  }

  print "# building sample graphs at k=$k\n";
  for my $sample (@samples) {
    # Create raw graph file
    my $sname = $sample->{'name'};
    my @files = get_all_sample_files($sample);

    print "$proj/k$k/graphs/$sname.raw.ctx: ".join(' ', @files)." | \$(DIRS)\n";
    print "\t$ctx build \$(CTX_ARGS) -k $k --sample $sname " .
          join(' ', (map {"--seq $_"}               @{$sample->{'se_files'}}),
                    (map {"--seq2 $_->[0]:$_->[1]"} @{$sample->{'pe_files'}}),
                    (map {"--seqi $_"}              @{$sample->{'i_files'}})) .
          ' $@ >& $@.log'."\n\n";
  }

  # Pop bubbles
  print "# Generate individual graphs for sample assembly with high covg indiv.\n";
  print "# Clean and pop bubbles at k=$k\n";
  print "$proj/k$k/graphs/%.pop.raw.covg.csv: $proj/k$k/graphs/%.pop.clean.ctx\n";
  print "$proj/k$k/graphs/%.pop.clean.ctx: $proj/k$k/graphs/%.raw.ctx\n";
  print "\t$ctx clean \$(CTX_ARGS) \$(KMER_CLEANING_ARGS) --covg-before $proj/k$k/graphs/\$*.pop.raw.covg.csv -o \$@ \$< >& \$@.log\n";
  print "$proj/k$k/graphs/%.pop.clean.ctx: $proj/k$k/graphs/%.pop.clean.ctx\n";
  print "\t$ctx popbubbles \$(CTX_ARGS) \$(POP_BUBBLES_ARGS) -o \$@ \$< >& \$@.log\n\n";

  # Clean graph files at k=$k
  print "# sample graph cleaning at k=$k\n";
  print "$proj/k$k/graphs/%.raw.covg.csv: $proj/k$k/graphs/%.clean.ctx\n";
  print "$proj/k$k/graphs/%.clean.ctx: $proj/k$k/graphs/%.raw.ctx\n";
  print "\t$ctx clean \$(CTX_ARGS) \$(KMER_CLEANING_ARGS) --covg-before $proj/k$k/graphs/\$*.raw.covg.csv -o \$@ \$< >& \$@.log\n";
  print "\t$ctx inferedges \$(CTX_ARGS) \$@ >& $proj/k$k/graphs/\$*.inferedges.ctx.log\n\n";

  # Dump unitigs
  print "# sample graph unitigs at k=$k\n";
  print "$proj/k$k/graphs/%.clean.unitigs.fa.gz: $proj/k$k/graphs/%.clean.ctx\n";
  print "\t($ctx unitigs \$(CTX_ARGS) \$< | gzip -c > \$@) 2> \$@.log\n\n";
}

# Create and clean link files
print "#\n# Generate link files\n#\n";
for my $k (@kmers) {
  print "# creating links at k=$k\n";
  my $ctx = get_mccortex($k);

  my @samples_with_pop = @samples;
  for my $s (@samples) {
    my %tmp = %$s;
    $tmp{'name'} .= '.pop';
    push(@samples_with_pop, \%tmp);
  }

  for my $sample (@samples_with_pop) {
    my $sname = $sample->{'name'};
    my @pe_files = (map {($_->[0], $_->[1])} @{$sample->{'pe_files'}},
                    @{$sample->{'i_files'}});
    my @se_files = (@{$sample->{'se_files'}}, @pe_files);

    my $ctx_clean_file    = "$proj/k$k/graphs/$sname.clean.ctx";
    my $ctp_se_raw_file   = "$proj/k$k/links/$sname.se.raw.ctp.gz";
    my $ctp_pe_raw_file   = "$proj/k$k/links/$sname.pe.raw.ctp.gz";
    my $ctp_se_clean_file = "$proj/k$k/links/$sname.se.clean.ctp.gz";
    my $ctp_pe_clean_file = "$proj/k$k/links/$sname.pe.clean.ctp.gz";

    # 1. Single ended threading: submit all reads as single ended first
    print "$ctp_se_raw_file: $ctx_clean_file @se_files | \$(DIRS)\n";
    print "\t$ctx thread \$(CTX_ARGS) \$(THREAD_ARGS) " .
          join(' ', (map {"--seq $_"}                    @{$sample->{'se_files'}}),
                    (map {"--seq $_->[0] --seq $_->[1]"} @{$sample->{'pe_files'}}),
                    (map {"--seq $_"}                    @{$sample->{'i_files'}})) .
          ' -o $@ $< >& $@.log'."\n\n";

    # 2. If we have any paired end reads, add that information in a second pass
    if(@{$sample->{'pe_files'}} > 0 || @{$sample->{'i_files'}} > 0) {
      print "$ctp_pe_raw_file: $ctx_clean_file $ctp_se_clean_file @pe_files | \$(DIRS)\n";
      print "\t$ctx thread \$(CTX_ARGS) \$(THREAD_ARGS) -p $ctp_se_clean_file " .
              join(' ', (map {"--seq2 $_->[0]:$_->[1]"} @{$sample->{'pe_files'}}),
                        (map {"--seqi $_"}              @{$sample->{'i_files'}})) .
              ' -o $@ $< >& $@.log'."\n\n";
    } else {
      print "$ctp_pe_clean_file: $ctp_se_clean_file\n";
      print "\tln \$< \$@\n\n";
    }
  }

  # Clean link files at k=$k
  my $ctp_raw_file     = "$proj/k$k/links/%.raw.ctp.gz";
  my $ctp_clean_file   = "$proj/k$k/links/%.clean.ctp.gz";
  my $ctp_thresh_file  = "$proj/k$k/links/%.thresh.txt";

  # Generate coverage CSV from first N kmers with links
  print "# link cleaning at k=$k\n";
  print "$ctp_thresh_file: $ctp_raw_file\n";
  print "\t$ctx links \$(LINK_CLEANING_ARGS) \$< > \$@ 2> \$@.log\n\n";

  print "$ctp_clean_file: $ctp_raw_file $ctp_thresh_file\n";
  print "\tTHRESH=`grep 'suggested_cutoff=' $proj/k$k/links/\$*.thresh.txt | grep -oE '[0-9,]+\$\$'`; \\\n";
  print "\t$ctx links -c \"\$\$THRESH\" -o \$@ \$< >& \$@.log\n\n";
}

# Assemble contigs
print "#\n# Assemble contigs\n#\n";
for my $k (@kmers) {
  my $ctx = get_mccortex($k);
  print "# assembly high covg sample k=$k\n";
  print "$proj/k$k/contigs/%.pop.raw.fa.gz: $proj/k$k/graphs/%.pop.clean.ctx $proj/k$k/links/%.pop.pe.clean.ctp.gz\n";
  print "\t( $ctx contigs \$(CTX_ARGS) \$(CONTIG_POP_ARGS) -o - -p $proj/k$k/links/\$*.pop.pe.clean.ctp.gz \$<               | gzip -c > \$@ ) >& \$@.log\n\n";
  print "# assembly k=$k\n";
  print "$proj/k$k/contigs/%.raw.fa.gz: $proj/k$k/graphs/%.clean.ctx $proj/k$k/links/%.pe.clean.ctp.gz \$(REF_GRAPH_K$k)\n";
  print "\t( $ctx contigs \$(CTX_ARGS) \$(CONTIG_ARGS) -o - -p $proj/k$k/links/\$*.pe.clean.ctp.gz \$< \$(REF_GRAPH_K$k) | gzip -c > \$@ ) >& \$@.log\n\n";
  print "# Remove redundant contigs\n";
  print "$proj/k$k/contigs/%.rmdup.fa.gz: $proj/k$k/contigs/%.raw.fa.gz\n";
  print "\t( $ctx rmsubstr -m \$(MEM) -k $k -o - \$< | gzip -c > \$@ ) >& \$@.log\n\n";
}

# Generate buble calls
print "#\n# Make bubble calls\n#\n";
for my $k (@kmers) {
  my $ctx = get_mccortex($k);
  my $link_args = get_p_args($k);
  my $hapcol = defined($ref_path) ? "--haploid ".scalar(@samples) : '';

  print "# bubble calls k=$k\n";
  print "$proj/k$k/bubbles/bubbles.txt.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) $hapcol -o \$@ $link_args \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";

  print "# bubble calls k=$k (without links)\n";
  print "$proj/k$k/bubbles_plain/bubbles.txt.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) $hapcol -o \$@            \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";
}

# Some things require a reference to be used
if(defined($ref_path))
{
  # Generate breakpoint calls
  print "#\n# Make breakpoint calls\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);
    my $link_args = get_p_args($k);

    print "# breakpoint calls k=$k\n";
    print "$proj/k$k/breakpoints/breakpoints.txt.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@ $link_args \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";

    print "# breakpoint calls k=$k\n";
    print "$proj/k$k/breakpoints_plain/breakpoints.txt.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@            \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";
  }

  # Generate buble VCFs
  print "#\n# Make bubble raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);

    print "# bubbles raw VCF k=$k\n";
    print "$proj/k$k/%/bubbles.flanks.fa.gz: $proj/k$k/%/bubbles.txt.gz\n";
    print "\t\$(CTXFLANKS) \$< > \$@\n\n";

    print "$proj/k$k/%/bubbles.raw.vcf: $proj/k$k/%/bubbles.txt.gz $proj/k$k/%/bubbles.flanks.sam \$(REF_FILE)\n";
    print "\t$ctx calls2vcf \$(CALL2VCF_ARGS) -F $proj/k$k/\$*/bubbles.flanks.sam -o \$@ \$< \$(REF_FILE) >& \$@.log\n\n";
  }

  # Generate breakpoint VCFs
  print "#\n# Make breakpoint raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);
    my $breakpoint_file    = "";
    my $breakpoint_raw_vcf = "";

    print "# breakpoints raw VCF k=$k\n";
    print "$proj/k$k/%/breakpoints.raw.vcf: $proj/k$k/%/breakpoints.txt.gz \$(REF_FILE)\n";
    print "\t$ctx calls2vcf \$(CALL2VCF_ARGS) -o \$@ \$< \$(REF_FILE) >& \$@.log\n\n";
  }

  # Post-processing rules for VCFs
  print "#\n# Post-processing for raw VCFs\n#\n";
  print "$proj/%.sort.vcf: $proj/%.raw.vcf\n";
  print "\t\$(VCFSORT) \$< > \$@\n\n";

  print "$proj/%.norm.vcf.gz: $proj/%.sort.vcf \$(REF_FILE)\n";
  print "\t\$(BCFTOOLS) norm --site-win 5000 --remove-duplicates --fasta-ref \$(REF_FILE) --multiallelics +both \$< | \\\n";
  print "\t\$(VCFRENAME) > $proj/\$*.norm.vcf\n";
  print "\t\$(BGZIP) -f $proj/\$*.norm.vcf\n\n";

  # Generate union VCF
  print "#\n# Create union compressed VCF\n#\n";
  print "$union_bubble_links_vcf: \$(BUBBLES_LINKS_VCFS) \$(BUBBLES_LINKS_CSIS)\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BUBBLES_LINKS_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_bubble_plain_vcf: \$(BUBBLES_PLAIN_VCFS) \$(BUBBLES_PLAIN_CSIS)\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BUBBLES_PLAIN_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_links_vcf: \$(BREAKPOINTS_LINKS_VCFS) \$(BREAKPOINTS_LINKS_CSIS)\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BREAKPOINTS_LINKS_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_plain_vcf: \$(BREAKPOINTS_PLAIN_VCFS) \$(BREAKPOINTS_PLAIN_CSIS)\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BREAKPOINTS_PLAIN_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";


  print "#\n# General VCF rules\n#\n";
  # Compress a VCF
  # This deletes the existing file which may upset Make
#  print "%.vcf.gz: %.vcf\n";
#  print "\t\$(BGZIP) -f \$<\n\n";

  # Create VCF index files .vcf.gz.csi
  print "%.vcf.gz.csi: %.vcf.gz\n";
  print "\t\$(BCFTOOLS) index -f \$<\n\n";

  my $ms = defined($stampy) ? "" : "# ";
  my $mb = defined($stampy) ? "# " : "";


  # BWA
  print "# Mapping with BWA\n";
  print $mb."\$(REF_FILE).bwt: \$(REF_FILE)\n";
  print $mb."\t\$(BWA) index \$(REF_FILE)\n\n";

  print $mb."%.sam: %.fa.gz \$(REF_FILE).bwt \$(REF_FILE)\n";
  print $mb."\t\$(BWA) mem \$(REF_FILE) \$< > \$@\n\n";

  # Stampy
  print "# Mapping with Stampy\n";
  print $ms."\$(STAMPY_BASE).stidx: \$(REF_FILE)\n";
  print $ms."\t\$(STAMPY) -G \$(STAMPY_BASE) \$(REF_FILE)\n\n";

  print $ms."\$(STAMPY_BASE).sthash: \$(STAMPY_BASE).stidx \$(REF_FILE)\n";
  print $ms."\t\$(STAMPY) -g \$(STAMPY_BASE) -H \$(STAMPY_BASE)\n\n";

  print $ms."%.sam: %.fa.gz \$(STAMPY_BASE).stidx \$(STAMPY_BASE).sthash\n";
  print $ms."\t\$(STAMPY) -g \$(STAMPY_BASE) -h \$(STAMPY_BASE) -M \$< > \$@\n\n";
}


print STDERR "Usage: make -f <script> [options] [target]\n";
print STDERR "  --always-run          Run/list all commands, inc. those already run\n";
print STDERR "  --dry-run             List commands, don't run them\n";
print STDERR "  CTXDIR=<mccortexdir>  Path to McCortex directory e.g. CTXDIR=~/mccortex\n";
print STDERR "  MEM=<MEM>             Maximum memory to use e.g. MEM=80G\n";
print STDERR "  NTHREADS=<N>          Maximum number of job threads to use\n";
print STDERR "\n";

# Done!
exit(0);


# Load sample file
# ## Comment lines
# <sample>  <se_file,...>  <pe_file1:file2,...>  <interleaved_files,...>
# returns array of ({'name','se_files','pe_files','i_files'}, ...)
sub load_samples_file
{
  my ($path) = @_;
  my @samples = ();
  my %sample_names = ();
  my $sfh = open_file($path);
  while(defined(my $line = <$sfh>)) {
    if($line !~ /^\s*$/ && $line !~ /^#/) {
      my @cols = split(/\s/, $line);
      if(@cols < 2 || @cols > 4) { die("Bad line"); }
      my ($sname, $se_txt, $pe_txt, $i_txt) = @cols;
      # Check sample name is sane and unique
      if($sname !~ /^[a-z0-9_\-\.]+$/i) { die("Bad name: $sname"); }
      if($sname =~ /\.pop$/) { die("sample name cannot end '.pop'"); }
      if(defined($sample_names{$sname})) { die("Duplicate sample name"); }
      # Parse file lists
      my @se_files = parse_file_list($se_txt);
      my @pe_files = parse_pe_file_list($pe_txt);
      my @i_files  = parse_file_list($i_txt);
      push(@samples, {'name'     => $sname,
                      'se_files' => \@se_files,
                      'pe_files' => \@pe_files,
                      'i_files'  => \@i_files});
    }
  }
  close($sfh);
  return @samples;
}

sub get_p_args
{
  my ($k) = @_;
  return join(' ', map {"-p $_:$proj/k$k/links/$samples[$_]->{'name'}.pe.clean.ctp.gz"} 0..$#samples);
}

sub get_all_sample_files
{
  my ($sample) = @_;
  return (@{$sample->{'se_files'}},
          (map {($_->[0], $_->[1])} @{$sample->{'pe_files'}}),
          @{$sample->{'i_files'}});
}

sub get_required_binaries
{
  return map { get_mccortex($_) } get_maxk_values(@_);
}

sub get_maxk_values
{
  my %maxks = ();
  for my $k (@_) { $maxks{(int(($k+31)/32) * 32 - 1)} = 1; }
  return keys %maxks;
}

sub get_mccortex
{
  my ($k) = @_;
  return "\$(CTXDIR)/bin/mccortex".(int(($k+31)/32) * 32 - 1);
}

# Split a comma separated, colon delimited list of PE files
# "A.1.fa:A.2.fa,B.1.fa:B.2.fa"
#  => (["A.1.fa","A.2.fa"],["B.1.fa","B.2.fa"])
sub parse_pe_file_list
{
  my ($txt) = @_;
  my @files = parse_file_list($txt);
  my @pe_files = ();
  for my $f (@files) {
    if($f =~ /^([^:]+):([^:]+)$/) { push(@pe_files, [$1, $2]); }
    else { die("Bad PE line: $txt"); }
  }
  return @pe_files;
}

# Split a list of files by commas. Check there are not empty entries
sub parse_file_list
{
  my ($txt) = @_;
  if(!defined($txt) || $txt eq "" || $txt eq "." || $txt eq "-") { return (); }
  my @files = split(',', $txt);
  for my $f (@files) { if($f eq "") { die("Empty file entry: $txt"); } }
  return @files;
}

# '31'      => (31)
# '31:41'   => (31,33,35,37,39,41)
# '31:41:4' => (31,35,39)
sub parse_kmer_list
{
  my ($txt) = @_;
  if($txt =~ /^(\d+(,\d+)*)$/) {
    my @ks = split(',',$txt);
    for my $k (@ks) { if($k % 2 == 0) { die("Kmers must not be odd: $txt"); }}
    return @ks;
  }
  elsif($txt =~ /^(\d+)(?::(\d+)(?::(\d+))?)?$/)
  {
    my ($start,$end,$step) = ($1,$1,2);
    if(defined($2)) { $end  = $2; }
    if(defined($3)) { $step = $3; }
    if($step  % 2 != 0) { die("Kmer step must be even: $txt"); }
    if($start % 2 == 0) { die("Kmer start must be odd: $txt"); }
    if($end   % 2 == 0) { die("Kmer end must be odd: $txt"); }
    my @ks = ();
    for(my $k=$start; $k <= $end; $k += $step) { push(@ks, $k); }
    return @ks;
  }
  else { die("Poorly formatted kmer list: $txt"); }
}
