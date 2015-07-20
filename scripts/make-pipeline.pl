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
# [ ] genotyping sites
# [ ] pass genome size / fetch from ref FASTA
# [x] 1-by-1 bubble/breakpoint calling for lower memory
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
"Usage: $0 [options] <list-of-kmers> <out-dir> <samples.txt>
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

# k=(29,31), kmerstr = "k29.k31"
my $kmerstr = join('.', map {"k$_"} @kmers);
my $union_bubble_joint_links_vcf = "$proj/vcfs/bubbles.joint.links.".$kmerstr.".vcf.gz";
my $union_bubble_1by1_links_vcf  = "$proj/vcfs/bubbles.1by1.links.".$kmerstr.".vcf.gz";
my $union_brkpnt_joint_links_vcf = "$proj/vcfs/breakpoints.joint.links.".$kmerstr.".vcf.gz";
my $union_brkpnt_1by1_links_vcf  = "$proj/vcfs/breakpoints.1by1.links.".$kmerstr.".vcf.gz";
my $union_bubble_joint_plain_vcf = "$proj/vcfs/bubbles.joint.plain.".$kmerstr.".vcf.gz";
my $union_bubble_1by1_plain_vcf  = "$proj/vcfs/bubbles.1by1.plain.".$kmerstr.".vcf.gz";
my $union_brkpnt_joint_plain_vcf = "$proj/vcfs/breakpoints.joint.plain.".$kmerstr.".vcf.gz";
my $union_brkpnt_1by1_plain_vcf  = "$proj/vcfs/breakpoints.1by1.plain.".$kmerstr.".vcf.gz";

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
#   vcfs            <- make bubble+breakpoint vcfs
#   contigs         <- assemble contigs for each sample
#   contigs-pop     <- assemble contigs after popping bubbles
#   unitigs         <- dump unitigs for each sample
#
# Options:
#   --always-make          List all commands even if dependencies exist.
#   --dry-run              Print commands but not run them
#   CTXDIR=<mccortex-dir>  e.g. CTXDIR=~/bin/mccortex
#   MEM=<mem-to-use>       e.g. MEM=80G
#   NTHREADS=<nthreads>
#   USE_LINKS=<B>          <B> is "yes" or "no"
#   JOINT_CALLING=<B>      <B> is "yes" or "no"
#

#
# File structure:
# ---------------
#
# <K> is kmer size
# <S> is sample name
# <P> set of sample names and "joint"
#
# Roughly listed in order of generation
#
# <outdir>/
#   -> k<K>/
#     -> graphs/
#       -> <S>.raw.ctx
#       -> <S>.raw.ctx.log
#       -> <S>.raw.covg.csv
#       -> <S>.clean.ctx
#       -> <S>.clean.ctx.log
#       -> <S>.inferedges.ctx.log
#       -> <S>.clean.unitigs.fa.gz
#     -> links/
#       -> <S>.se.raw.ctp.gz
#       -> <S>.se.raw.ctp.gz.log
#       -> <S>.se.clean.ctp.gz
#       -> <S>.se.clean.ctp.gz.log
#       -> <S>.se.thresh.txt
#       -> <S>.se.thresh.txt.log
#       -> <S>.pe.raw.ctp.gz
#       -> <S>.pe.raw.ctp.gz.log
#       -> <S>.pe.clean.ctp.gz
#       -> <S>.pe.clean.ctp.gz.log
#       -> <S>.pe.thresh.txt
#       -> <S>.pe.thresh.txt.log
#     -> bubbles/
#       -> <P>.bub.gz
#       -> <P>.bub.gz.log
#       -> <P>.flanks.fa.gz
#       -> <P>.flanks.sam
#       -> <P>.bub.raw.vcf
#       -> <P>.bub.raw.vcf.log
#       -> <P>.bub.sort.vcf
#       -> <P>.bub.norm.vcf.gz
#       -> <P>.bub.norm.vcf.gz.csi
#     -> bubbles_plain/
#       -> SAME AS ../bubbles/
#     -> breakpoints/
#       -> <P>.brk.gz
#       -> <P>.brk.gz.log
#       -> <P>.brk.raw.vcf
#       -> <P>.brk.raw.vcf.log
#       -> <P>.brk.sort.vcf
#       -> <P>.brk.norm.vcf.gz
#       -> <P>.brk.norm.vcf.gz.csi
#     -> breakpoints_plain/
#       -> SAME AS ../breakpoints/
#     -> contigs/
#       -> <S>.raw.fa.gz
#       -> <S>.raw.fa.gz.log
#       -> <S>.rmdup.fa.gz
#       -> <S>.rmdup.fa.gz.log
#     -> ref/
#       -> ref.ctx
#       -> ref.ctx.log
#   -> vcfs/
#     -> 1by1_samples/
#       -> <S>.k29.k31.brk.norm.vcf.gz
#       -> <S>.k29.k31.brk.norm.vcf.gz.csi
#     -> <breakpoints|bubbles>.<joint|1by1>.<plain|links>.k<K>.vcf.gz
#     -> <breakpoints|bubbles>.<joint|1by1>.<plain|links>.k<K>.vcf.gz.csi
#     e.g.
#     -> breakpoints.joint.plain.k29.k31.vcf.gz
#     -> breakpoints.joint.plain.k29.k31.vcf.gz.csi
#     -> breakpoints.1by1.plain.k29.k31.vcf.gz
#     -> breakpoints.1by1.plain.k29.k31.vcf.gz.csi
#     -> breakpoints.joint.links.k29.k31.vcf.gz
#     -> breakpoints.joint.links.k29.k31.vcf.gz.csi
#     -> breakpoints.1by1.links.k29.k31.vcf.gz
#     -> breakpoints.1by1.links.k29.k31.vcf.gz.csi
#     -> bubbles.joint.plain.k29.k31.vcf.gz
#     -> bubbles.joint.plain.k29.k31.vcf.gz.csi
#     -> bubbles.1by1.plain.k29.k31.vcf.gz
#     -> bubbles.1by1.plain.k29.k31.vcf.gz.csi
#     -> bubbles.joint.links.k29.k31.vcf.gz
#     -> bubbles.joint.links.k29.k31.vcf.gz.csi
#     -> bubbles.1by1.links.k29.k31.vcf.gz
#     -> bubbles.1by1.links.k29.k31.vcf.gz.csi
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
MIN_FRAG_LEN=150
MAX_FRAG_LEN=1000

# Command arguments
BUILD_ARGS=--fq-cutoff 10 --cut-hp 10 --keep-pcr --matepair $(MATEPAIR)
KMER_CLEANING_ARGS=--fallback 2
POP_BUBBLES_ARGS=--max-diff 50 --max-covg 5
THREAD_ARGS=--min-frag-len $(MIN_FRAG_LEN) --max-frag-len $(MAX_FRAG_LEN) --fq-cutoff 5 --matepair $(MATEPAIR) --one-way --gap-diff-const 5 --gap-diff-coeff 0.1
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

#
# Parse USE_LINKS and JOINT_CALLING Makefile options
#

# Use links is default on
ifeq ($(USE_LINKS),yes)
  LINKS=1
endif
ifeq ($(USE_LINKS),true)
  LINKS=1
endif
ifeq ($(USE_LINKS),1)
  LINKS=1
endif
ifndef USE_LINKS
  LINKS=1
endif

# Joint calling is default on
ifeq ($(JOINT_CALLING),yes)
  JOINT=1
endif
ifeq ($(JOINT_CALLING),true)
  JOINT=1
endif
ifeq ($(JOINT_CALLING),1)
  JOINT=1
endif
ifndef JOINT_CALLING
  JOINT=1
endif

# LINKS is defined iff we are using links
# JOINT is defined iff we are doing joint calling

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
  print "CONTIGS_K$k=".join(' ', map {"$proj/k$k/contigs/$_->{'name'}.rmdup.fa.gz"} @samples)."\n";
  print "CONTIGS_POP_K$k=".join(' ', map {"$proj/k$k/contigs/$_->{'name'}.pop.rmdup.fa.gz"} @samples)."\n";
  print "\n";
}

sub refonly { return defined($ref_path) ? $_[0] : ""; }

print "
ifdef LINKS
\tifdef JOINT\n";
    # links+joint calling
    for my $k (@kmers) {
      print "\t\tBUBBLES_K$k=$proj/k$k/bubbles/joint.bub.gz\n";
      print "\t\tBREAKPOINTS_K$k=".refonly("$proj/k$k/breakpoints/joint.brk.gz") . "\n";
    }
    print "\t\tBUBBLES_UNION_VCFS=".refonly("$union_bubble_joint_links_vcf $union_bubble_joint_links_vcf.csi") . "\n";
    print "\t\tBREAKPOINTS_UNION_VCFS=".refonly("$union_brkpnt_joint_links_vcf $union_brkpnt_joint_links_vcf.csi") . "\n";
print "\telse\n";
    # links+1by1 calling
    for my $k (@kmers) {
      print "\t\tBUBBLES_K$k=".join(' ', map {"$proj/k$k/bubbles/$_->{'name'}.bub.gz"} @samples)."\n";
      print "\t\tBREAKPOINTS_K$k=".join(' ', map {"$proj/k$k/breakpoints/$_->{'name'}.brk.gz"} @samples)."\n";
    }
    print "\t\tBUBBLES_UNION_VCFS=".refonly("$union_bubble_1by1_links_vcf $union_bubble_1by1_links_vcf.csi") . "\n";
    print "\t\tBREAKPOINTS_UNION_VCFS=".refonly("$union_brkpnt_1by1_links_vcf $union_brkpnt_1by1_links_vcf.csi") . "\n";
print "\tendif
else
\tifdef JOINT\n";
    # plain+joint calling
    for my $k (@kmers) {
      print "\t\tBUBBLES_K$k=$proj/k$k/bubbles_plain/joint.bub.gz\n";
      print "\t\tBREAKPOINTS_K$k=".refonly("$proj/k$k/breakpoints_plain/joint.brk.gz") . "\n";
    }
    print "\t\tBUBBLES_UNION_VCFS=".refonly("$union_bubble_joint_plain_vcf $union_bubble_joint_plain_vcf.csi") . "\n";
    print "\t\tBREAKPOINTS_UNION_VCFS=".refonly("$union_brkpnt_joint_plain_vcf $union_brkpnt_joint_plain_vcf.csi") . "\n";
print "\telse\n";
    # plain+1by1 calling
    for my $k (@kmers) {
      print "\t\tBUBBLES_K$k=".join(' ', map {"$proj/k$k/bubbles_plain/$_->{'name'}.bub.gz"} @samples)."\n";
      print "\t\tBREAKPOINTS_K$k=".join(' ', map {"$proj/k$k/breakpoints_plain/$_->{'name'}.brk.gz"} @samples)."\n";
    }
    print "\t\tBUBBLES_UNION_VCFS=".refonly("$union_bubble_1by1_plain_vcf $union_bubble_1by1_plain_vcf.csi") . "\n";
    print "\t\tBREAKPOINTS_UNION_VCFS=".refonly("$union_brkpnt_1by1_plain_vcf $union_brkpnt_1by1_plain_vcf.csi") . "\n";
print "\tendif
endif\n\n";

print "RAW_GRAPHS=" .join(' ', map {"\$(RAW_GRAPHS_K$_)"}  @kmers)."\n";
print "CLEAN_GRAPHS=\$(RAW_GRAPHS:.raw.ctx=.clean.ctx)\n";
print "CLEAN_UNITIGS=\$(CLEAN_GRAPHS:.ctx=.unitigs.fa.gz)\n";
print "RAW_LINKS="  .join(' ', map {"\$(RAW_SE_LINKS_K$_) \$(RAW_PE_LINKS_K$_) "} @kmers)."\n";
print "CLEAN_SE_LINKS=".join(' ', map {"\$(CLEAN_SE_LINKS_K$_)"} @kmers)."\n";
print "CLEAN_PE_LINKS=".join(' ', map {"\$(CLEAN_PE_LINKS_K$_)"} @kmers)."\n";
print "FINAL_LINKS=\$(CLEAN_PE_LINKS)\n";
print "BUBBLES="    .join(' ', map {"\$(BUBBLES_K$_)"}        @kmers)."\n";
print "BREAKPOINTS=".join(' ', map {"\$(BREAKPOINTS_K$_)"}    @kmers)."\n";
print "CONTIGS="    .join(' ', map {"\$(CONTIGS_K$_)"}        @kmers)."\n";
print "CONTIGS_POP=".join(' ', map {"\$(CONTIGS_POP_K$_)"}    @kmers)."\n";

print "\n# Files to merge to create various union VCFs\n";
print "# .csi are index files (for VCF in this case)\n";

sub merge_vcf_list
{
  my ($isbubble,$isjoint,$use_links) = @_;
  my $dir;
  my @r = ();
  if($isbubble) { $dir = ($use_links ? "bubbles"     : "bubbles_plain"); }
  else          { $dir = ($use_links ? "breakpoints" : "breakpoints_plain"); }
  my $ext = ($isbubble ? "bub" : "brk");
  for my $k (@kmers) {
    if($isjoint) {
      push(@r, "$proj/k$k/$dir/joint.$ext.norm.vcf.gz");
    } else {
      for my $s (@samples) {
        push(@r, "$proj/k$k/$dir/$s->{'name'}.$ext.norm.vcf.gz");
      }
    }
  }
  return @r;
}

print "BUBBLES_JOINT_PLAIN_VCFS=".    join(' ', merge_vcf_list(1,1,0))."\n";
print "BUBBLES_JOINT_LINKS_VCFS=".    join(' ', merge_vcf_list(1,1,1))."\n";
print "BREAKPOINTS_JOINT_PLAIN_VCFS=".join(' ', merge_vcf_list(0,1,0))."\n";
print "BREAKPOINTS_JOINT_LINKS_VCFS=".join(' ', merge_vcf_list(0,1,1))."\n";

print "BUBBLES_1BY1_PLAIN_VCFS=".     join(' ', merge_vcf_list(1,0,0))."\n";
print "BUBBLES_1BY1_LINKS_VCFS=".     join(' ', merge_vcf_list(1,0,1))."\n";
print "BREAKPOINTS_1BY1_PLAIN_VCFS=". join(' ', merge_vcf_list(0,0,0))."\n";
print "BREAKPOINTS_1BY1_LINKS_VCFS=". join(' ', merge_vcf_list(0,0,1))."\n";

print "BUBBLES_JOINT_LINKS_CSIS=\$(BUBBLES_JOINT_LINKS_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BUBBLES_JOINT_PLAIN_CSIS=\$(BUBBLES_JOINT_PLAIN_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BREAKPOINTS_JOINT_LINKS_CSIS=\$(BREAKPOINTS_JOINT_LINKS_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BREAKPOINTS_JOINT_PLAIN_CSIS=\$(BREAKPOINTS_JOINT_PLAIN_VCFS:.vcf.gz=.vcf.gz.csi)\n";

print "BUBBLES_1BY1_LINKS_CSIS=\$(BUBBLES_1BY1_LINKS_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BUBBLES_1BY1_PLAIN_CSIS=\$(BUBBLES_1BY1_PLAIN_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BREAKPOINTS_1BY1_LINKS_CSIS=\$(BREAKPOINTS_1BY1_LINKS_VCFS:.vcf.gz=.vcf.gz.csi)\n";
print "BREAKPOINTS_1BY1_PLAIN_CSIS=\$(BREAKPOINTS_1BY1_PLAIN_VCFS:.vcf.gz=.vcf.gz.csi)\n";

print "\n";

my @dirlist = ();
for my $k (@kmers) {
  my $dirs = join(' ', "$proj/k$k/graphs/", "$proj/k$k/links/",
                       "$proj/k$k/contigs/",
                       "$proj/k$k/bubbles/", "$proj/k$k/breakpoints/",
                       "$proj/k$k/bubbles_plain/", "$proj/k$k/breakpoints_plain/",
                       "$proj/k$k/ref/");
  push(@dirlist, $dirs);
}
push(@dirlist, "$proj/vcfs/1by1_samples/");

print 'DIRS='.join(" \\\n     ", @dirlist).'

COVG_CSV_FILES=$(RAW_GRAPHS:.raw.ctx=.raw.covg.csv)

# Referece Graphs
';

for my $k (@kmers) {
  print "REF_GRAPH_K$k=".refonly("$proj/k$k/ref/ref.ctx")."\n";
}

print '
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

bubbles: $(BUBBLES) | checks
';

if(defined($ref_path)) {
  print 'bubbles-vcf: $(BUBBLES_UNION_VCFS) | checks
breakpoints: $(BREAKPOINTS) | checks
breakpoints-vcf: $(BREAKPOINTS_UNION_VCFS) | checks
';
} else {
  for my $tgt (qw(bubbles-vcf breakpoints breakpoints-vcf)) {
    print "$tgt:\n\t\@echo 'Need to give make-pipeline.pl --ref <r.fa> to run $tgt 2>1 && false\n\n";
  }
}

print 'vcfs: bubbles-vcf breakpoints-vcf

contigs: $(CONTIGS) | checks
contigs-pop: $(CONTIGS_POP) | checks

checks:'."\n";
my @ctx_maxks = get_maxk_values(@kmers);
for my $maxk (@ctx_maxks) {
  print "\t@[ -x \$(CTXDIR)/bin/mccortex$maxk ] || ( echo 'Error: Please compile cortex with `make MAXK=$maxk` or pass CTXDIR=<path/to/mccortex/>' 1>&2 && false )\n";
}

print "
\$(DIRS):
\tmkdir -p \$@

clean:
\t\@echo To delete: rm -rf $proj

.PHONY: all clean checks graphs links unitigs contigs contigs-pop
.PHONY: bubbles breakpoints bubbles-vcf breakpoints-vcf vcfs

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
      print "\t$ctx thread \$(CTX_ARGS) \$(THREAD_ARGS) -p $ctp_se_clean_file --zero-paths " .
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

  # joint bubble calling
  print "# bubble calls k=$k joint+links\n";
  print "$proj/k$k/bubbles/joint.bub.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) $hapcol -o \$@ $link_args \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";

  print "# bubble calls k=$k joint+nolinks\n";
  print "$proj/k$k/bubbles_plain/joint.bub.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) | \$(DIRS)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) $hapcol -o \$@ \$< \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";

  # 1by1 bubble calling
  print "# bubble calls k=$k 1by1+links\n";
  print "$proj/k$k/bubbles/%.bub.gz: $proj/k$k/graphs/%.clean.ctx \$(REF_GRAPH_K$k) $proj/k$k/links/%.pe.clean.ctp.gz\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) --haploid 1 -o \$@ -p $proj/k$k/links/\$*.pe.clean.ctp.gz \$< \$(REF_GRAPH_K$k) >& \$@.log\n\n";

  print "# bubble calls k=$k 1by1+nolinks\n";
  print "$proj/k$k/bubbles_plain/%.bub.gz: $proj/k$k/graphs/%.clean.ctx \$(REF_GRAPH_K$k)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) \$(BUBBLES_ARGS) --haploid 1 -o \$@ \$< \$(REF_GRAPH_K$k) >& \$@.log\n\n";
}

# Some things require a reference to be used
if(defined($ref_path))
{
  # Generate breakpoint calls
  print "#\n# Make breakpoint calls\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);
    my $link_args = get_p_args($k);

    # joint breakpoint calling
    print "# breakpoint calls k=$k joint+links\n";
    print "$proj/k$k/breakpoints/joint.brk.gz: \$(CLEAN_GRAPHS_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@ $link_args \$(CLEAN_GRAPHS_K$k) >& \$@.log\n\n";

    print "# breakpoint calls k=$k joint+nolinks\n";
    print "$proj/k$k/breakpoints_plain/joint.brk.gz: \$(CLEAN_GRAPHS_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@            \$(CLEAN_GRAPHS_K$k) >& \$@.log\n\n";

    # 1by1 breakpoint calling
    print "# breakpoint calls k=$k 1by1+links\n";
    print "$proj/k$k/breakpoints/%.brk.gz: $proj/k$k/graphs/%.clean.ctx $proj/k$k/links/%.pe.clean.ctp.gz\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@ -p $proj/k$k/links/\$*.pe.clean.ctp.gz \$< >& \$@.log\n\n";

    print "# breakpoint calls k=$k 1by1+nolinks\n";
    print "$proj/k$k/breakpoints_plain/%.brk.gz: $proj/k$k/graphs/%.clean.ctx\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) \$(BREAKPOINTS_ARGS) -s \$(REF_FILE) -o \$@ \$< >& \$@.log\n\n";
  }

  # Generate buble VCFs
  print "#\n# Make bubble raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);

    print "# bubbles raw VCF k=$k\n";
    print "$proj/k$k/%.flanks.fa.gz: $proj/k$k/%.bub.gz\n";
    print "\t\$(CTXFLANKS) \$< | gzip -c > \$@\n\n";

    print "$proj/k$k/%.bub.raw.vcf: $proj/k$k/%.bub.gz $proj/k$k/%.flanks.sam \$(REF_FILE)\n";
    print "\t$ctx calls2vcf \$(CALL2VCF_ARGS) -F $proj/k$k/\$*.flanks.sam -o \$@ \$< \$(REF_FILE) >& \$@.log\n\n";
  }

  # Generate breakpoint VCFs
  print "#\n# Make breakpoint raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_mccortex($k);
    my $breakpoint_file    = "";
    my $breakpoint_raw_vcf = "";

    print "# breakpoints raw VCF k=$k\n";
    print "$proj/k$k/%.brk.raw.vcf: $proj/k$k/%.brk.gz \$(REF_FILE)\n";
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

  print "VCF_CONCAT=\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates\n";
  print "VCF_MERGE=\$(BCFTOOLS) merge\n\n";

  my @brkpnt_1by1_links_vcfs;
  my @brkpnt_1by1_plain_vcfs;

  # need to concat all kmers for each sample before merging all samples
  if(@kmers > 1) {
    print "# Merge 1-by-1 breakpoint calls for each sample first\n";
    print "# Bubble calls are sites only, not sample so do not need merging\n";
    for my $s (map {$_->{'name'}} @samples) {
      for(my $i = 0; $i < 2; $i++) {
        my $dir = ($i == 0 ? "breakpoints" : "breakpoints_plain");
        my $type = ($i == 0 ? "links" : "plain");
        my @files_vcfs = map {"$proj/k$_/$dir/$s.brk.norm.vcf.gz"} @kmers;
        my @files_csis = map {$_.".csi"} @files_vcfs;
        print "$proj/vcfs/1by1_samples/$s.$type.$kmerstr.brk.norm.vcf.gz: @files_vcfs @files_csis | \$(DIRS)\n";
        print "\t\$(VCF_CONCAT) --output-type z --output \$@ @files_vcfs\n\n";
      }
    }

    @brkpnt_1by1_links_vcfs = map {"$proj/vcfs/1by1_samples/$_->{'name'}.links.$kmerstr.brk.norm.vcf.gz"} @samples;
    @brkpnt_1by1_plain_vcfs = map {"$proj/vcfs/1by1_samples/$_->{'name'}.plain.$kmerstr.brk.norm.vcf.gz"} @samples;
  }
  else {
    # Only one kmer -- no need to concat across kmers for each sample
    my $k = $kmers[0];
    @brkpnt_1by1_links_vcfs = map {"$proj/k$k/breakpoints/$_->{'name'}.brk.norm.vcf.gz"} @samples;
    @brkpnt_1by1_plain_vcfs = map {"$proj/k$k/breakpoints_plain/$_->{'name'}.brk.norm.vcf.gz"} @samples;
  }

  my @brkpnt_1by1_links_csis = map {$_.".csi"} @brkpnt_1by1_links_vcfs;
  my @brkpnt_1by1_plain_csis = map {$_.".csi"} @brkpnt_1by1_plain_vcfs;

  # Generate union VCF
  print "#\n# Create union compressed VCF\n#\n";
  print "$union_bubble_joint_links_vcf: \$(BUBBLES_JOINT_LINKS_VCFS) \$(BUBBLES_JOINT_LINKS_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BUBBLES_JOINT_LINKS_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_bubble_joint_plain_vcf: \$(BUBBLES_JOINT_PLAIN_VCFS) \$(BUBBLES_JOINT_PLAIN_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BUBBLES_JOINT_PLAIN_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_joint_links_vcf: \$(BREAKPOINTS_JOINT_LINKS_VCFS) \$(BREAKPOINTS_JOINT_LINKS_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BREAKPOINTS_JOINT_LINKS_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_joint_plain_vcf: \$(BREAKPOINTS_JOINT_PLAIN_VCFS) \$(BREAKPOINTS_JOINT_PLAIN_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BREAKPOINTS_JOINT_PLAIN_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_bubble_1by1_links_vcf: \$(BUBBLES_1BY1_LINKS_VCFS) \$(BUBBLES_1BY1_LINKS_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BUBBLES_1BY1_LINKS_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_bubble_1by1_plain_vcf: \$(BUBBLES_1BY1_PLAIN_VCFS) \$(BUBBLES_1BY1_PLAIN_CSIS)\n";
  print "\t\$(VCF_CONCAT) \$(BUBBLES_1BY1_PLAIN_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_1by1_links_vcf: @brkpnt_1by1_links_vcfs @brkpnt_1by1_links_csis\n";
  print "\t\$(VCF_MERGE) @brkpnt_1by1_links_vcfs | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "$union_brkpnt_1by1_plain_vcf: @brkpnt_1by1_plain_vcfs @brkpnt_1by1_plain_csis\n";
  print "\t\$(VCF_MERGE) @brkpnt_1by1_plain_vcfs | \\\n";
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
print STDERR "  --always-make         Run/list all commands, inc. those already run\n";
print STDERR "  --dry-run             List commands, don't run them\n";
print STDERR "  CTXDIR=<mccortexdir>  Path to McCortex directory e.g. CTXDIR=~/mccortex\n";
print STDERR "  MEM=<MEM>             Maximum memory to use e.g. MEM=80G\n";
print STDERR "  NTHREADS=<N>          Maximum number of job threads to use\n";
print STDERR "  USE_LINKS=<B>         <B> is 'yes' or 'no'\n";
print STDERR "  JOINT_CALLING=<B>     <B> is 'yes' or 'no'\n";
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
      if($sname eq "joint") { die("sample name cannot be 'joint'"); }
      if($sname eq "1by1") { die("sample name cannot be '1by1'"); }
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
