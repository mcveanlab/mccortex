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
# * Merge info fields when merging VCF files
# * Add pooled cleaning (for low coverage samples)
# * 1-by-1 bubble/breakpoint calling for lower memory
# * genotyping
# * use genome size
# * use stampy to map
# * take paths to ref resources (ref_fa, ref_stampy, ref_bwa, ref_ctx)
#

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./make-pipeline.pl [options] <list-of-kmers> <out-dir> <samples.txt>
  Generate a Makefile to run common McCortex pipelines

  Options:
    -r,--ref <ref.fa>  Reference sequence

  Example:
    ./make-pipeline.pl 31:39:2 my_proj samples.txt > job.mk
    make -f job.mk bubblevcf

  To list all the commands without running:
    make -f job.mk --always-make --dry-run bubblevcf

  <kmers> specifies which kmers are to be used. It must be a comma separated
  list e.g. '21,33', or of the form <firstK[:lastK[:stepK]]>. Examples:
    '27' => 27;  '27:31' => 27,29,31; '27:39:4' => 27,31,35,49

  <samples.txt> should space or tab separated with 2-4 columns of the format:
    # comment
    <sample-name> <se_file,...> <pefile1:pefile2,...> <interleaved_file,...>
    ...
";
  exit(-1);
}

my $args = "$0 @ARGV";

my $default_mem = "1G";
my $default_ctxdir = "~/mccortex/";
my $default_nthreads = 2;
# Sample 5000 kmers to pick link threshold limit
my $default_link_clean_nkmers = 5000;

my $ref_path; # path to reference FASTA if available

# Parse command line args
while(@ARGV > 3) {
  my $arg = shift;
  if($arg =~ /^(-r|--ref)$/) { $ref_path = shift; }
  else { print_usage("Unknown argument: $arg"); }
}

if(@ARGV != 3) { print_usage(); }

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

my $union_bubble_vcf = "$proj/vcfs/bubbles.".join('.',map {"k$_"} @kmers).".vcf.gz";
my $union_brkpnt_vcf = "$proj/vcfs/breakpoints.".join('.',map {"k$_"} @kmers).".vcf.gz";

print '# '.strftime("%F %T", localtime($^T)).'
#
# Generated with:
#     make-pipeline.pl $args
#
# To use this file:
#     make -f <thisfile> graphs        <- build and clean graphs
#     make -f <thisfile> links         <- build and clean links
#     make -f <thisfile> bubbles       <- make bubble calls
#     make -f <thisfile> breakpoints   <- make breakpoint calls
#     make -f <thisfile> bubblevcf     <- make bubble vcf
#     make -f <thisfile> breakpointvcf <- make breakpoint vcf
#     make -f <thisfile> vcfs          <- make all vcfs including union
#     make -f <thisfile> contigs       <- assemble contigs for each sample
#     make -f <thisfile> <outdir>/k<K>/contigs/<S>.rmdup.fa.gz
#                          ^- assemble contigs for sample <S> with k=<K>
#
# Make will automatically generate dependencies.
# Add option --dry-run to print commands but not run them. Include option
# --always-make to list all commands even if dependencies exist.
#
# Other options:
#    CTXDIR=<path-to-ctx-dir>
#    MEM=<mem-to-use>
#    NTHREADS=<nthreads>
#

SHELL=/bin/bash -eou pipefail

# Override these values when running
CTXDIR='.$default_ctxdir.'
MEM='.$default_mem.'
NTHREADS='.$default_nthreads.'
LINK_CLEAN_FDR=0.001
CLEANING_ARGS=
LINK_CLEAN_NKMERS='.$default_link_clean_nkmers.'
REF_FILE='.(defined($ref_path) ? $ref_path : '').'

ifdef NKMERS
  CTX_ARGS=-m $(MEM) -t $(NTHREADS) -n $(NKMERS)
else
  CTX_ARGS=-m $(MEM) -t $(NTHREADS)
endif

# Paths to scripts
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
VCFSORT=$(CTXDIR)/scripts/bash/vcf-sort
VCFRENAME=$(CTXDIR)/scripts/bash/vcf-rename

# Third lib paths
BWA=$(CTXDIR)/libs/bwa/bwa
BGZIP=$(CTXDIR)/libs/htslib/bgzip
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools

';

for my $k (@kmers) {
  print "# Files at k=$k\n";
  print "RAW_GRAPHS_K$k=".join(' ', map {"$proj/k$k/graphs/$_->{'name'}.raw.ctx"} @samples)."\n";
  print "CLEAN_GRAPHS_K$k=\$(RAW_GRAPHS_K$k:.raw.ctx=.clean.ctx)\n";
  print "RAW_SE_LINKS_K$k=". join(' ', map {"$proj/k$k/links/$_->{'name'}.se.raw.ctp.gz"} @samples)."\n";
  print "RAW_PE_LINKS_K$k=". join(' ', map {"$proj/k$k/links/$_->{'name'}.pe.raw.ctp.gz"} @samples)."\n";
  print "CLEAN_SE_LINKS_K$k=\$(RAW_SE_LINKS_K$k:.raw.ctp.gz=.clean.ctp.gz)\n";
  print "CLEAN_PE_LINKS_K$k=\$(RAW_PE_LINKS_K$k:.raw.ctp.gz=.clean.ctp.gz)\n";
  print "BUBBLES_K$k=$proj/k$k/bubbles/bubbles.txt.gz\n";
  print "CONTIGS_K$k=".join(' ', map {"$proj/k$k/contigs/$_->{'name'}.rmdup.fa.gz"} @samples)."\n";
  if(defined($ref_path)) {
    print "BREAKPOINTS_K$k=$proj/k$k/breakpoints/breakpoints.txt.gz\n";
  } else {
    print "BREAKPOINTS_K$k=\n";
  }
  print "\n";
}

print "RAW_GRAPHS=" .join(' ', map {"\$(RAW_GRAPHS_K$_)"}  @kmers)."\n";
print "CLEAN_GRAPHS=\$(RAW_GRAPHS:.raw.ctx=.clean.ctx)\n";
print "RAW_LINKS="  .join(' ', map {"\$(RAW_SE_LINKS_K$_) \$(RAW_PE_LINKS_K$_) "} @kmers)."\n";
print "CLEAN_SE_LINKS=".join(' ', map {"\$(CLEAN_SE_LINKS_K$_)"} @kmers)."\n";
print "CLEAN_PE_LINKS=".join(' ', map {"\$(CLEAN_PE_LINKS_K$_)"} @kmers)."\n";
print "FINAL_LINKS=\$(CLEAN_PE_LINKS)\n";
print "BUBBLES="    .join(' ', map {"\$(BUBBLES_K$_)"}        @kmers)."\n";
print "BREAKPOINTS=".join(' ', map {"\$(BREAKPOINTS_K$_)"}    @kmers)."\n";
print "CONTIGS="    .join(' ', map {"\$(CONTIGS_K$_)"}        @kmers)."\n";

my @dirlist = ();
for my $k (@kmers) {
  my $dirs = join(' ', "$proj/k$k/graphs/", "$proj/k$k/links/",
                       "$proj/k$k/contigs/",
                       "$proj/k$k/bubbles/", "$proj/k$k/breakpoints/",
                       "$proj/k$k/ref/");
  push(@dirlist, $dirs);
}
push(@dirlist, "$proj/vcfs/");

print 'DIRS='.join(" \\\n     ", @dirlist).'

COVG_CSV_FILES=$(RAW_GRAPHS:.raw.ctx=.raw.covg.csv)

# .csi are index files (for VCF in this case)
BUBBLE_VCFS=$(BUBBLES:.txt.gz=.norm.vcf.gz)
BUBBLE_CSIS=$(BUBBLE_VCFS:=.csi)
BREAKPOINT_VCFS=$(BREAKPOINTS:.txt.gz=.norm.vcf.gz)
BREAKPOINT_CSIS=$(BREAKPOINT_VCFS:=.csi)
CALL_FILES=$(BUBBLES) $(BREAKPOINTS)
RAW_VCFS=$(CALL_FILES:.txt.gz=.raw.vcf)
CALL_VCFS=$(CALL_FILES:.txt.gz=.norm.vcf.gz)
CALL_CSIS=$(BUBBLE_CSIS) $(BREAKPOINT_CSIS)
VCF_TMP_FILES=$(BUBBLES:.txt.gz=.flanks.fa.gz) $(BUBBLES:.txt.gz=.flanks.sam) \
              $(CALL_FILES:.txt.gz=.sort.vcf) $(CALL_FILES:.txt.gz=.norm.vcf)

# Referece Graphs
';

if(defined($ref_path)) {
  for my $k (@kmers) { print "REF_GRAPH_K$k=$proj/k$k/ref/ref.ctx\n"; }
} else {
  for my $k (@kmers) { print "REF_GRAPH_K$k=\n"; }
}

print 'REF_GRAPHS='.join(' ', map {'$(REF_GRAPH_K'.$_.')'} @kmers).'

HAVE_LOGS=$(RAW_GRAPHS) $(CLEAN_GRAPHS) $(REF_GRAPHS) $(RAW_LINKS) $(CLEAN_SE_LINKS) $(CLEAN_PE_LINKS) $(LINK_TMP_FILES) $(CALL_FILES) $(RAW_VCFS)
LOG_FILES=$(HAVE_LOGS:=.log)

# Mark all dependencies as secondary
# It means don\'t re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:

# Delete files if their recipe fails
.DELETE_ON_ERROR:

all: ' .(defined($ref_path) ? 'bubblevcf breakpointvcf' : 'bubbles').' | checks

graphs: $(CLEAN_GRAPHS) | checks

links: $(FINAL_LINKS) | checks

bubbles: $(BUBBLES) | checks

contigs: $(CONTIGS) | checks

checks:'."\n";
my @ctx_maxks = get_maxk_values(@kmers);
for my $maxk (@ctx_maxks) {
  print "\t@[ -x \$(CTXDIR)/bin/mccortex$maxk ] || ( echo 'Error: Please compile cortex with `make MAXK=$maxk` or pass CTXDIR=<path/to/mccortex/>' 1>&2 && false )\n";
}
print "\n";

# Can only create VCFs if we have a reference
if(defined($ref_path)) {
  print "breakpoints: \$(BREAKPOINTS)\n\n";
  print "bubblevcf: $union_bubble_vcf $union_bubble_vcf.csi\n\n";
  print "breakpointvcf: $union_brkpnt_vcf $union_brkpnt_vcf.csi\n\n";
} else {
  for my $tgt (qw(breakpoints bubblevcf breakpointvcf)) {
    print "$tgt:\n\t\@echo 'Need to give make-pipeline.pl --ref <r.fa> to run $tgt 2>1 && false\n\n";
  }
}

print "\$(DIRS):
\tmkdir -p \$@

clean:
\t\@echo To delete: rm -rf $proj

.PHONY: all clean checks graphs links contigs bubbles breakpoints bubblevcf breakpointvcf

";

# Create and clean graph files
print "#\n# Build graph files\n#\n";
for my $k (@kmers) {
  my $ctx = get_ctx($k);

  # Build reference
  if(defined($ref_path)) {
    print "# reference at k=$k\n";
    print "$proj/k$k/ref/ref.ctx: $ref_path | \$(DIRS)\n";
    print "\t$ctx build \$(CTX_ARGS) -k $k --sample ref --seq \$< \$@ >& \$@.log\n\n";
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

  # Clean graph files at k=$k
  print "# sample graph cleaning at k=$k\n";
  print "$proj/k$k/graphs/%.raw.covg.csv $proj/k$k/graphs/%.clean.ctx: $proj/k$k/graphs/%.raw.ctx\n";
  print "\t($ctx clean \$(CTX_ARGS) --covg-before $proj/k$k/graphs/\$*.raw.covg.csv -o $proj/k$k/graphs/\$*.clean.ctx \$(CLEANING_ARGS) \$<; \\\n";
  print "\t $ctx inferedges \$(CTX_ARGS) $proj/k$k/graphs/\$*.clean.ctx) >& $proj/k$k/graphs/\$*.clean.ctx.log\n\n";
}

# Create and clean link files
print "#\n# Generate link files\n#\n";
for my $k (@kmers) {
  print "# creating links at k=$k\n";
  my $ctx = get_ctx($k);

  for my $sample (@samples) {
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
    print "\t$ctx thread \$(CTX_ARGS) " .
          join(' ', (map {"--seq $_"}                    @{$sample->{'se_files'}}),
                    (map {"--seq $_->[0] --seq $_->[1]"} @{$sample->{'pe_files'}}),
                    (map {"--seq $_"}                    @{$sample->{'i_files'}})) .
          ' -o $@ $< >& $@.log'."\n\n";

    # 2. If we have any paired end reads, add that information in a second pass
    if(@{$sample->{'pe_files'}} > 0 || @{$sample->{'i_files'}} > 0) {
      print "$ctp_pe_raw_file: $ctx_clean_file $ctp_se_clean_file @pe_files | \$(DIRS)\n";
      print "\t$ctx thread \$(CTX_ARGS) -p $ctp_se_clean_file " .
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
  print "\t$ctx links -L \$(LINK_CLEAN_NKMERS) -T \$(LINK_CLEAN_FDR) \$< > \$@ 2> \$@.log\n\n";

  print "$ctp_clean_file: $ctp_raw_file $ctp_thresh_file\n";
  print "\tTHRESH=`grep 'suggested_cutoff=' $proj/k$k/links/\$*.thresh.txt | grep -oE '[0-9,]+\$\$'`; \\\n";
  print "\t$ctx links -c \"\$\$THRESH\" -o \$@ \$< >& \$@.log\n\n";
}

# Assemble contigs
print "#\n# Assemble contigs\n#\n";
for my $k (@kmers) {
  my $ctx = get_ctx($k);
  print "# assembly k=$k\n";
  print "$proj/k$k/contigs/%.raw.fa.gz: $proj/k$k/graphs/%.clean.ctx $proj/k$k/links/%.pe.clean.ctp.gz \$(REF_GRAPH_K$k) | \$(DIRS)\n";
  print "\t$ctx contigs \$(CTX_ARGS) -o \$@ -p $proj/k$k/links/\$*.pe.clean.ctp.gz \$< \$(REF_GRAPH_K$k) >& \$@.log\n\n";
  print "$proj/k$k/contigs/%.rmdup.fa.gz: $proj/k$k/contigs/%.raw.fa.gz\n";
  print "\t$ctx rmsubstr -m \$(MEM) -k $k -o \$@ \$< >& \$@.log\n\n";
}

# Generate buble calls
print "#\n# Make bubble calls\n#\n";
for my $k (@kmers) {
  my $ctx = get_ctx($k);
  my $ctp_txt = get_p_args($k);
  print "# bubble calls k=$k\n";
  print "$proj/k$k/bubbles/bubbles.txt.gz: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
  print "\t$ctx bubbles \$(CTX_ARGS) -o \$@ $ctp_txt \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";
}

# Some things require a reference to be used
if(defined($ref_path))
{
  # Generate breakpoint calls
  print "#\n# Make breakpoint calls\n#\n";
  for my $k (@kmers) {
    my $ctx = get_ctx($k);
    my $ctp_txt = get_p_args($k);
    my $brkpnt_file = "$proj/k$k/breakpoints/breakpoints.txt.gz";

    print "# breakpoint calls k=$k\n";
    print "$brkpnt_file: \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) \$(CLEAN_PE_LINKS_K$k) | \$(DIRS)\n";
    print "\t$ctx breakpoints \$(CTX_ARGS) -s \$(REF_FILE) -o \$@ $ctp_txt \$(CLEAN_GRAPHS_K$k) \$(REF_GRAPH_K$k) >& \$@.log\n\n";
  }

  # Generate buble VCFs
  print "#\n# Make bubble raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_ctx($k);
    my $bubbles_file = "$proj/k$k/bubbles/bubbles.txt.gz";

    my $bubble_flanks_fa_file  = "$proj/k$k/bubbles/bubbles.flanks.fa.gz";
    my $bubble_flanks_sam_file = "$proj/k$k/bubbles/bubbles.flanks.sam";
    my $raw_bubble_vcf         = "$proj/k$k/bubbles/bubbles.raw.vcf";

    print "# bubbles raw VCF k=$k\n";
    print "$bubble_flanks_fa_file: $bubbles_file\n";
    print "\t\$(CTXFLANKS) \$< > \$@\n\n";

    print "$bubble_flanks_sam_file: $bubble_flanks_fa_file \$(REF_FILE)\n";
    print "\t\$(BWA) index \$(REF_FILE)\n";
    print "\t\$(BWA) mem \$(REF_FILE) \$< > \$@\n\n";

    print "$raw_bubble_vcf: $bubbles_file $bubble_flanks_sam_file \$(REF_FILE)\n";
    print "\t$ctx calls2vcf -F $bubble_flanks_sam_file -o \$@ \$< \$(REF_FILE) >& \$@.log\n\n";
  }

  # Generate breakpoint VCFs
  print "#\n# Make breakpoint raw VCFs\n#\n";
  for my $k (@kmers) {
    my $ctx = get_ctx($k);
    my $breakpoint_file    = "$proj/k$k/breakpoints/breakpoints.txt.gz";
    my $breakpoint_raw_vcf = "$proj/k$k/breakpoints/breakpoints.raw.vcf";

    print "# breakpoints raw VCF k=$k\n";
    print "$breakpoint_raw_vcf: $breakpoint_file \$(REF_FILE)\n";
    print "\t$ctx calls2vcf -o \$@ \$< \$(REF_FILE) >& \$@.log\n\n";
  }

  # Post-processing rules for VCFs
  print "#\n# Post-processing for raw VCFs\n#\n";
  print "$proj/%.sort.vcf: $proj/%.raw.vcf\n";
  print "\t\$(VCFSORT) \$< > \$@\n\n";

  print "$proj/%.norm.vcf: $proj/%.sort.vcf \$(REF_FILE)\n";
  print "\t\$(BCFTOOLS) norm --remove-duplicates --fasta-ref \$(REF_FILE) --multiallelics +both \$< | \\\n";
  print "\t\$(VCFRENAME) > \$@\n\n";

  # Generate union VCF
  print "#\n# Create union compressed VCF\n#\n";
  print "$union_bubble_vcf: \$(BUBBLE_VCFS) \$(BUBBLE_CSIS)\n";
  # print "\t\$(BCFTOOLS) concat --allow-overlaps --output-type z --output \$@ \$(BUBBLE_VCFS)\n\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BUBBLE_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";


  print "$union_brkpnt_vcf: \$(BREAKPOINT_VCFS) \$(BREAKPOINT_CSIS)\n";
  # print "\t\$(BCFTOOLS) concat --allow-overlaps --output-type z --output \$@ \$(BREAKPOINT_VCFS)\n\n";
  print "\t\$(BCFTOOLS) concat --allow-overlaps --remove-duplicates \$(BREAKPOINT_VCFS) | \\\n";
  print "\t\$(VCFRENAME) | \$(BCFTOOLS) view --output-type z --output-file \$@ -\n\n";

  print "#\n# General VCF rules\n#\n";
  # Compress a VCF
  print "%.vcf.gz: %.vcf\n";
  print "\t\$(BGZIP) -f \$<\n\n";

  # Create VCF index files .vcf.gz.csi
  print "%.vcf.gz.csi: %.vcf.gz\n";
  print "\t\$(BCFTOOLS) index -f \$<\n\n";
}


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
      if($sname !~ /^[a-z0-9_\-\.]+$/i) { print STDERR "Bad name: $sname"; exit(-1); }
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
  return map { get_ctx($_) } get_maxk_values(@_);
}

sub get_maxk_values
{
  my %maxks = ();
  for my $k (@_) { $maxks{(int(($k+31)/32) * 32 - 1)} = 1; }
  return keys %maxks;
}

sub get_ctx
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
