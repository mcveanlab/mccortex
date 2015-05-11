SHELL=/bin/bash -euo pipefail

#
# Smash up a genome and count how many breakpoints we recover
#

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
READSIM=$(CTXDIR)/libs/readsim/readsim
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
SIMBRK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_breakpoints.pl
CHKBRK=$(CTXDIR)/scripts/check-breakpoints.pl
MKCIRCOS=$(CTXDIR)/scripts/make-circos.pl
CIRCOS=circos

K=31
NBREAKS=100
HAPDEPTH=30
READLEN=150
ERRRATE=0.01
MEM=1G
OUTDIR=proj
REF=../data/chr22/chr22_17M_18M.fa

DIRS=reads diploid
READS=reads/chrom0.$(HAPDEPTH)X.1.fa.gz reads/chrom0.$(HAPDEPTH)X.2.fa.gz \
      reads/chrom1.$(HAPDEPTH)X.1.fa.gz reads/chrom1.$(HAPDEPTH)X.2.fa.gz

# Mark all dependencies as secondary
# It means don't re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:
# Delete files if their recipe fails
.DELETE_ON_ERROR:
# Remove implicit rules for certain suffixes
.SUFFIXES:

all: run-mccortex

# 1 Mb human diploid
# Generate a diploid genome from a diploid reference
diploid/truth.txt: diploid/chrom1.fa
diploid/chrom1.fa:
	mkdir -p diploid
	$(SIMBRK) $(NBREAKS) $(REF) $@ diploid/truth.txt

# Simulate PE reads of each chrom each 50X
reads/chrom1.$(HAPDEPTH)X.1.fa.gz: reads/chrom1.$(HAPDEPTH)X.2.fa.gz
reads/chrom1.$(HAPDEPTH)X.2.fa.gz: diploid/chrom1.fa reads
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom1.$(HAPDEPTH)X

reads/chrom0.$(HAPDEPTH)X.1.fa.gz: reads/chrom0.$(HAPDEPTH)X.2.fa.gz
reads/chrom0.$(HAPDEPTH)X.2.fa.gz: $(REF) reads
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom0.$(HAPDEPTH)X

samples.haploid.txt:
	printf "MissSample . %s:%s,%s:%s\n" \
	  reads/chrom0.$(HAPDEPTH)X.1.fa.gz reads/chrom0.$(HAPDEPTH)X.2.fa.gz \
    reads/chrom1.$(HAPDEPTH)X.1.fa.gz reads/chrom1.$(HAPDEPTH)X.2.fa.gz > $@

samples.diploid.txt:
	printf "MissSample . %s:%s,%s:%s\n" \
	  reads/chrom0.$(HAPDEPTH)X.1.fa.gz reads/chrom0.$(HAPDEPTH)X.2.fa.gz \
    reads/chrom1.$(HAPDEPTH)X.1.fa.gz reads/chrom1.$(HAPDEPTH)X.2.fa.gz > $@

%.k$(K).mk: samples.%.txt
	$(CTXPIPELINE) -r $(REF) $(K) run_$* $< > $@

run-mccortex: haploid.k$(K).mk diploid.k$(K).mk $(READS) diploid/truth.txt
	$(MAKE) -f haploid.k$(K).mk CTXDIR=$(CTXDIR) MEM=$(MEM) breakpoints
	$(MAKE) -f diploid.k$(K).mk CTXDIR=$(CTXDIR) MEM=$(MEM) breakpoints
	$(CHKBRK) diploid/truth.txt <(gzip -fcd run_haploid/k$(K)/breakpoints/breakpoints.txt.gz)
	$(CHKBRK) diploid/truth.txt <(gzip -fcd run_diploid/k$(K)/breakpoints/breakpoints.txt.gz)

$(DIRS):
	mkdir -p $@

plot-circos: run-mccortex
	gzip -fcd run_haploid/k$(K)/breakpoints/breakpoints.txt.gz | $(MKCIRCOS) circos-plot -
	cd circos-plot && $(CIRCOS)

clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt haploid.k$(K).mk diploid.k$(K).mk circos-plot

.PHONY: all clean run-mccortex run-cortex plot-circos
