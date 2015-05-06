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
HAPDEPTH=50
READLEN=150
ERRRATE=0.01
MEM=1G
OUTDIR=proj
REF=../data/chr22/chr22_17M_18M.fa

DIRS=reads haploid
READS=reads/chrom.$(HAPDEPTH)X.1.fa.gz reads/chrom.$(HAPDEPTH)X.2.fa.gz

# Mark all dependencies as secondary
# It means don't re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:
# Delete files if their recipe fails
.DELETE_ON_ERROR:
# Remove implicit rules for certain suffixes
.SUFFIXES:

all: run-mccortex $(TRUTH_FILES)

# 1 Mb human haploid
# Generate a haploid genome from a haploid reference
haploid/truth.txt: haploid/chrom.fa
haploid/chrom.fa:
	mkdir -p haploid
	$(SIMBRK) $(NBREAKS) $(REF) $@ haploid/truth.txt

# Simulate PE reads of each chrom each 50X
reads/chrom.$(HAPDEPTH)X.1.fa.gz: reads/chrom.$(HAPDEPTH)X.2.fa.gz
reads/chrom.$(HAPDEPTH)X.2.fa.gz: haploid/chrom.fa
	mkdir -p reads
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom.$(HAPDEPTH)X

samples.txt:
	echo "MissSample . reads/chrom.$(HAPDEPTH)X.1.fa.gz:reads/chrom.$(HAPDEPTH)X.2.fa.gz" > $@

task.k$(K).mk: samples.txt
	$(CTXPIPELINE) -r $(REF) $(K) $(OUTDIR) $< > $@

run-mccortex: task.k$(K).mk $(READS) haploid/truth.txt
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) breakpoints
	$(CHKBRK) haploid/truth.txt <(gzip -fcd proj/k$(K)/breakpoints/breakpoints.txt.gz)

$(DIRS):
	mkdir -p $@

plot-circos: run-mccortex
	gzip -fcd proj/k$(K)/breakpoints/breakpoints.txt.gz | $(MKCIRCOS) circos-plot -
	cd circos-plot && $(CIRCOS)


clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt task.k$(K).mk circos-plot

.PHONY: all clean run-mccortex run-cortex plot-circos
