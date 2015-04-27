SHELL=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
READSIM=$(CTXDIR)/libs/readsim/readsim
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

K=31
HAPDEPTH=30
READLEN=100
ERRRATE=0.01
MEM=1G
OUTDIR=proj
REF=../data/chr22/chr22_17M_18M.fa

CHROMS=diploid/genome0.fa diploid/genome1.fa
READS=reads/chrom0.$(HAPDEPTH)X.1.fa.gz reads/chrom0.$(HAPDEPTH)X.2.fa.gz \
			reads/chrom1.$(HAPDEPTH)X.1.fa.gz reads/chrom1.$(HAPDEPTH)X.2.fa.gz
DIRS=reads diploid

.INTERMEDIATE:

all: pipeline

test:
	@echo $(READS)
	@echo reads/chrom0.$(HAPDEPTH)X.1.fa.gz

# 1 Mb human diploid

# Generate a diploid genome
diploid/genome1.fa: diploid/genome0.fa | $(DIRS)
diploid/genome0.fa: $(REF) | $(DIRS)
	$(SIMMUT) --snps 1000 --indels 100 --invs 0 diploid 2 $<

# Remove deletion marks (-) and convert to uppercase
diploid/chrom%.fa: diploid/genome%.fa
	cat $< | tr -d '-' | $(DNACAT) -u -F - > $@

# Simulate PE reads of each chrom each 50X
reads/chrom0.$(HAPDEPTH)X.1.fa.gz: reads/chrom0.$(HAPDEPTH)X.2.fa.gz
reads/chrom1.$(HAPDEPTH)X.1.fa.gz: reads/chrom1.$(HAPDEPTH)X.2.fa.gz

reads/chrom%.$(HAPDEPTH)X.2.fa.gz: diploid/chrom%.fa
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom$*.$(HAPDEPTH)X

samples.txt:
	echo "MissSample . reads/chrom0.$(HAPDEPTH)X.1.fa.gz:reads/chrom1.$(HAPDEPTH)X.2.fa.gz" > $@

bubbles.mk: samples.txt
	$(CTXPIPELINE) -r $(REF) $(K) $(OUTDIR) $< > $@

pipeline: bubbles.mk $(READS)
	$(MAKE) -f $< CTXDIR=$(CTXDIR)

# 5. Compare

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt bubbles.mk

.PHONY: all clean pipeline
