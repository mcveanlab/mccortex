SHELL:=/bin/bash -euo pipefail

K=9
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

GRAPHS=genome.k9.ctx
LINKS=reads.pe.one.ctp reads.pe.two.ctp
LOGS=$(addsuffix .log,$(GRAPHS) $(LINKS))
TGTS=genome.fa read.1.fa read.2.fa $(GRAPHS) $(LINKS)
# non-default target: genome.k9.pdf

all: $(TGTS)

clean:
	rm -rf $(TGTS) $(LOGS) gap_sizes.*.csv mp_sizes.*.csv genome.k$(K).pdf

plots: genome.k$(K).pdf

genome.fa:
	echo gCATCAGTGGCCttggcgactcgc > genome.fa
	echo TCATCAGTGGCCATGACGCTAACT >> genome.fa

read.1.fa:
	echo TCATCAGTGG > read.1.fa

read.2.fa:
	# echo ACGCTAACT > read.2.fa # Actually revcmp read2
	echo AGTTAGCGT > read.2.fa

genome.k$(K).ctx: genome.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample Genome --seq genome.fa genome.k$(K).ctx >& $@.log

reads.pe.one.ctp: genome.k$(K).ctx read.1.fa read.2.fa
	$(MCCORTEX) thread -m 1M --print-contigs --one-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(K).ctx >& $@.log

reads.pe.two.ctp: genome.k$(K).ctx read.1.fa read.2.fa
	$(MCCORTEX) thread -m 1M --print-contigs --two-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(K).ctx >& $@.log

genome.k$(K).pdf:
	$(CTX2DOT) genome.k$(K).ctx | dot -Tpdf > genome.k$(K).pdf

.PHONY: all clean plots
