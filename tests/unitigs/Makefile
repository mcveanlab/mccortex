SHELL:=/bin/bash -euo pipefail

#
# Test unitigs command by generating 200 random DNA bases, building cortex graph
# then generating untigs with various output options
#

K=7
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

FILES=genome.fa genome.k$(K).ctx
UNITIGS=genome.k$(K).unitigs.fa genome.k$(K).unitigs.dot genome.k$(K).unitigs.gfa
PLOTS=genome.k$(K).unitigs.dot genome.k$(K).kmers.dot
PDFS=$(PLOTS:.dot=.pdf)

TGTS=$(FILES) $(UNITIGS) $(PLOTS)

all: $(TGTS)

clean:
	rm -rf $(TGTS) $(PDFS)

# Sample random genome
genome.fa:
	$(DNACAT) -F -n 200 > genome.fa

genome.k$(K).ctx: genome.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample MssrGenome --seq $< $@

genome.k$(K).unitigs.fa: genome.k$(K).ctx
	$(MCCORTEX) unitigs -q -m 1M -o $@ $<

genome.k$(K).unitigs.dot: genome.k$(K).ctx
	$(MCCORTEX) unitigs -q -m 1M --dot --points $< > $@

genome.k$(K).unitigs.gfa: genome.k$(K).ctx
	$(MCCORTEX) unitigs -q -m 1M --gfa $< > $@

genome.k$(K).kmers.dot: genome.k$(K).ctx
	$(CTX2DOT) $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

plots: $(PDFS)

.PHONY: all clean plots
