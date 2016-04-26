SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31
PLOTCTX=$(CTXDIR)/scripts/cortex_to_graphviz.pl

KMER=17

GRAPHS=seq.k$(KMER).raw.ctx seq.k$(KMER).clean.ctx
DOTS=$(GRAPHS:.ctx=.dot) $(GRAPHS:.ctx=.unitigs.dot)
PLOTS=$(DOTS:.dot=.pdf)

all: seq.fa $(GRAPHS) $(DOTS)

plots: $(PLOTS)

seq.fa: Makefile
	# Repeat of CAAAGGCCTCACGGGTA
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTACAAAGGCCTCACGGGTAGAATCCCCTTTG > seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo AAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA >> seq.fa

seq.k$(KMER).raw.ctx: seq.fa
	$(CTX) build -m 10M -k $(KMER) --sample SeqJr --seq $< $@
	$(CTX) check -q $@

seq.k$(KMER).clean.ctx: seq.k$(KMER).raw.ctx
	$(CTX) clean --unitigs=2 -o $@ $<
	$(CTX) check -q $@

%.unitigs.dot: %.ctx
	$(CTX) unitigs -m 1M --dot $< > $@
	# $(PLOTCTX) --simplify $< > $@

%.dot: %.ctx
	$(PLOTCTX) $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

clean:
	rm -rf seq.fa $(GRAPHS) $(DOTS) $(PLOTS)

.PHONY: all plots clean
