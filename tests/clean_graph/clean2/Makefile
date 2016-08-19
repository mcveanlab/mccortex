SHELL:=/bin/bash -euo pipefail

K=17
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
PLOTCTX=$(CTXDIR)/scripts/cortex_to_graphviz.pl

GRAPHS=seq.k$(K).raw.ctx seq.k$(K).clean.ctx
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

seq.k$(K).raw.ctx: seq.fa
	$(MCCORTEX) build -q -m 10M -k $(K) --sample SeqJr --seq $< $@
	$(MCCORTEX) check -q $@

seq.k$(K).clean.ctx: seq.k$(K).raw.ctx
	$(MCCORTEX) clean -q --unitigs=2 -o $@ $<
	$(MCCORTEX) check -q $@

%.unitigs.dot: %.ctx
	$(MCCORTEX) unitigs -q -m 1M --dot $< > $@
	# $(PLOTCTX) --simplify $< > $@

%.dot: %.ctx
	$(PLOTCTX) $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

clean:
	rm -rf seq.fa $(GRAPHS) $(DOTS) $(PLOTS)

.PHONY: all plots clean
