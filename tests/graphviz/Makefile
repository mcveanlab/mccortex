SHELL:=/bin/bash -euo pipefail

K=15
CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)
CTX2DOT=$(CTXDIR)/scripts/perl/mccortex-graph-to-graphviz.pl

all: seq.k$(K).ctx seq.k$(K).unitigs.dot seq.k$(K).kmers.dot

plots: seq.k$(K).unitigs.pdf seq.k$(K).kmers.pdf

seq.fa:
	echo TACATTCCCCCATAGTCGTAGGCGTTAAATACA     > seq.fa
	echo TACATTCCCCCATAGTAGTAGGCGTTAAATACA    >> seq.fa
	echo GGCGTTCGCTTATCCGGATAAGCGAACGCC       >> seq.fa
	echo ATATATATATATATATATATATATATATATATATAT >> seq.fa
	echo ACTTCTTCGAAAAAAAAAAAAAAATACTGGCCCTAACTTCTTCGAAAAAA >> seq.fa

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -q -k $(K) --sample MasterSeq --seq seq.fa seq.k$(K).ctx

seq.k$(K).unitigs.dot: seq.k$(K).ctx
	$(MCCORTEX) unitigs -q --dot $< > $@

seq.k$(K).kmers.dot: seq.k$(K).ctx
	$(CTX2DOT) $< > $@

%.pdf: %.dot
	cat $< | dot -Tpdf > $@

clean:
	rm -rf seq.fa seq.k$(K).ctx
	rm -rf seq.k$(K).unitigs.dot seq.k$(K).kmers.dot
	rm -rf seq.k$(K).unitigs.pdf seq.k$(K).kmers.pdf

.PHONY: all clean plots
