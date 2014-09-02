SHELL:=/bin/bash

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
BIOINF=$(CTXDIR)/libs/bioinf-perl
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
CTX=$(CTXDIR)/bin/ctx31

K=9
NSAMPLES=3
LAST_SAMP=2

SEQS=$(shell echo seq.{0..$(LAST_SAMP)}.fa)
POP_GRAPHS=$(SEQS:.fa=.k$(K).ctx)
POP_PATHS=$(SEQS:.fa=.k$(K).ctp.gz)
PLOTS=$(shell echo {supernodes,kmers}.{0..$(LAST_SAMP)}.k$(K).pdf)
CONTIGS=$(shell echo contigs.{0..$(LAST_SAMP)}.fa)
RMDUP_CONTIGS=$(shell echo rmdup.contigs.{0..$(LAST_SAMP)}.fa)

all: test

seq.%.fa:
	$(DNACAT) -F -n 1001 > $@

seq.%.k$(K).ctx: seq.%.fa
	$(CTX) build -m 1M -k $(K) --sample MrsSample$* --seq $< $@
	$(CTX) inferedges $@

seq.%.k$(K).ctp.gz: seq.%.k$(K).ctx seq.%.fa
	$(CTX) thread -m 1M --seq seq.$*.fa --out $@ seq.$*.k$(K).ctx

pop.ctx: $(POP_GRAPHS)
	$(CTX) join --out $@ $(POP_GRAPHS)

pop.ctp.gz: $(POP_PATHS)
	$(CTX) pjoin --out $@ $(POP_PATHS)

contigs.%.fa: pop.ctx pop.ctp.gz
	$(CTX) contigs --out $@ --colour $* -p pop.ctp.gz pop.ctx

rmdup.contigs.%.fa: contigs.%.fa
	$(CTX) rmsubstr -k $(K) $< > $@

supernodes.%.k$(K).pdf: seq.%.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

kmers.%.k$(K).pdf: seq.%.k$(K).ctx
	$(CTX2DOT) $< | dot -Tpdf > $@

plots: $(PLOTS)

test: $(CONTIGS) $(RMDUP_CONTIGS) $(SEQS)
	set -euo pipefail; \
	for i in {0..$(LAST_SAMP)}; do \
		echo \# Sample $$i; \
		$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.$$i.fa seq.$$i.fa; \
		$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 rmdup.contigs.$$i.fa seq.$$i.fa; \
	done;

clean:
	rm -rf $(SEQS) $(POP_GRAPHS) $(POP_PATHS)
	rm -rf pop.ctx pop.ctp.gz $(CONTIGS) $(RMDUP_CONTIGS)

.PHONY: all clean test plots
