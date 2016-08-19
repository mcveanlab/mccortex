SHELL:=/bin/bash -euo pipefail

# Might work with diff K values

K=9
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
BIOINF=$(CTXDIR)/libs/bioinf-perl
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

NSAMPLES=3
LAST_SAMP=2

SEQS=$(shell echo seq.{0..$(LAST_SAMP)}.fa)
POP_GRAPHS=$(SEQS:.fa=.k$(K).ctx)
POP_PATHS=$(SEQS:.fa=.k$(K).ctp.gz)
POP_PATHS_CSV=$(SEQS:.fa=.k$(K).contig-hist.csv)
CONFID_CSV=$(SEQS:.fa=.k$(K).confid.csv)
PLOTS=$(shell echo {unitigs,kmers}.{0..$(LAST_SAMP)}.k$(K).pdf)
CONTIGS=$(shell echo contigs.{0..$(LAST_SAMP)}.fa)
RMDUP_CONTIGS=$(shell echo rmdup.contigs.{0..$(LAST_SAMP)}.fa)

GENOME=1001

all: test

seq.%.fa:
	$(DNACAT) -F -n $(GENOME) > $@

seq.%.k$(K).ctx: seq.%.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample MrsSample$* --seq $< $@ >& $@.log
	$(MCCORTEX) inferedges -q $@

seq.%.k$(K).ctp.gz: seq.%.k$(K).ctx seq.%.fa
	$(MCCORTEX) thread -m 1M --seq seq.$*.fa --out $@ seq.$*.k$(K).ctx >& $@.log

pop.k$(K).ctx: $(POP_GRAPHS)
	$(MCCORTEX) join -q --out $@ $(POP_GRAPHS)

pop.k$(K).ctp.gz: $(POP_PATHS)
	$(MCCORTEX) pjoin -q --out $@ $(POP_PATHS)

contigs.%.fa: pop.k$(K).ctx pop.k$(K).ctp.gz
	$(MCCORTEX) contigs --use-seed-paths --no-missing-check --out $@ --colour $* --genome $(GENOME) --confid-csv seq.$*.k$(K).confid.csv -p pop.k$(K).ctp.gz pop.k$(K).ctx >& $@.log

rmdup.contigs.%.fa: contigs.%.fa
	$(MCCORTEX) rmsubstr -q -k $(K) $< > $@

unitigs.%.k$(K).pdf: seq.%.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

kmers.%.k$(K).pdf: seq.%.k$(K).ctx
	$(CTX2DOT) $< | dot -Tpdf > $@

plots: $(PLOTS)

test: $(CONTIGS) $(RMDUP_CONTIGS) $(SEQS)
	for i in {0..$(LAST_SAMP)}; do \
		echo \# Sample $$i; \
		$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.$$i.fa seq.$$i.fa; \
		$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 rmdup.contigs.$$i.fa seq.$$i.fa; \
	done;

clean:
	rm -rf $(SEQS) $(POP_GRAPHS) $(POP_PATHS) $(POP_PATHS_CSV) $(CONFID_CSV)
	rm -rf pop.k$(K).ctx pop.k$(K).ctp.gz $(CONTIGS) $(RMDUP_CONTIGS) *.log

.PHONY: all clean test plots
