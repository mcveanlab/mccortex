SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

K=7
PLOTS=genome.k$(K).ctx.pdf genome.k$(K).perl.pdf
KEEP=genome.fa genome.k$(K).ctx $(PLOTS:.pdf=.dot)

all: $(KEEP)

clean:
	rm -rf $(KEEP) $(PLOTS) $(PLOTS:.pdf=.dot)

# Sample random genome
genome.fa:
	$(DNACAT) -F -n 200 > genome.fa

genome.k$(K).ctx: genome.fa
	$(CTX) build -m 1M -k $(K) --sample MssrGenome --seq $< $@

genome.k$(K).ctx.dot: genome.k$(K).ctx
	$(CTX) supernodes -m 1M --graphviz --points $< > $@

genome.k$(K).perl.dot: genome.k$(K).ctx
	$(CTX2DOT) --simplify $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

plots: $(PLOTS)

.PHONY: all clean plots
