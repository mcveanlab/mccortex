SHELL:=/bin/bash -euo pipefail

CTXDIR=../../..
CTX=$(CTXDIR)/bin/mccortex31
PLOTCTX=$(CTXDIR)/scripts/cortex_to_graphviz.pl
KMER=31

KEEP=seq.fa seq.k$(KMER).raw.ctx
PLOTS=seq.k$(KMER).raw.pdf

all: $(KEEP)

seq.fa: Makefile
	# echo GGCCTCGTGATCCGCCCGCCTCGGCCTCCCACA > seq.fa
	echo GATATATATCTAGATATATATCTATATATAAAA > seq.fa

seq.k$(KMER).raw.ctx: seq.fa
	$(CTX) build -m 10M -k $(KMER) --sample SeqJr --seq $< $@

seq.k$(KMER).%.dot: seq.k$(KMER).%.ctx
	$(PLOTCTX) $< > $@

seq.k$(KMER).%.pdf: seq.k$(KMER).%.dot
	dot -Tpdf $< > $@

plots: $(PLOTS)

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all plot clean
