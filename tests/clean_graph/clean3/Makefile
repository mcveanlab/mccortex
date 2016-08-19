SHELL:=/bin/bash -euo pipefail

K=31
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
PLOTCTX=$(CTXDIR)/scripts/cortex_to_graphviz.pl

KEEP=seq.fa seq.k$(K).raw.ctx
PLOTS=seq.k$(K).raw.pdf

all: $(KEEP)

seq.fa: Makefile
	# echo GGCCTCGTGATCCGCCCGCCTCGGCCTCCCACA > seq.fa
	echo GATATATATCTAGATATATATCTATATATAAAA > seq.fa

seq.k$(K).raw.ctx: seq.fa
	$(MCCORTEX) build -q -m 10M -k $(K) --sample SeqJr --seq $< $@

seq.k$(K).%.dot: seq.k$(K).%.ctx
	$(PLOTCTX) $< > $@

seq.k$(K).%.pdf: seq.k$(K).%.dot
	dot -Tpdf $< > $@

plots: $(PLOTS)

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all plot clean
