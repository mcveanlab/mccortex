SHELL:=/bin/bash -euo pipefail

#
# Isaac Turner 2014-11-25
#
# In this test we thread a read through a loop to see how links allow us
# to reconstruct repetitive regions.
#

K=9
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

TGTS=seq.fa seq.fa seq.k$(K).ctx seq.k$(K).ctp.gz seq.regen.fa

all: $(TGTS) check

# CGTTAATCCG (ATTGAGTTCA x 4) TTGGTCACGG
seq.fa:
	echo "CGTTAATCCT""ATTGAGTTCA""ATTGAGTTCA""ATTGAGTTCA""ATTGAGTTCA""CTGGTCACGG" > $@

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -m 10M -k $(K) --sample LoopTest --seq $< $@ >& $@.log

seq.k$(K).ctp.gz: seq.k$(K).ctx seq.fa
	$(MCCORTEX) thread --seq seq.fa --out $@ seq.k$(K).ctx >& $@.log

seq.regen.fa: seq.k$(K).ctx seq.k$(K).ctp.gz
	( $(MCCORTEX) contigs --no-missing-check -o - -p seq.k$(K).ctp.gz seq.k$(K).ctx | \
	  $(MCCORTEX) rmsubstr -n 1000 -k $(K) - > $@ ) >& $@.log

check: seq.fa seq.regen.fa
	@echo Checking reconstructed sequence matches...
	diff -q <($(DNACAT) -P -k -r seq.fa) <($(DNACAT) -P -k -r seq.regen.fa)
	@echo '=> Sequences match.'

view: seq.k$(K).ctp.gz seq.fa seq.regen.fa
	gzip -dc seq.k$(K).ctp.gz
	cat seq.fa
	cat seq.regen.fa

seq.k$(K).pdf: seq.k$(K).ctx
	$(CTX2DOT) $< | dot -Tpdf > $@

plots: seq.k$(K).pdf

clean:
	rm -rf $(TGTS) seq.k$(K).pdf *.log

.PHONY: all clean plots check view
