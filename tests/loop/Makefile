SHELL:=/bin/bash -euo pipefail

#
# Isaac Turner 2014-11-25
#
# In this test we thread a read through a loop to see how links allow us
# to reconstruct repetitive regions.
#

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

K=9
TGTS=seq.fa seq.fa seq.k$(K).ctx seq.k$(K).ctp.gz seq.regen.fa

all: $(TGTS) check

# CGTTAATCCG (ATTGAGTTCA x 4) TTGGTCACGG
seq.fa:
	echo "CGTTAATCCT""ATTGAGTTCA""ATTGAGTTCA""ATTGAGTTCA""ATTGAGTTCA""CTGGTCACGG" > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -m 10M -k $(K) --sample LoopTest --seq $< $@

seq.k$(K).ctp.gz: seq.k$(K).ctx seq.fa
	$(CTX) thread --seq seq.fa --out $@ seq.k$(K).ctx

seq.regen.fa: seq.k$(K).ctx seq.k$(K).ctp.gz
	$(CTX) contigs --no-missing-check -o - -p seq.k$(K).ctp.gz seq.k$(K).ctx | \
	$(CTX) rmsubstr -n 1000 -k $(K) - > $@

check:
	gzip -dc seq.k$(K).ctp.gz
	@echo; echo Checking reconstructed sequence matches
	diff -q <(dnacat -P -k -r seq.fa) <(dnacat -P -k -r seq.regen.fa)

seq.k$(K).pdf: seq.k$(K).ctx
	$(CTX2DOT) $< | dot -Tpdf > $@

plots: seq.k$(K).pdf

clean:
	rm -rf $(TGTS) seq.k$(K).pdf

.PHONY: all clean plots check
