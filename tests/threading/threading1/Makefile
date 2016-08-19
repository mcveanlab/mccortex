SHELL:=/bin/bash -euo pipefail

K=31
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CTX2GV=$(CTXDIR)/scripts/cortex_to_graphviz.pl

GRAPHS=ref.k$(K).ctx ref.ctp.gz
LOGS=$(addsuffix .log,$(GRAPHS))
PLOTS=ref.pdf ref.dot
TGTS=ref.fa $(GRAPHS)

all: $(TGTS)

plots: ref.dot ref.pdf

ref.fa:
	echo CGATTGAATTCCACCGATAATGCAGATGTGAGCCTCAGCATCTACTGCTTCCTCGTCGTCGGGGACTTTTGTTGACCTACCACATGATACATGCGGCCAT >  $@
	echo ACAAGCTAAAGAAGCTAGCCAGTGCAGGCTCCCTTCAGCATCTACTGCTTCCTCGTCGTCGGGGACTAGAAACGTGACCATCGGCCACCGAAAGATAAGG >> $@
	echo ACAAGCTAAAGAAGCTAGCCAGTGCAGGCTCCCTTCAGCATCTACTGCTGGCTCGTCGTCGGGGACTAGAAACGTGACCATCGGCCACCGAAAGATAAGG >> $@

ref.k$(K).ctx: ref.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample MsSample --seq ref.fa ref.k$(K).ctx >& $@.log

ref.ctp.gz: ref.k$(K).ctx
	$(MCCORTEX) thread -m 1M --seq ref.fa -o $@ $< >& $@.log

ref.dot: ref.k$(K).ctx
	$(CTX2GV) --simplify $< > $@

ref.pdf: ref.dot
	dot -Tpdf $< > $@

clean:
	rm -rf $(TGTS) $(PLOTS) $(LOGS)

.PHONY: all plots clean
