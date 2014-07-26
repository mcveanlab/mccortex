
#
# 1 het site in s0, 1 het between s0,s1
#
#                 x                  y
# s0: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
# s0: CCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
#
# s1: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGAT
#

SHELL:=/bin/bash

K=11
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

SEQS=seq0.fa seq1.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
BUBBLES=bubbles.seq.txt bubbles.join.txt
PLOTS=seq.k$(K).pdf
TGTS=$(SEQS) $(GRAPHS) join.k$(K).ctx $(BUBBLES)

all: $(TGTS) compare

test:
	echo $(SEQS)
	echo $(GRAPHS)

seq0.fa:
	echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA > seq0.fa
	echo CCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA >> seq0.fa

seq1.fa:
	echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGAT > seq1.fa

seq%.k$(K).ctx: seq%.fa
	$(CTX) build -m 10M -k $(K) --sample $< --seq $< $@

join.k$(K).ctx: $(GRAPHS)
	$(CTX) join --flatten -o $@ $(GRAPHS)

bubbles.seq.txt: $(GRAPHS)
	$(CTX) bubbles -t 1 -m 10M -o bubbles.seq.txt.gz $(GRAPHS)
	gzip -d bubbles.seq.txt.gz
	cat bubbles.seq.txt

bubbles.join.txt: join.k$(K).ctx
	$(CTX) bubbles -t 1 -m 10M -o bubbles.join.txt.gz $(GRAPHS)
	gzip -d bubbles.join.txt.gz
	cat bubbles.join.txt

compare: bubbles.join.txt bubbles.seq.txt
	set -euo pipefail; L1=`wc -l bubbles.join.txt | grep -oE '[0-9]+'`; \
	L2=`wc -l bubbles.seq.txt | grep -oE '[0-9]+'`; \
	echo join=$$L1; echo seq=$$L2; [ $$L1 == $$L2 ]
	@echo bubbles.join.txt, bubbles.seq.txt have same number of lines

seq.k$(K).pdf: join.k$(K).ctx
	set -euo pipefail; $(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: $(PLOTS)

clean:
	rm -rf $(TGTS) $(PLOTS)

.PHONY: all clean plots compare
