
#
# 1 het site in s0, 1 het between s0,s1
#
#                 x                  y
# s0: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
# s0: CCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
#
# s1: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGAT
#

K=11
CTXDIR=../../
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

SEQS=seq0.fa seq1.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
TGTS=bubbles.txt $(GRAPHS) join.k$(K).ctx

all: $(TGTS)

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

bubbles.txt.gz: $(GRAPHS)
	$(CTX) bubbles -t 1 -m 10M -o bubbles.txt.gz $(GRAPHS)

bubbles.txt: bubbles.txt.gz
	gzip -df bubbles.txt.gz
	cat bubbles.txt

join.k$(K).ctx: $(GRAPHS)
	$(CTX) join --flatten -o $@ $(GRAPHS)

seq.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: seq.k$(K).pdf

clean:
	rm -rf $(TGTS) $(SEQS) bubbles.txt.gz seq.k$(K).pdf

.PHONY: all clean plots
