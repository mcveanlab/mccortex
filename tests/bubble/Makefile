SHELL:=/bin/bash -euo pipefail

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
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

SEQS=seq0.fa seq1.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
FLAT_GRAPHS=$(shell printf "0:%s\n" $(GRAPHS))
BUBBLES=bubbles.seq.txt bubbles.join.txt
PLOTS=seq.k$(K).pdf
TGTS=$(SEQS) $(GRAPHS) join.k$(K).ctx $(BUBBLES)

all: $(TGTS) compare

test:
	echo $(SEQS)
	echo $(GRAPHS)

ref.fa:
	echo '>chr1'; echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA > $@

seq0.fa:
	echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA > $@
	echo CCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA >> $@

seq1.fa:
	echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGAT > $@

seq%.k$(K).ctx: seq%.fa
	$(CTX) build -m 10M -k $(K) --sample $< --seq $< $@

join.k$(K).ctx: $(GRAPHS)
	$(CTX) join -o $@ $(FLAT_GRAPHS)

bubbles.seq.txt: $(GRAPHS)
	$(CTX) bubbles -t 1 -m 10M -o bubbles.seq.txt.gz $(GRAPHS)
	gzip -d bubbles.seq.txt.gz
	cat bubbles.seq.txt

bubbles.join.txt: join.k$(K).ctx
	$(CTX) bubbles -t 1 -m 10M -o bubbles.join.txt.gz $(GRAPHS)
	gzip -d bubbles.join.txt.gz
	cat bubbles.join.txt

compare: bubbles.join.txt bubbles.seq.txt
	L1=`wc -l bubbles.join.txt | grep -oE '[0-9]+'`; \
	L2=`wc -l bubbles.seq.txt | grep -oE '[0-9]+'`; \
	echo join=$$L1; echo seq=$$L2; [ $$L1 == $$L2 ]
	@echo bubbles.join.txt, bubbles.seq.txt have same number of lines

seq.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

bubbles.vcf: ref.fa bubbles.join.txt
	bwa index ref.fa
	cat bubbles.join.txt | awk -F '[ \t]' 'm{print $$0;m=0;} /^>bubble\..*\.5pflank/{print $$1; m=1;}' > bubbles.join.txt.flanks
	bwa bwasw ref.fa bubbles.join.txt.flanks > bubbles.join.txt.sam
	$(CTX) calls2vcf -F bubbles.join.txt.sam -o $@ bubbles.join.txt ref.fa

decompose: bubbles.vcf

plots: $(PLOTS)

clean:
	rm -rf $(TGTS) $(PLOTS) bubbles.seq.txt.gz bubbles.join.txt.gz
	rm -rf ref.fa* bubbles.vcf bubbles.join.txt.flanks bubbles.join.txt.sam

.PHONY: all clean compare plots decompose
