
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
PLOTCTX=$(CTXDIR)/scripts/cortex_to_graphviz.pl

KMER=17

GRAPHS=seq.k$(KMER).raw.ctx seq.k$(KMER).clean.ctx
DOTS=$(GRAPHS:.ctx=.dot) $(GRAPHS:.ctx=.supernodes.dot)
PLOTS=$(DOTS:.dot=.pdf)

all: seq.fa $(GRAPHS) $(DOTS)

plots: $(PLOTS)

seq.fa: Makefile
	# Repeat of CAAAGGCCTCACGGGTA
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTACAAAGGCCTCACGGGTAGAATCCCCTTTG > seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo AAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA >> seq.fa

seq.k$(KMER).raw.ctx: seq.fa
	$(CTX) build -m 10M -k $(KMER) --sample SeqJr --seq $< $@
	$(CTX) check -q $@

seq.k$(KMER).clean.ctx: seq.k$(KMER).raw.ctx
	$(CTX) clean --supernodes=2 -o $@ $<
	$(CTX) check -q $@

%.supernodes.dot: %.ctx
	$(CTX) supernodes -m 1M --graphviz $< > $@
	# $(PLOTCTX) --simplify $< > $@

%.dot: %.ctx
	$(PLOTCTX) $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

clean:
	rm -rf seq.fa $(GRAPHS) $(DOTS) $(PLOTS)

.PHONY: all plots clean
