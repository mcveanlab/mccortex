
CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
PLOTCTX=$(CTXPATH)/scripts/cortex_to_graphviz.pl

KMER=17

KEEP=seq.fa seq.k$(KMER).raw.ctx seq.k$(KMER).clean.ctx \
     seq.k$(KMER).supernodes.clean.dot seq.k$(KMER).supernodes.raw.dot

PLOTS=seq.k$(KMER).clean.pdf seq.k$(KMER).raw.pdf \
      seq.k$(KMER).supernodes.clean.pdf seq.k$(KMER).supernodes.raw.pdf

all: $(KEEP)

seq.fa: Makefile
	# Repeat of CAAAGGCCTCACGGGTA
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTACAAAGGCCTCACGGGTAGAATCCCCTTTG > seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo GTGAGGCCAAGCAAAGGCCTCACGGGTAGAATCCCCTTTG >> seq.fa
	echo AAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA >> seq.fa

seq.k$(KMER).raw.ctx: seq.fa
	$(CTX) build -m 10M -k $(KMER) --sample SeqJr --seq $< $@

seq.k$(KMER).clean.ctx: seq.k$(KMER).raw.ctx
	$(CTX) clean --supernodes=2 -o $@ $<

seq.k$(KMER).supernodes.%.dot: seq.k$(KMER).%.ctx
	$(PLOTCTX) --simplify $< > $@

seq.k$(KMER).%.dot: seq.k$(KMER).%.ctx
	$(PLOTCTX) $< > $@

seq.k$(KMER).%.pdf: seq.k$(KMER).%.dot
	dot -Tpdf $< > $@

plots: $(PLOTS)

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all plots clean
