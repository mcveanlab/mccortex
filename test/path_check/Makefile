
SHELL:=/bin/bash

CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
SEQRND=$(CTXPATH)/libs/seq_file/bin/seqrnd
K=7
GENOME=200
KEEP=$(shell echo genome.fa genome.k$(K).{ctx,ctp} reads.1.fa.gz reads.2.fa.gz reads.{se,pe}.k$(K).ctp)
PLOTS=genome.k$(K).dot genome.k$(K).pdf

all: $(KEEP)

plots: $(PLOTS)

clean:
	rm -rf $(KEEP) $(PLOTS)

# Sample random genome
genome.fa:
	$(SEQRND) $(GENOME) > genome.fa

genome.k$(K).ctx: genome.fa
	$(CTX) build -m 10M -k $(K) --sample MssrGenome --seq $< $@

genome.k$(K).dot: genome.k$(K).ctx
	$(CTX) supernodes --graphviz $< > $@

genome.k$(K).pdf: genome.k$(K).dot
	dot -Tpdf $< > $@

reads.1.fa.gz reads.2.fa.gz:
	../../libs/readsim/readsim -r genome.fa -l 10 -i 20 -v 0.1 -d 2 reads

genome.k$(K).ctp: genome.k$(K).ctx
	$(CTX) thread -m 10M -t 1 --seq genome.fa -o $@ $<

reads.se.k$(K).ctp: genome.k$(K).ctx reads.1.fa.gz reads.2.fa.gz
	$(CTX) thread -m 10M -t 1 --seq reads.1.fa.gz --seq reads.2.fa.gz -o reads.se.k$(K).ctp genome.k$(K).ctx

reads.pe.k$(K).ctp: genome.k$(K).ctx reads.1.fa.gz reads.2.fa.gz
	$(CTX) thread -m 10M -t 1 --seq2 reads.1.fa.gz:reads.2.fa.gz -o reads.pe.k$(K).ctp genome.k$(K).ctx

.PHONY: all plots clean
