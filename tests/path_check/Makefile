SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
K=7
KEEP=$(shell echo genome.fa genome.k$(K).{ctx,ctp} reads.1.fa.gz reads.2.fa.gz reads.{se,pe}.k$(K).ctp)
PLOTS=genome.k$(K).dot genome.k$(K).pdf

all: $(KEEP)

plots: $(PLOTS)

clean:
	rm -rf $(KEEP) $(PLOTS)

# Sample random genome
genome.fa:
	$(DNACAT) -F -n 200 > genome.fa

genome.k$(K).ctx: genome.fa
	$(MCCORTEX) build -m 10M -k $(K) --sample MssrGenome --seq $< $@

genome.k$(K).dot: genome.k$(K).ctx
	$(MCCORTEX) unitigs --dot $< > $@

genome.k$(K).pdf: genome.k$(K).dot
	dot -Tpdf $< > $@

reads.1.fa.gz reads.2.fa.gz:
	../../libs/readsim/readsim -r genome.fa -l 10 -i 20 -v 0.1 -d 2 reads

genome.k$(K).ctp: genome.k$(K).ctx
	$(MCCORTEX) thread -m 10M -t 1 --seq genome.fa -o $@ $<

reads.se.k$(K).ctp: genome.k$(K).ctx reads.1.fa.gz reads.2.fa.gz
	$(MCCORTEX) thread -m 10M -t 1 --seq reads.1.fa.gz --seq reads.2.fa.gz -o reads.se.k$(K).ctp genome.k$(K).ctx

reads.pe.k$(K).ctp: genome.k$(K).ctx reads.1.fa.gz reads.2.fa.gz
	$(MCCORTEX) thread -m 10M -t 1 --seq2 reads.1.fa.gz:reads.2.fa.gz -o reads.pe.k$(K).ctp genome.k$(K).ctx

.PHONY: all plots clean
