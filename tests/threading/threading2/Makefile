SHELL:=/bin/bash -euo pipefail

K=9
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

LINKS=reads.se.one.k$(K).ctp.gz reads.se.two.k$(K).ctp.gz \
      reads.pe.one.k$(K).ctp.gz reads.pe.two.k$(K).ctp.gz \
      join.k$(K).ctp.gz

CSV_FILES=$(shell echo {frag,gap,contig}-hist.{se,pe}{1,2}.csv)

TGTS=genome.fa read.se.fa read.1.fa read.2.fa genome.k$(K).ctx $(LINKS)

all: $(TGTS) check

plots: genome.k$(K).pdf

clean:
	rm -rf $(TGTS) $(CSV_FILES) genome.k$(K).pdf tmp.k$(K).ctx

genome.fa:
	echo TCGGCATCAGTGGCCATA > genome.fa
	echo TCGTCATCAGTGGCCGTA >> genome.fa

read.se.fa:
	echo TCGGCATCAtGTGGCCATA > $@

read.1.fa:
	echo GTCATCAGTG > $@

read.2.fa:
	# echo CAGTGGCCGT > $@ # Actually revcmp read2
	echo ACGGCCACTG > $@

genome.k$(K).ctx: genome.fa
	$(MCCORTEX) build -k $(K) --sample Genome --seq genome.fa --sample Gen2 --seq genome.fa genome.k$(K).ctx

reads.se.one.k$(K).ctp.gz: genome.k$(K).ctx read.se.fa
	$(MCCORTEX) thread -m 1M --frag-hist frag-hist.se1.csv --gap-hist gap-hist.se1.csv --print-contigs --print-paths --one-way --seq read.se.fa -o $@ genome.k$(K).ctx:0
	$(MCCORTEX) check -p $@ genome.k$(K).ctx:0

reads.se.two.k$(K).ctp.gz: genome.k$(K).ctx read.se.fa
	$(MCCORTEX) thread -m 1M --frag-hist frag-hist.se2.csv --gap-hist gap-hist.se2.csv --print-contigs --print-paths --two-way --seq read.se.fa -o $@ genome.k$(K).ctx:0
	$(MCCORTEX) check -p $@ genome.k$(K).ctx:0

reads.pe.one.k$(K).ctp.gz: genome.k$(K).ctx read.1.fa read.2.fa
	$(MCCORTEX) thread -m 1M --frag-hist frag-hist.pe1.csv --gap-hist gap-hist.pe1.csv --print-contigs --print-paths --one-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(K).ctx:0
	$(MCCORTEX) check -p $@ genome.k$(K).ctx:0

reads.pe.two.k$(K).ctp.gz: genome.k$(K).ctx read.1.fa read.2.fa
	$(MCCORTEX) thread -m 1M --frag-hist frag-hist.pe2.csv --gap-hist gap-hist.pe2.csv --print-contigs --print-paths --two-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(K).ctx:0
	$(MCCORTEX) check -p $@ genome.k$(K).ctx:0

join.k$(K).ctp.gz: reads.se.one.k$(K).ctp.gz reads.se.two.k$(K).ctp.gz reads.pe.one.k$(K).ctp.gz reads.pe.two.k$(K).ctp.gz
	$(MCCORTEX) pjoin -m 1M --out $@ reads.{se,pe}.{one,two}.k$(K).ctp.gz

genome.k$(K).pdf:
	$(CTX2DOT) genome.k$(K).ctx:0 | dot -Tpdf > genome.k$(K).pdf

tmp.k$(K).ctx:
	$(MCCORTEX) join -m 1M -o $@ 0-3:genome.k$(K).ctx:0,0,0,0

check: join.k$(K).ctp.gz tmp.k$(K).ctx
	$(MCCORTEX) check -p join.k$(K).ctp.gz tmp.k$(K).ctx

# seq-gaps -> gap-hist
# mp-gaps -> frag-hist

.INTERMEDIATE: tmp.k$(K).ctx

.PHONY: all plots clean check
