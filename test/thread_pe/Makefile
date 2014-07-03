
CTX=../../bin/ctx31
PATHS=reads.se.one.ctp.gz reads.se.two.ctp.gz \
      reads.pe.one.ctp.gz reads.pe.two.ctp.gz \
      join.ctp.gz

TGTS=genome.fa read.se.fa read.1.fa read.2.fa genome.k9.ctx $(PATHS)

SHELL:=/bin/bash

all: $(TGTS) check

plots: genome.k9.pdf

clean:
	rm -rf $(TGTS) *.csv genome.k9.pdf

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

genome.k9.ctx: genome.fa
	$(CTX) build -k 9 --sample Genome --seq genome.fa --sample Gen2 --seq genome.fa genome.k9.ctx

reads.se.one.ctp.gz: genome.k9.ctx read.se.fa
	$(CTX) thread -m 1M --frag-hist frag-hist.se1.csv --gap-hist gap-hist.se1.csv --print-contigs --print-paths --one-way --seq read.se.fa -o $@ genome.k9.ctx:0
	$(CTX) check -p $@ genome.k9.ctx:0

reads.se.two.ctp.gz: genome.k9.ctx read.se.fa
	$(CTX) thread -m 1M --frag-hist frag-hist.se2.csv --gap-hist gap-hist.se2.csv --print-contigs --print-paths --two-way --seq read.se.fa -o $@ genome.k9.ctx:0
	$(CTX) check -p $@ genome.k9.ctx:0

reads.pe.one.ctp.gz: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --frag-hist frag-hist.pe1.csv --gap-hist gap-hist.pe1.csv --print-contigs --print-paths --one-way --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0
	$(CTX) check -p $@ genome.k9.ctx:0

reads.pe.two.ctp.gz: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --frag-hist frag-hist.pe2.csv --gap-hist gap-hist.pe2.csv --print-contigs --print-paths --two-way --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0
	$(CTX) check -p $@ genome.k9.ctx:0

join.ctp.gz: reads.se.one.ctp.gz reads.se.two.ctp.gz reads.pe.one.ctp.gz reads.pe.two.ctp.gz
	$(CTX) pjoin -m 1M --out $@ reads.{se,pe}.{one,two}.ctp.gz

genome.k9.pdf:
	../../scripts/cortex_to_graphviz.pl genome.k9.ctx:0 | dot -Tpdf > genome.k9.pdf

tmp.k9.ctx: 
	$(CTX) join -m 1M -o $@ genome.k9.ctx:0 genome.k9.ctx:0 genome.k9.ctx:0 genome.k9.ctx:0

check: join.ctp.gz tmp.k9.ctx
	$(CTX) check -p join.ctp.gz tmp.k9.ctx

# seq-gaps -> gap-hist
# mp-gaps -> frag-hist

.INTERMEDIATE: tmp.k9.ctx

.PHONY: all plots clean check
