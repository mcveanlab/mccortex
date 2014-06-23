
CTX=../../bin/ctx31
PATHS=reads.se.ctp.gz reads.pe.one.ctp.gz reads.pe.two.ctp.gz join.ctp.gz
TGTS=genome.fa read.1.fa read.2.fa genome.k9.ctx $(PATHS)

# non-default target: genome.k9.pdf

all: $(TGTS) check

plots: genome.k9.pdf

clean:
	rm -rf $(TGTS) gap_sizes.*.csv mp_sizes.*.csv genome.k9.pdf

genome.fa:
	echo GGCATCAGTGGCCAT > genome.fa
	echo GTCATCAGTGGCCGT >> genome.fa

read.1.fa:
	echo GTCATCAGTG > read.1.fa

read.2.fa:
	# echo CAGTGGCCGT > read.2.fa # Actually revcmp read2
	echo ACGGCCACTG > read.2.fa

genome.k9.ctx: genome.fa
	$(CTX) build -k 9 --sample Genome --seq genome.fa --sample Gen2 --seq genome.fa genome.k9.ctx

reads.se.ctp.gz: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths --seq read.1.fa --seq read.2.fa -o $@ genome.k9.ctx:0

reads.pe.one.ctp.gz: reads.se.ctp.gz genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths -p reads.se.ctp.gz -C --oneway --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0

reads.pe.two.ctp.gz: reads.se.ctp.gz genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths -p reads.se.ctp.gz --twoway --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0

join.ctp.gz: reads.se.ctp.gz reads.pe.one.ctp.gz reads.pe.two.ctp.gz
	$(CTX) pjoin -m 1M --out $@ reads.se.ctp.gz reads.pe.one.ctp.gz reads.pe.two.ctp.gz

genome.k9.pdf:
	../../scripts/cortex_to_graphviz.pl genome.k9.ctx:0 | dot -Tpdf > genome.k9.pdf

check: $(TGTS)
	$(CTX) check -p reads.se.ctp.gz -p reads.pe.one.ctp.gz -p reads.pe.two.ctp.gz genome.k9.ctx

.PHONY: all plots clean check
