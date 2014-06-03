
CTX=../../bin/ctx31
PATHS=reads.se.ctp reads.pe.one.ctp reads.pe.two.ctp
TGTS=genome.fa read.1.fa read.2.fa genome.k9.ctx $(PATHS)

# non-default target: genome.k9.pdf

all: $(TGTS)

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
	$(CTX) build -m 1M -k 9 --sample Genome --seq genome.fa --sample Gen2 --seq genome.fa genome.k9.ctx

reads.se.ctp: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths --seq read.1.fa --seq read.2.fa -o $@ genome.k9.ctx:0

reads.pe.one.ctp: reads.se.ctp genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths -p reads.se.ctp -C --oneway --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0

reads.pe.two.ctp: reads.se.ctp genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --print-paths -p reads.se.ctp --twoway --seq2 read.1.fa:read.2.fa -o $@ genome.k9.ctx:0

genome.k9.pdf:
	../../scripts/cortex_to_graphviz.pl genome.k9.ctx:0 | dot -Tpdf > genome.k9.pdf

.PHONY: all plots clean
