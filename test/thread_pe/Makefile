
CTX=../../bin/ctx31
PATHS=reads.pe.one.ctp reads.pe.two.ctp
TGTS=genome.fa read.1.fa read.2.fa genome.k9.ctx $(PATHS)

# non-default target: genome.k9.pdf

all: $(TGTS)

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

reads.pe.one.ctp: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --printcontigs --oneway --col 0 --seq2 read.1.fa read.2.fa --col 1 --seq genome.fa $@ genome.k9.ctx

reads.pe.two.ctp: genome.k9.ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --printcontigs --twoway --col 0 --seq2 read.1.fa read.2.fa --col 1 --seq genome.fa $@ genome.k9.ctx

genome.k9.pdf:
	../../scripts/cortex_to_graphviz.pl genome.k9.ctx | dot -Tpdf > genome.k9.pdf
