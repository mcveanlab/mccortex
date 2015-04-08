SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

PATHS=reads.pe.one.ctp reads.pe.two.ctp
TGTS=genome.fa read.1.fa read.2.fa genome.k9.ctx $(PATHS)
KMER=9
# non-default target: genome.k9.pdf

all: $(TGTS)

clean:
	rm -rf $(TGTS) gap_sizes.*.csv mp_sizes.*.csv genome.k$(KMER).pdf

plots: genome.k$(KMER).pdf

genome.fa:
	echo gCATCAGTGGCCttggcgactcgc > genome.fa
	echo TCATCAGTGGCCATGACGCTAACT >> genome.fa

read.1.fa:
	echo TCATCAGTGG > read.1.fa

read.2.fa:
	# echo ACGCTAACT > read.2.fa # Actually revcmp read2
	echo AGTTAGCGT > read.2.fa

genome.k$(KMER).ctx: genome.fa
	$(CTX) build -m 1M -k $(KMER) --sample Genome --seq genome.fa genome.k$(KMER).ctx

reads.pe.one.ctp: genome.k$(KMER).ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --one-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(KMER).ctx

reads.pe.two.ctp: genome.k$(KMER).ctx read.1.fa read.2.fa
	$(CTX) thread -m 1M --print-contigs --two-way --seq2 read.1.fa:read.2.fa -o $@ genome.k$(KMER).ctx

genome.k$(KMER).pdf:
	$(CTX2DOT) genome.k$(KMER).ctx | dot -Tpdf > genome.k$(KMER).pdf

.PHONY: all clean plots
