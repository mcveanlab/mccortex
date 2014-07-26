
SHELL:=/bin/bash

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
K=9

GRAPHS=$(shell echo genome.k$(K).{ctx,dot})
PLOTS=$(shell echo genome.k$(K).{pdf,kmers.pdf})
READSSE=$(shell echo read{1,2}{a,b}.fa)
READSPE=$(shell echo read{3,4}.{1,2}.fa)
PATHS=$(shell echo reads.{se,pe}.k$(K).ctp.gz)

KEEP=genome.fa $(GRAPHS) $(READSSE) $(READSPE) $(PATHS)

# ctx thread arguments
SEQSE=$(shell printf " --seq %s" $(READSSE))
SEQPE=--seq2 read3.1.fa:read3.2.fa --seq2 read4.1.fa:read4.2.fa

all: $(KEEP)

clean:
	rm -rf $(KEEP) $(PLOTS)
	rm -rf gap_sizes.*.csv

plots: $(PLOTS)

# Sample random genome
genome.fa:
	echo GCTAAGGAATACATATACGGTACAGTG > genome.fa
	echo GCTCAGGAATACAGATACGGTACTGTG >> genome.fa

read%.fa:
	echo GCTAAGGAATACAT > read1a.fa
	echo TATACGGTACAGTG > read1b.fa
	echo GCTCAGGAATACAG > read2a.fa
	echo GATACGGTACTGTG > read2b.fa
	set -euo pipefail; echo GCTAAGGAA > read3.1.fa; echo GGTACAGTG | $(DNACAT) --revcmp - > read3.2.fa
	set -euo pipefail; echo GCTCAGGAA > read4.1.fa; echo GGTACTGTG | $(DNACAT) --revcmp - > read4.2.fa

genome.k$(K).ctx:
	$(CTX) build -m 1M -k $(K) --sample MssrGenome --seq genome.fa genome.k$(K).ctx

genome.k$(K).dot:
	$(CTX) supernodes --graphviz genome.k$(K).ctx > genome.k$(K).dot

genome.k$(K).pdf:
	dot -Tpdf genome.k$(K).dot > genome.k$(K).pdf

genome.k$(K).kmers.pdf:
	set -euo pipefail; $(CTX2DOT) genome.k$(K).ctx | dot -Tpdf > genome.k$(K).kmers.pdf

reads.se.k$(K).ctp.gz:
	$(CTX) thread -m 1M -t 6 $(SEQSE) -o reads.se.k$(K).ctp.gz genome.k$(K).ctx
	zcat reads.se.k$(K).ctp.gz

reads.pe.k$(K).ctp.gz:
	$(CTX) thread -m 1M -t 6 $(SEQSE) $(SEQPE) -o reads.pe.k$(K).ctp.gz genome.k$(K).ctx
	zcat reads.pe.k$(K).ctp.gz

.PHONY: all clean plots
