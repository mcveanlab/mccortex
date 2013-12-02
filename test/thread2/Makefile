
SHELL:=/bin/bash
K=9
CTX=../../bin/ctx31
REVCMP=../../libs/seq_file/bin/revcmp

GRAPHS=$(shell echo genome.k$(K).{ctx,dot})
PLOTS=$(shell echo genome.k$(K).{pdf,kmers.pdf})
READSSE=$(shell echo read{1,2}{a,b}.fa)
READSPE=$(shell echo read{3,4}.{1,2}.fa)
PATHS=$(shell echo reads.{se,pe}.k$(K).ctp)

KEEP=genome.fa $(GRAPHS) $(READSSE) $(READSPE) $(PATHS)

# ctx thread arguments
SEQSE=$(shell printf " --seq %s" $(READSSE))
SEQPE=$(shell echo --seq2 read3.{1,2}.fa --seq2 read4.{1,2}.fa)

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
	echo GCTAAGGAA > read3.1.fa; echo GGTACAGTG | $(REVCMP) > read3.2.fa
	echo GCTCAGGAA > read4.1.fa; echo GGTACTGTG | $(REVCMP) > read4.2.fa

genome.k$(K).ctx:
	$(CTX) build -k $(K) --sample MssrGenome --seq genome.fa genome.k$(K).ctx

genome.k$(K).dot:
	$(CTX) supernodes --graphviz genome.k$(K).ctx > genome.k$(K).dot

genome.k$(K).pdf:
	dot -Tpdf genome.k$(K).dot > genome.k$(K).pdf

genome.k$(K).kmers.pdf:
	../../scripts/cortex_to_graphviz.pl genome.k$(K).ctx | dot -Tpdf > genome.k$(K).kmers.pdf

reads.se.k$(K).ctp:
	$(CTX) thread -t 6 --col 0 $(SEQSE) reads.se.k$(K).ctp genome.k$(K).ctx
	$(CTX) pview reads.se.k$(K).ctp

reads.pe.k$(K).ctp:
	$(CTX) thread -t 6 --col 0 $(SEQSE) $(SEQPE) reads.pe.k$(K).ctp genome.k$(K).ctx
	$(CTX) pview reads.pe.k$(K).ctp

.PHONY: all clean plots
