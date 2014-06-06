
CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
SEQRND=$(CTXPATH)/libs/seq_file/bin/seqrnd
K=7
GENOME=200
KEEP=$(shell echo genome.fa genome.k$(K).{ctx,dot,pdf} genome.k$(K).perl.pdf)

all: $(KEEP)

clean:
	rm -rf $(KEEP)

# Sample random genome
genome.fa:
	$(SEQRND) $(GENOME) > genome.fa

genome.k$(K).ctx:
	$(CTX) build -m 1M -k $(K) --sample MssrGenome --seq genome.fa genome.k$(K).ctx

genome.k$(K).dot:
	$(CTX) supernodes -m 1M --graphviz --points genome.k$(K).ctx > genome.k$(K).dot

genome.k$(K).pdf:
	dot -Tpdf genome.k$(K).dot > genome.k$(K).pdf

# Generate a dot file using perl
genome.k$(K).perl.pdf:
	../../scripts/cortex_to_graphviz.pl --simplify genome.k$(K).ctx | dot -Tpdf > $@

.PHONY: all clean
