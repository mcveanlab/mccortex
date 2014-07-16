
CTXDIR=../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dna
BIOINF=$(CTXDIR)/libs/bioinf-perl
CTX=$(CTXDIR)/bin/ctx31

K=9

all: test

seq.fa:
	$(DNACAT) -n 1001 -f > seq.fa

seq.k9.ctx: seq.fa
	$(CTX) build -k $(K) --sample MrsSample --seq seq.fa seq.k$(K).ctx

test: seq.fa seq.k$(K).ctx
	 $(CTX) contigs --ncontigs 10 --out - seq.k$(K).ctx | \
	 $(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 - seq.fa

clean:
	rm -rf seq.fa seq.k9.ctx

.PHONY: all clean test
