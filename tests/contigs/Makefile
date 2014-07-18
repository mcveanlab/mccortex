
CTXDIR=../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dna
BIOINF=$(CTXDIR)/libs/bioinf-perl
CTX=$(CTXDIR)/bin/ctx31

K=9

all: test

seq.fa:
	$(DNACAT) -n 1001 -f > seq.fa

seq.k9.ctx: seq.fa
	$(CTX) build -m 1M -k $(K) --sample MrsSample --seq seq.fa $@

seq.k9.ctp.gz: seq.k9.ctx seq.fa
	$(CTX) thread -m 1M --seq seq.fa --out $@ seq.k$(K).ctx

contigs.fa: seq.k9.ctx seq.k9.ctp.gz
	$(CTX) contigs --ncontigs 10 --out $@ -p seq.k9.ctp.gz seq.k$(K).ctx

test: seq.fa contigs.fa
	$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.fa seq.fa

clean:
	rm -rf seq.fa seq.k9.ctx seq.k$(K).ctp.gz contigs.fa

.PHONY: all clean test
