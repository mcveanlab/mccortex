
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
	$(CTX) contigs --no-reseed --out $@ -p seq.k9.ctp.gz seq.k$(K).ctx

contigs.rmdup.fa: contigs.fa
	$(CTX) rmsubstr $< > $@

test: seq.fa contigs.fa contigs.rmdup.fa
	$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.fa seq.fa
	$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.rmdup.fa seq.fa

clean:
	rm -rf seq.fa seq.k9.ctx seq.k$(K).ctp.gz contigs.fa contigs.rmdup.fa

.PHONY: all clean test
