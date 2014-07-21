
CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
BIOINF=$(CTXDIR)/libs/bioinf-perl
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
CTX=$(CTXDIR)/bin/ctx31

K=9

all: test

seq.fa:
	$(DNACAT) -F -n 1001 > seq.fa

seq.k$(K).ctx: seq.fa
	$(CTX) build -m 1M -k $(K) --sample MrsSample --seq seq.fa $@

seq.k$(K).ctp.gz: seq.k$(K).ctx seq.fa
	$(CTX) thread -m 1M --seq seq.fa --out $@ seq.k$(K).ctx

contigs.fa: seq.k$(K).ctx seq.k$(K).ctp.gz
	$(CTX) contigs --no-reseed --out $@ -p seq.k$(K).ctp.gz seq.k$(K).ctx

contigs.rmdup.fa: contigs.fa
	$(CTX) rmsubstr $< > $@

seq.k$(K).pdf: seq.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: seq.k$(K).pdf

test: seq.fa contigs.fa contigs.rmdup.fa
	$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.fa seq.fa
	$(BIOINF)/sim_mutations/sim_substrings.pl $(K) 0.1 contigs.rmdup.fa seq.fa

clean:
	rm -rf seq.fa seq.k$(K).ctx seq.k$(K).ctp.gz contigs.fa contigs.rmdup.fa seq.k$(K).pdf

.PHONY: all clean test plots
