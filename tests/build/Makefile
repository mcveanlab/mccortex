SHELL:=/bin/bash

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx31
K=9

all: seq.k$(K).ctx

clean:
	rm -rf seq.fa seq.k*.ctx test.ctx A:2.fa

seq.fa:
	$(DNACAT) -F -n 50 > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample Wallace \
	                     --sample Gromit --seq seq.fa \
	                     --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(CTX) check -q $@
	$(CTX) view $@

test.ctx:
	$(SEQRND) 20 > A:2.fa
	$(CTX) build -k $(K) --sample Colon -1 A:2.fa $@
	$(CTX) check -q $@
	rm A:2.fa

.PHONY: all clean
