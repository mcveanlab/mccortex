SHELL:=/bin/bash

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx31
K=5

TGTS=seq.fa rnd.fa seq.k$(K).ctx coverage.txt

all: $(TGTS)

clean:
	rm -rf $(TGTS)

%.fa:
	$(DNACAT) -F -n 50 > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample Wallace \
	                     --sample Gromit --seq seq.fa \
	                     --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(CTX) view --kmers $@

coverage.txt: seq.k$(K).ctx rnd.fa
	$(CTX) coverage --seq rnd.fa -1 seq.fa seq.k$(K).ctx > coverage.txt
	cat coverage.txt

.PHONY: all clean
