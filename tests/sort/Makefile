SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/mccortex63
K=51

TGTS=seq.fa seq.k$(K).ctx sort.k$(K).ctx

all: $(TGTS)

clean:
	rm -rf $(TGTS)

seq.fa:
	$(DNACAT) -F -n 100 > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample Jimmy --seq $< $@
	$(CTX) check -q $@

sort.k$(K).ctx: seq.k$(K).ctx
	$(CTX) sort -o $@ $<
	$(CTX) check -q $@

.PHONY: all clean
