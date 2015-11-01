SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex63
K=51

TGTS=seq.fa seq.k$(K).ctx sort.k$(K).ctx sort.k$(K).ctx.idx

all: $(TGTS)

clean:
	rm -rf $(TGTS)

seq.fa:
	$(DNACAT) -F -n 100 > $@

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -k $(K) --sample Jimmy --seq $< $@
	$(MCCORTEX) check -q $@

sort.k$(K).ctx: seq.k$(K).ctx
	$(MCCORTEX) sort -o $@ $<
	$(MCCORTEX) check -q $@

sort.k$(K).ctx.idx: sort.k$(K).ctx
	$(MCCORTEX) index --out $@ --block-kmers 11 $<

.PHONY: all clean
