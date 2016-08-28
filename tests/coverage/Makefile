SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex31
K=5

TGTS=seq.fa rnd.fa seq.k$(K).ctx coverage.txt

all: $(TGTS)

clean:
	rm -rf $(TGTS)

%.fa:
	$(DNACAT) -F -n 50 > $@

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -q -k $(K) --sample Wallace \
	                             --sample Gromit --seq seq.fa \
	                             --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(MCCORTEX) view -qk $@

coverage.txt: seq.k$(K).ctx rnd.fa
	$(MCCORTEX) coverage -q --seq rnd.fa -1 seq.fa seq.k$(K).ctx > coverage.txt
	cat coverage.txt

.PHONY: all clean
