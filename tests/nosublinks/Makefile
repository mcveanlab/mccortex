SHELL:=/bin/bash -euo pipefail

CTXDIR=../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex31
K=5

SEQS=seq.a.fa seq.b.fa seq.c.fa
LINKS=$(SEQS:.fa=.k$(K).ctp.gz)
TGTS=$(SEQS) seq.k$(K).ctx $(LINKS)

all: $(TGTS)

clean:
	rm -rf $(TGTS)

# a: CGTGCCCTACAGTCACTCCCG
# b: tcTGCCCTACAGTCACTCCac
# c: tgctcgCTACAGTCAgacacg

seq.a.fa:
	echo CGTGCCCTACAGTCACTCCCG > $@

seq.b.fa:
	echo tcTGCCCTACAGTCACTCCac > $@

seq.c.fa:
	echo tgctcgCTACAGTCAgacacg > $@

seq.k$(K).ctx: $(SEQS)
	$(MCCORTEX) build -q -m 1M -k $(K) --sample Plato -1 seq.a.fa -1 seq.b.fa -1 seq.c.fa $@
	$(MCCORTEX) inferedges -q $@

%.k$(K).ctp.gz: seq.k$(K).ctx $(SEQS)
	$(MCCORTEX) thread -q -o $*.k5.ctp.gz -1 $*.fa --no-sublinks seq.k$(K).ctx

.PHONY: all clean
