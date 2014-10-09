SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx31
CONTIGSTATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl
K=9

TGTS=seq.fa seq.k$(K).ctx sort.k$(K).ctx sort.k$(K).ctx.idx

all: $(TGTS) test_assemble

clean:
	rm -rf $(TGTS)

seq.fa:
	$(DNACAT) -F -n 50 > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample Wallace \
	                     --sample Gromit --seq seq.fa \
	                     --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(CTX) check -q $@
	$(CTX) view $@
	$(CTX) view -k 1,3,5:$@:2,1,0

sort.k$(K).ctx: seq.k$(K).ctx
	cp $< $@
	$(CTX) view -k $<
	$(CTX) sort $@
	$(CTX) check -q $@
	$(CTX) view -k $@

%.ctx.idx: %.ctx
	$(CTX) index $< > $@
	cat $@

test.ctx:
	$(SEQRND) 20 > A:2.fa
	$(CTX) build -k $(K) --sample Goofy -1 A:2.fa $@
	$(CTX) check -q $@
	rm A:2.fa

test_assemble:
	$(CTX) contigs -q -o - seq.k9.ctx:1 | \
	  $(CTX) rmsubstr -q -n 1M - | \
	  $(CONTIGSTATS) -

.PHONY: all clean test_assemble
