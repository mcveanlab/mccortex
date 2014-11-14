SHELL:=/bin/bash -euo pipefail

#
# Isaac Turner 2014-10-08
#
# In this test we generate a random 1001 base sequence, add 9bp unique kmers
# either side and build graph and path files. Then, using k=9, we reassemble
# contigs from the graph+path, remove duplicates, trim off 9bp added each side
# and compare to original.
#

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx31

# 4**7 = 16,384
K=9
GENOME=10000
TGTS=seq.fa seq.uniq.fa seq.uniq.k$(K).ctx seq.uniq.k$(K).ctp.gz seq.regen.fa

all: $(TGTS) check

seq.fa:
	$(DNACAT) -n $(GENOME) -F -M <(echo seq) > $@

seq.uniq.fa: seq.fa
	$(CTX) uniqkmers -m 10M -k $$[$(K)-2] -F $< 0 > $@

seq.uniq.k$(K).ctx: seq.uniq.fa
	$(CTX) build -m 10M -k $(K) --sample SeqUniq --seq $< $@
	# $(CTX) inferedges $@

seq.uniq.k$(K).ctp.gz: seq.uniq.k$(K).ctx seq.uniq.fa
	$(CTX) thread --seq seq.uniq.fa --out $@ seq.uniq.k$(K).ctx

seq.regen.fa: seq.uniq.k$(K).ctx seq.uniq.k$(K).ctp.gz
	$(CTX) contigs --no-missing-check -G $(GENOME) -o - -p seq.uniq.k$(K).ctp.gz seq.uniq.k$(K).ctx | \
	$(CTX) rmsubstr -n $$[1000+2*$(GENOME)] -k $(K) - | \
	$(DNACAT) -P - | awk '{print substr($$0,$(K)-2+1,length($$0)-($(K)-2)*2)}'> $@

check: seq.fa seq.regen.fa seq.uniq.k$(K).ctx
	@echo Check kmer count distribution...
	$(CTX) view -q -k seq.uniq.k$(K).ctx | awk '{print $$2}' | sort -n | uniq -c
	@echo; echo Check only one contig assembled...
	[ `$(DNACAT) -s seq.regen.fa | grep 'Number of reads' | grep -oE '[0-9]+'` -eq 1 ]
	@echo; echo Checking if regenerated file matches original...
	diff -q <($(DNACAT) -r -k -P seq.fa) <($(DNACAT) -r -k -P seq.regen.fa)
	@echo

clean:
	rm -rf $(TGTS)

.PHONY: all clean
