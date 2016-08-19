SHELL:=/bin/bash -euo pipefail

#
# Isaac Turner 2014-10-08
#
# In this test we generate a random 1001 base sequence, add 9bp unique kmers
# either side and build graph and link files. Links are generated using the
# entire sequence. Then, using k=9, we reassemble contigs from the graph+links,
# remove duplicates, trim off 9bp added each side and compare to original.
# It should match perfectly.
#

#
# Inferring edges breaks lossless recovery because forks can be induced
# before we have a chance to pick up links:
#
# start__/ ... \__end
#        \     /
#

K=9
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

# 4**7 = 16,384
GENOME=10000
TGTS=seq.fa seq.uniq.fa seq.uniq.k$(K).ctx seq.uniq.k$(K).ctp.gz \
     seq.regen.fa seq.regen.txt

all: $(TGTS) check

seq.fa:
	$(DNACAT) -n $(GENOME) -F -M <(echo seq) > $@

seq.uniq.fa: seq.fa
	$(MCCORTEX) uniqkmers -q -m 10M -k $$[$(K)-2] -F $< 0 > $@

seq.uniq.k$(K).ctx: seq.uniq.fa
	$(MCCORTEX) build -q -m 10M -k $(K) --sample SeqUniq --seq $< $@
	# Inferring edges breaks lossless recovery, see not above
	# $(MCCORTEX) inferedges $@

seq.uniq.k$(K).ctp.gz: seq.uniq.k$(K).ctx seq.uniq.fa
	$(MCCORTEX) thread -q --seq seq.uniq.fa --out $@ seq.uniq.k$(K).ctx

seq.regen.fa: seq.uniq.k$(K).ctx seq.uniq.k$(K).ctp.gz
	$(MCCORTEX) contigs -q --no-missing-check -G $(GENOME) -o - -p seq.uniq.k$(K).ctp.gz seq.uniq.k$(K).ctx | \
	$(MCCORTEX) rmsubstr -q -n $$[1000+2*$(GENOME)] -k $(K) - > $@

seq.regen.txt: seq.regen.fa
	$(DNACAT) -P $< | awk '{print substr($$0,$(K)-2+1,length($$0)-($(K)-2)*2)}'> $@

check: seq.fa seq.regen.txt seq.uniq.k$(K).ctx
	@echo Check kmer count distribution...
	$(MCCORTEX) view -q -k seq.uniq.k$(K).ctx | awk '{print $$2}' | sort -n | uniq -c
	@echo; echo Check only one contig assembled...
	[ `$(DNACAT) -s seq.regen.txt | grep 'Number of reads' | grep -oE '[0-9]+'` -eq 1 ]
	@echo; echo Checking if regenerated file matches original...
	diff -q <($(DNACAT) -r -k -P seq.fa) <($(DNACAT) -r -k -P seq.regen.txt)
	@echo

clean:
	rm -rf $(TGTS)

.PHONY: all clean
