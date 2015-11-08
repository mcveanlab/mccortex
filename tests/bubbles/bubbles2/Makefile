SHELL=/bin/bash -euo pipefail

#
# Test we call all bubbles when filtering ref bubbles
#


CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFNALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
VCFSORT=$(CTXDIR)/libs/biogrok/vcf-sort
BWA=$(CTXDIR)/libs/bwa/bwa

K=31

SAMPLES=itchy scratchy ref
FASTAS=$(SAMPLES:=.fa)
GRAPHS=$(SAMPLES:=.k$(K).ctx)

all: bubbles.txt test

itchy.fa: scratchy.fa
scratchy.fa: ref.fa
ref.fa:
	python make-exhaustive-alleles.py itchy.fa scratchy.fa ref.fa

%.k$(K).ctx: %.fa
	$(MCCORTEX31) build -k $(K) --sample "$*" --seq $< $@ >& $@.log

bubbles.txt.log: bubbles.txt
bubbles.txt: $(GRAPHS)
	$(MCCORTEX31) bubbles -o $@.gz --haploid 2 $(GRAPHS) >& $@.log
	gzip -fd $@.gz

# expect 60 bubbles called = (2**4-1)*2*2
# expect 64 haploid bubbles = (2**4)*2*2 (i.e. bubbles in ref ignored)
test: bubbles.txt.log
	grep -q 'Bubble Caller called 60 bubbles' $<
	grep -q 'Haploid bubbles dropped: 64' $<
	grep -q 'Serial bubbles dropped: 0' $<
	@echo "=> all looks good."

clean:
	rm -rf $(GRAPHS) bubbles.txt *.log

.PHONY: all clean test
