SHELL=/bin/bash -euo pipefail

#
# Test vcfcov with an empty VCF
# Works with any kmer size (K)
#

K=21
CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)
VCFNALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles

REF=../ref/ref.fa

VCFS=calls.cov.k$(K).vcf lowmem.cov.k$(K).vcf
GRAPH=graph.k$(K).ctx
LOGS=$(VCFS:=.log) $(GRAPH:=.log)

all: test

clean:
	rm -rf $(VCFS) $(LOGS) $(GRAPH)

calls.cov.k$(K).vcf: $(REF) calls.vcf $(GRAPH)
	$(MCCORTEX) vcfcov -m 10M -o $@ -r $(REF) --high-mem calls.vcf $(GRAPH) >& $@.log

lowmem.cov.k$(K).vcf: $(REF) calls.vcf $(GRAPH)
	$(MCCORTEX) vcfcov -m 10M -o $@ -r $(REF) --low-mem calls.vcf graph.k$(K).ctx >& $@.log

graph.k$(K).ctx:
	$(MCCORTEX) build -k $(K) \
	  --sample John --seq <(echo '') \
	  --sample Jane --seq <(echo '') \
	  $@ >& $@.log

test: $(VCFS)
	[[ `$(VCFNALLELES) calls.cov.k$(K).vcf` -eq 0 ]]
	[[ `$(VCFNALLELES) lowmem.cov.k$(K).vcf` -eq 0 ]]
	@echo "=> Empty VCF works."

view: calls.cov.k$(K).vcf
	gzip -fcd $<

.PHONY: all clean view test
