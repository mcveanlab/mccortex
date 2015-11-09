SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFNALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles

REF=../ref/ref.fa
K=21

#
# Test vcfcov with an empty VCF
#

all: test

clean:
	rm -rf calls.cov.vcf lowmem.cov.vcf graph.k$(K).ctx

calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) calls.vcf graph.k$(K).ctx >& $@.log

lowmem.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --low-mem calls.vcf graph.k$(K).ctx >& $@.log

graph.k$(K).ctx:
	$(MCCORTEX31) build -k $(K) \
	  --sample John --seq <(echo '') \
	  --sample Jane --seq <(echo '') \
	  $@ >& $@.log

test: calls.cov.vcf lowmem.cov.vcf
	[[ `$(VCFNALLELES) calls.cov.vcf` -eq 0 ]]
	[[ `$(VCFNALLELES) lowmem.cov.vcf` -eq 0 ]]
	@echo "=> Empty VCF works."

view: calls.cov.vcf
	gzip -fcd calls.cov.vcf

.PHONY: all clean view test
