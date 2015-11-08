SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

#
# Test vcfcov with too many overlapping variants
#
# Test VCF of a SNP per base generated with:
#   ./fake-vcf.py ../ref/ref.fa > calls.vcf
#

all: test

clean:
	rm -rf calls.cov.vcf* lowmem.cov.vcf* graph.k$(K).ctx

calls.cov.vcf.log: calls.cov.vcf
calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 4 calls.vcf graph.k$(K).ctx >& $@.log

lowmem.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 4 --low-mem calls.vcf graph.k$(K).ctx >& $@.log

graph.k$(K).ctx: sample.fa
	$(MCCORTEX31) build -m 10M -k $(K) --sample John --seq sample.fa $@ >& $@.log

test: calls.cov.vcf lowmem.cov.vcf truth.cov.vcf calls.cov.vcf.log
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	diff -q <($(VCFENTRIES) lowmem.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."
	[[ `grep -o 'max alleles in buffer:.*' calls.cov.vcf.log | grep -o '[0-9][0-9]*'` -lt 50 ]]
	@echo "=> Buffer kept below 50 VCF entries."

view: calls.cov.vcf truth.cov.vcf
	gzip -fcd calls.cov.vcf
	gzip -fcd truth.cov.vcf

.PHONY: all clean view test
