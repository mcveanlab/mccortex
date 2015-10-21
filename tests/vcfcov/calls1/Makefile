SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

#
# Test vcfcov with three groups of overlapping SNPs at positions 1, 50 and 199
# on a reference 200bp long
#

all: check

clean:
	rm -rf calls.cov.vcf graph.k$(K).ctx

calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) calls.vcf graph.k$(K).ctx

graph.k$(K).ctx: john.fa jane.fa
	$(MCCORTEX31) build -m 10M -k $(K) \
	  --sample John --seq john.fa \
	  --sample Jane --seq jane.fa \
	  --sample Empty --seq <(echo '') \
	  $@

check: calls.cov.vcf truth.cov.vcf
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."

view: calls.cov.vcf truth.cov.vcf
	zcat -f calls.cov.vcf
	zcat -f truth.cov.vcf

.PHONY: all clean view check
