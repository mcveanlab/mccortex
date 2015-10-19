SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

#
# Test vcfcov with too many overlapping variants
#

all: check

clean:
	rm -rf calls.cov.vcf graph.k$(K).ctx

calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 4 calls.vcf graph.k$(K).ctx

graph.k$(K).ctx: sample.fa
	$(MCCORTEX31) build -m 10M -k $(K) --sample John --seq sample.fa $@

check: calls.cov.vcf truth.cov.vcf
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)

view: calls.cov.vcf truth.cov.vcf
	zcat -f calls.cov.vcf
	zcat -f truth.cov.vcf

.PHONY: all clean view check
