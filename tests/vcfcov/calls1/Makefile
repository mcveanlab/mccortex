SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

all: check

clean:
	rm -rf calls.cov.vcf graph.k$(K).ctx

calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX31) vcfcov -o $@ -r $(REF) calls.vcf graph.k$(K).ctx

graph.k$(K).ctx: sample.1.fa sample.2.fa
	$(MCCORTEX31) build -k $(K) \
	  --sample John --seq sample.1.fa \
	  --sample Jane --seq sample.2.fa \
	  --sample Empty --seq <(echo '') \
	  $@

check: calls.cov.vcf truth.cov.vcf
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)

view: calls.cov.vcf truth.cov.vcf
	zcat -f calls.cov.vcf
	zcat -f truth.cov.vcf

.PHONY: all clean view check
