SHELL=/bin/bash -euo pipefail

#
# Test vcfcov with three groups of overlapping SNPs at positions ref:1,50,199
# and chr1:30. Length of chromosome is ref=200, chr1=100.
# We also test that we don't crash if we encounter a contig that was not defined
# in the header.
#

K=21
CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa

all: check

clean:
	rm -rf calls.cov.vcf lowmem.cov.vcf graph.k$(K).ctx

calls.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX) vcfcov -m 10M -o $@ -r $(REF) --high-mem calls.vcf graph.k$(K).ctx >& $@.log

lowmem.cov.vcf: $(REF) calls.vcf graph.k$(K).ctx
	$(MCCORTEX) vcfcov -m 10M -o $@ -r $(REF) --low-mem calls.vcf graph.k$(K).ctx >& $@.log

graph.k$(K).ctx: john.fa jane.fa
	$(MCCORTEX) build -m 10M -k $(K) \
	  --sample John --seq john.fa \
	  --sample Jane --seq jane.fa \
	  --sample Empty --seq <(echo '') \
	  $@ >& $@.log

check: calls.cov.vcf lowmem.cov.vcf truth.cov.vcf
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	diff -q <($(VCFENTRIES) lowmem.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."

view: calls.cov.vcf truth.cov.vcf
	gzip -fcd calls.cov.vcf
	gzip -fcd truth.cov.vcf

.PHONY: all clean view check
