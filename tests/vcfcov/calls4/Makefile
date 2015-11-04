SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

#
# Test variants exactly k-1 bp apart (mix of SNPs, insertions and deletions)
#
# Test files (calls.vcf pluto.fa) made with:
#   python3 fake-spaced-vcf.py ../ref/ref.fa 21 pluto 2 3 > calls.vcf 2> pluto.fa
#

all: test

clean:
	rm -rf calls.cov.vcf* lowmem.cov.vcf* pluto.k$(K).ctx

%.k$(K).ctx: %.fa
	$(MCCORTEX31) build -m 10M -k $(K) --sample $* --seq $< $@

calls.cov.vcf.log: calls.cov.vcf
calls.cov.vcf: calls.vcf $(REF) pluto.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 1 $< pluto.k$(K).ctx >& $@.log

lowmem.cov.vcf: calls.vcf $(REF) pluto.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 1 --low-mem $< pluto.k$(K).ctx >& $@.log

test: calls.cov.vcf lowmem.cov.vcf truth.cov.vcf calls.cov.vcf.log
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	diff -q <($(VCFENTRIES) lowmem.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."
	[[ `grep -o 'max alleles in buffer:.*' calls.cov.vcf.log | grep -o '[0-9][0-9]*'` -lt 3 ]]
	@echo "=> Buffer kept below 3 VCF entries."

view: calls.cov.vcf truth.cov.vcf
	gzip -fcd calls.cov.vcf
	gzip -fcd truth.cov.vcf

.PHONY: all clean view test
