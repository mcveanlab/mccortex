SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=21

#
# Test two samples on several overlapping variants
#
# max vars=4, max-allele=20
#
# large ref insertion to ignore
# 1 that can be done
# block of 5 that can't be done
# 2 that can be done
# <gap>
# 3 snps
# <gap>
# large ref deletion to ignore
# 2 that can be done
# block  of 8 that can't
# 1 that can be done
#

# ref:80C>T num kmers => 85-80=5, k=21, 5+(k-5)*2 = 37, mean cov = (5*2+35*1)/37 = 1.21

all: test

clean:
	rm -rf calls.cov.vcf* lowmem.cov.vcf* laurel.k$(K).ctx hardy.k$(K).ctx

%.k$(K).ctx: %.fa
	$(MCCORTEX31) build -m 10M -k $(K) --sample $* --seq $< $@

calls.cov.vcf: calls.vcf $(REF) laurel.k$(K).ctx hardy.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 4 --max-var-len 20 $< laurel.k$(K).ctx hardy.k$(K).ctx >& $@.log

lowmem.cov.vcf: calls.vcf $(REF) laurel.k$(K).ctx hardy.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 4 --max-var-len 20 --low-mem $< laurel.k$(K).ctx hardy.k$(K).ctx >& $@.log

test: calls.cov.vcf lowmem.cov.vcf truth.cov.vcf
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	diff -q <($(VCFENTRIES) lowmem.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."
	[[ `grep -o 'max alleles in buffer:.*' calls.cov.vcf.log | grep -o '[0-9][0-9]*'` -lt 9 ]]
	@echo "=> Buffer kept below 9 VCF entries."

view: calls.cov.vcf truth.cov.vcf
	gzip -cdf calls.cov.vcf
	gzip -cdf truth.cov.vcf

.PHONY: all clean view test
