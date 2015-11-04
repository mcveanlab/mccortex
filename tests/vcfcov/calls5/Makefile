SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries

REF=../ref/ref.fa
K=11

#
# Test that large indels are used and the variants they overlap are not grouped
# with them to stop genotyping of near by variants
#

all: test

alt.fa:
	printf '>alt del:15-35 40T>A\n\
ACTATGGCCAAAGAgCTAGGaTGTTTTTCGGCTCAAGACTCTATCCTGCGGACCGTTCCGCAGGCGTGCCCAGCACCAGGGTCCGTACATTAATACCGTCGCGACTTACTTATTAAGCGTAGGGCACAGCAATATTTCCGCTGGCCCTTACAACCTAGTTTGTCCATAGAGCCATCATAGG\n' > $@

clean:
	rm -rf calls.cov.vcf* lowmem.cov.vcf* wally.k$(K).ctx alt.fa

%.k$(K).ctx: $(REF) alt.fa
	$(MCCORTEX31) build -m 10M -k $(K) --sample $* -1 $(REF) -1 $(REF) -1 alt.fa -1 alt.fa  $@

calls.cov.vcf.log: calls.cov.vcf
calls.cov.vcf: calls.vcf $(REF) wally.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 5 $< wally.k$(K).ctx >& $@.log

lowmem.cov.vcf: calls.vcf $(REF) wally.k$(K).ctx
	$(MCCORTEX31) vcfcov -m 10M -o $@ -r $(REF) --max-nvars 5 --low-mem $< wally.k$(K).ctx >& $@.log

test: calls.cov.vcf lowmem.cov.vcf truth.cov.vcf calls.cov.vcf.log
	diff -q <($(VCFENTRIES) calls.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	diff -q <($(VCFENTRIES) lowmem.cov.vcf) <($(VCFENTRIES) truth.cov.vcf)
	@echo "=> VCF files match."
	[[ `grep -o 'max alleles in buffer:.*' calls.cov.vcf.log | grep -o '[0-9][0-9]*'` -lt 7 ]]
	@echo "=> Buffer kept below 7 VCF entries."

view: calls.cov.vcf truth.cov.vcf
	gzip -fcd calls.cov.vcf
	gzip -fcd truth.cov.vcf

.PHONY: all clean view test
