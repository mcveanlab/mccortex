SHELL=/bin/bash -euo pipefail

#
# Test we can place a het-non-ref SNP (triallelic),
# when it's flanks don't match the ref
#

CTXDIR=../../..
MCCORTEX31=$(CTXDIR)/bin/mccortex31
MAKEJOB=$(CTXDIR)/scripts/make-pipeline.pl
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip
VCFCOUNT=$(CTXDIR)/libs/biogrok/vcf-count

K=31
OUTDIR=proj
VCF=$(OUTDIR)/vcfs/bubbles.joint.links.k$(K).vcf.gz
RAWVCF=$(OUTDIR)/k$(K)/bubbles_links/joint.bub.raw.vcf

all: $(VCF) test

ref/ref.fa:
	mkdir -p ref
	( echo ">chr1"; \
	  echo "GGAAGGCGTACTTCGGCATTGTAGCGACCTAAACAAAGGTCACAGATTGGGTTCCTCCGCGCAGGGTCAAGGAACACATA"; ) > $@

sample/sample.fa:
	mkdir -p sample
	( echo ">read1"; \
	  echo "GGAAGGCGTACTTCGGCATTGTAGCGtCCTAAACAAAGGaCACAGATTGGGTTCCT-CGCGCAGGGTCAAGGAACACATA"; \
	  echo ">read2"; \
	  echo "GGAAGGCGTACTTCGGCATTGTAGCGtCCTAAACAAAGGcCACAGATTGGGTTCCT-CGCGCAGGGTCAAGGAACACATA"; ) | tr -d '-' > $@

truth.vcf.gz:
	( DATE=$$(date '+%Y%m%d'); \
	  echo '##fileformat=VCFv4.2'; \
	  echo '##fileDate='"$$DATE"''; \
		echo '##FILTER=<ID=PASS,Description="All filters passed">'; \
		echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'; \
		echo '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT' | tr ' ' '\t'; \
		echo "chr1 40 . T A . PASS . GT" | tr ' ' '\t'; \
	  echo "chr1 40 . T C . PASS . GT" | tr ' ' '\t'; ) > truth.vcf
	$(BGZIP) truth.vcf
	$(BCFTOOLS) index truth.vcf.gz

samples.txt:
	printf "Betty\tsample/sample.fa\n" > $@

job.mk: samples.txt ref/ref.fa sample/sample.fa
	$(MAKEJOB) -r ref/ref.fa $(K) $(OUTDIR) $< >$@

$(RAWVCF) $(VCF): job.mk
	$(MAKE) -f job.mk CTXDIR=$(CTXDIR) KMER_CLEANING_ARGS="--unitigs=0" bub-vcf

test: $(VCF) truth.vcf.gz
	[ `$(VCFCOUNT) $(VCF)` -eq 2 ]
	[ `$(VCFCOUNT) $(RAWVCF)` -eq 4 ]
	[ `$(BCFTOOLS) isec -w 1 -n =2 truth.vcf.gz $(VCF) | $(VCFCOUNT)` -eq 2 ]
	@echo "=> all looks good."

clean:
	rm -rf $(OUTDIR) ref sample samples.txt job.mk truth.vcf*

.PHONY: all clean test
