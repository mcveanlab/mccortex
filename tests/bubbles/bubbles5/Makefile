SHELL=/bin/bash -euo pipefail

CTXDIR=../../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
READSIM=$(CTXDIR)/libs/readsim/readsim
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries
BGZIP=$(CTXDIR)/libs/htslib/bgzip
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools

MKPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl

BRKVCF=proj/vcfs/breakpoints.joint.links.k31.geno.vcf.gz
BUBVCF=proj/vcfs/bubbles.joint.links.k31.geno.vcf.gz
CALLS=$(BRKVCF) $(BUBVCF)

DEPTH=50
ERRRATE=0.01
READLEN=50

all: test

$(CALLS): job.mk ref/ref.fa reads.fa.gz
	$(MAKE) -f job.mk CTXDIR=$(CTXDIR) MEM=50M NKMERS=1M vcfs

ref/ref.fa:
	mkdir -p ref
	( echo '>ref'; echo \
'CCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCAtCCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTGAC' | tr -d '-' ) > $@

sample.fa:
	( echo '>sample'; echo \
'CCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCAACCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTGAC' | tr -d '-' ) > $@

reads.fa.gz: sample.fa ref/ref.fa
	$(READSIM) -d $(DEPTH) -l $(READLEN) -s -e $(ERRRATE) -r <($(DNACAT) sample.fa ref/ref.fa) reads

samples.txt:
	echo "MrBond reads.fa.gz" > $@

job.mk: samples.txt ref/ref.fa
	$(MKPIPELINE) -r ref/ref.fa 31 proj $< > $@

view: | calls
	$(VCFENTRIES) $(BRKVCF)
	$(VCFENTRIES) $(BUBVCF)

truth.vcf.gz:
	$(BGZIP) -c truth.vcf > $@

truth.vcf.gz.csi: truth.vcf.gz
	$(BCFTOOLS) index $<

test: truth.vcf.gz truth.vcf.gz.csi $(CALLS)
	$(BCFTOOLS) isec truth.vcf.gz $(BUBVCF) -p bub_isec
	diff -q <($(VCFENTRIES) truth.vcf.gz) <($(VCFENTRIES) bub_isec/0002.vcf)
	@echo '=> Genotyped VCFs match!'

clean:
	rm -rf samples.txt job.mk ref sample.fa reads.fa.gz proj truth.vcf.gz* bub_isec

.PHONY: all clean calls view test
