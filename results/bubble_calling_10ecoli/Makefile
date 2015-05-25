SHELL=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
READSIM=$(CTXDIR)/libs/readsim/readsim
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
SIMTRUTH=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_vcf.pl
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip
VCFALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles

K=21
HAPDEPTH=50
READLEN=150
ERRRATE=0.01
MEM=2G
OUTDIR=proj
REF=../data/ecoli/ecoli.fa

MASKS=$(shell echo genomes/mask{0..9}.fa)
CHROMS_ALN=$(subst mask,genome,$(MASKS))
READS=$(shell echo reads/chrom{0..9}.$(HAPDEPTH)X.1.fa.gz reads/chrom{0..9}.$(HAPDEPTH)X.2.fa.gz)
DIRS=reads genomes

TRUTH_FILES=truth.k$(K).norm.vcf.gz truth.k$(K).norm.vcf.gz.csi

# Mark all dependencies as secondary
# It means don't re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:
# Delete files if their recipe fails
.DELETE_ON_ERROR:
# Remove implicit rules for certain suffixes
.SUFFIXES:

all: run-mccortex $(TRUTH_FILES)

# ecoli 4.6Mb, 1SNP per 100 => 46000 SNPs
# Generate a 10 ecoli genomes
$(MASKS) $(CHROMS_ALN): genomes/about.txt
genomes/about.txt: $(REF) | $(DIRS)
	$(SIMMUT) --snps 46000 --indels 0 --invs 0 genomes 10 $<
	echo "10 Ecoli genomes generated from ref $(REF)" > $@

$(CHROMS): $(CHROMS_ALN)

# Remove deletion marks (-) and convert to uppercase
genomes/chrom%.fa: genomes/genome%.fa
	cat $< | tr -d '-' | $(DNACAT) -u -F - > $@

truth.k$(K).vcf: $(CHROMS) $(MASKS)
	files=$$(for i in {0..9}; do echo genomes/genome$$i.fa genomes/mask$$i.fa; done); \
	$(SIMTRUTH) $(REF) $$files > $@

truth.k$(K).norm.vcf: truth.k$(K).vcf $(REF)
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref $(REF) --multiallelics +both $< > $@

%.vcf.gz: %.vcf
	$(BGZIP) -f $<

%.vcf.gz.csi: %.vcf.gz
	$(BCFTOOLS) index -f $<

# Simulate PE reads of each chrom each 50X
reads/chrom%.$(HAPDEPTH)X.1.fa.gz reads/chrom%.$(HAPDEPTH)X.2.fa.gz: genomes/chrom%.fa
	mkdir -p reads
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom$*.$(HAPDEPTH)X

samples.txt:
	for i in {0..9}; do \
		echo "Ecoli$$i . reads/chrom$$i.$(HAPDEPTH)X.1.fa.gz:reads/chrom$$i.$(HAPDEPTH)X.2.fa.gz"; \
  done > $@

task.k$(K).mk: samples.txt
	$(CTXPIPELINE) -r $(REF) $(K) $(OUTDIR) $< > $@

run-mccortex: task.k$(K).mk $(TRUTH_FILES) $(READS)
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) bubblevcf
	$(BCFTOOLS) isec truth.k$(K).norm.vcf.gz proj/vcfs/bubbles.k$(K).vcf.gz -p truthisec
	@A=`$(VCFALLELES) truthisec/0000.vcf`; \
	 B=`$(VCFALLELES) truthisec/0001.vcf`; \
	 C=`$(VCFALLELES) truthisec/0002.vcf`; \
	MC=`$(VCFALLELES) proj/vcfs/bubbles.k$(K).vcf.gz`; \
	NT=`$(VCFALLELES) truth.k$(K).norm.vcf.gz`; \
	awk 'BEGIN{printf("McCortex Missed: %4d / %4d (%5.2f%%)\n",'$$A','$$NT',100*'$$A'/'$$NT')}'; \
	awk 'BEGIN{printf("McCortex FP:     %4d / %4d (%5.2f%%)\n",'$$B','$$MC',100*'$$B'/'$$MC')}'; \
	awk 'BEGIN{printf("McCortex Found:  %4d / %4d (%5.2f%%)\n",'$$C','$$NT',100*'$$C'/'$$NT')}';

run-cortex:
	cd cortex && $(MAKE) K=$(K)
	$(BCFTOOLS) isec truth.k$(K).norm.vcf.gz cortex/cortex.k$(K).norm.vcf.gz -p truthisec2
	@A=`$(VCFALLELES) truthisec2/0000.vcf`; \
	 B=`$(VCFALLELES) truthisec2/0001.vcf`; \
	 C=`$(VCFALLELES) truthisec2/0002.vcf`; \
	CT=`$(VCFALLELES) cortex/cortex.k$(K).norm.vcf.gz`; \
	NT=`$(VCFALLELES) truth.k$(K).norm.vcf.gz`; \
	awk 'BEGIN{printf("Cortex Missed: %4d / %4d (%5.2f%%)\n",'$$A','$$NT',100*'$$A'/'$$NT')}'; \
	awk 'BEGIN{printf("Cortex FP:     %4d / %4d (%5.2f%%)\n",'$$B','$$CT',100*'$$B'/'$$CT')}'; \
	awk 'BEGIN{printf("Cortex Found:  %4d / %4d (%5.2f%%)\n",'$$C','$$NT',100*'$$C'/'$$NT')}';

run-mccortex-brkpt: task.k$(K).mk $(TRUTH_FILES) $(READS)
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) breakpointvcf
	$(BCFTOOLS) isec truth.k$(K).norm.vcf.gz proj/vcfs/breakpoints.k$(K).vcf.gz -p truthisec
	@A=`$(VCFALLELES) truthisec/0000.vcf`; \
	 B=`$(VCFALLELES) truthisec/0001.vcf`; \
	 C=`$(VCFALLELES) truthisec/0002.vcf`; \
	MC=`$(VCFALLELES) proj/vcfs/breakpoints.k$(K).vcf.gz`; \
	NT=`$(VCFALLELES) truth.k$(K).norm.vcf.gz`; \
	awk 'BEGIN{printf("McCortex-brkpt Missed: %4d / %4d (%5.2f%%)\n",'$$A','$$NT',100*'$$A'/'$$NT')}'; \
	awk 'BEGIN{printf("McCortex-brkpt FP:     %4d / %4d (%5.2f%%)\n",'$$B','$$MC',100*'$$B'/'$$MC')}'; \
	awk 'BEGIN{printf("McCortex-brkpt Found:  %4d / %4d (%5.2f%%)\n",'$$C','$$NT',100*'$$C'/'$$NT')}';

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt task.k$(K).mk truth.*vcf* truthisec truthisec2
	cd cortex && $(MAKE) clean k=$(K)

.PHONY: all clean run-mccortex run-cortex
