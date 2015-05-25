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

K=31
HAPDEPTH=30
READLEN=150
ERRRATE=0.01
MEM=1G
OUTDIR=proj
REF=../data/chr22/chr22_17M_18M.fa

MASKS=diploid/mask0.fa diploid/mask1.fa
CHROMS_ALN=diploid/genome0.fa diploid/genome1.fa
CHROMS=diploid/chrom0.fa diploid/chrom1.fa
READS=reads/chrom0.$(HAPDEPTH)X.1.fa.gz reads/chrom0.$(HAPDEPTH)X.2.fa.gz \
			reads/chrom1.$(HAPDEPTH)X.1.fa.gz reads/chrom1.$(HAPDEPTH)X.2.fa.gz
DIRS=reads diploid

TRUTH_FILES=truth.k$(K).norm.vcf.gz truth.k$(K).norm.vcf.gz.csi

# Mark all dependencies as secondary
# It means don't re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:
# Delete files if their recipe fails
.DELETE_ON_ERROR:
# Remove implicit rules for certain suffixes
.SUFFIXES:

all: run-mccortex run-mccortex-brkpt $(TRUTH_FILES)

# 1 Mb human diploid
# Generate a diploid genome from a haploid reference
diploid/mask0.fa: diploid/genome0.fa
diploid/mask1.fa: diploid/genome1.fa
diploid/genome0.fa: diploid/genome1.fa
diploid/genome1.fa: $(REF) | $(DIRS)
	$(SIMMUT) --snps 1000 --indels 100 --invs 0 diploid 2 $<

$(CHROMS): $(CHROMS_ALN)

# Remove deletion marks (-) and convert to uppercase
diploid/chrom%.fa: diploid/genome%.fa
	cat $< | tr -d '-' | $(DNACAT) -u -F - > $@

truth.k$(K).vcf: $(CHROMS) $(MASKS)
	$(SIMTRUTH) $(REF) diploid/genome0.fa diploid/mask0.fa diploid/genome1.fa diploid/mask1.fa > $@

truth.k$(K).norm.vcf: truth.k$(K).vcf $(REF)
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref $(REF) --multiallelics +both $< > $@

%.vcf.gz: %.vcf
	$(BGZIP) -f $<

%.vcf.gz.csi: %.vcf.gz
	$(BCFTOOLS) index -f $<

# Simulate PE reads of each chrom each 50X
reads/chrom0.$(HAPDEPTH)X.1.fa.gz: reads/chrom0.$(HAPDEPTH)X.2.fa.gz
reads/chrom1.$(HAPDEPTH)X.1.fa.gz: reads/chrom1.$(HAPDEPTH)X.2.fa.gz

reads/chrom%.$(HAPDEPTH)X.2.fa.gz: diploid/chrom%.fa
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -e $(ERRRATE) reads/chrom$*.$(HAPDEPTH)X

samples.txt:
	echo "MissSample . "\
"reads/chrom0.$(HAPDEPTH)X.1.fa.gz:reads/chrom0.$(HAPDEPTH)X.2.fa.gz,"\
"reads/chrom1.$(HAPDEPTH)X.1.fa.gz:reads/chrom1.$(HAPDEPTH)X.2.fa.gz" > $@

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
