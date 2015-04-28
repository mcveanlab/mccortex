SHELL=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
READSIM=$(CTXDIR)/libs/readsim/readsim
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
SIMTRUTH=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_vcf.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip
VCFALLELES=$(CTXDIR)/scripts/bash/vcf-count-alleles

K=31
HAPDEPTH=30
READLEN=100
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

TRUTH_FILES=truth.norm.vcf.gz truth.norm.vcf.gz.csi

# Mark all dependencies as secondary
# It means don't re-run if the dependency file disappears -- allows us to delete unused files
.SECONDARY:

# Delete files if their recipe fails
.DELETE_ON_ERROR:

all: run-mccortex $(TRUTH_FILES)

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

truth.vcf: $(CHROMS) $(MASKS)
	$(SIMTRUTH) $(REF) diploid/genome0.fa diploid/mask0.fa diploid/genome1.fa diploid/mask1.fa > $@

truth.norm.vcf: truth.vcf $(REF)
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
	echo "MissSample . reads/chrom0.$(HAPDEPTH)X.1.fa.gz:reads/chrom1.$(HAPDEPTH)X.2.fa.gz" > $@

task.mk: samples.txt
	$(CTXPIPELINE) -r $(REF) $(K) $(OUTDIR) $< > $@

run-mccortex: task.mk $(TRUTH_FILES) $(READS)
	$(MAKE) -f $< CTXDIR=$(CTXDIR) bubblevcf
	$(BCFTOOLS) isec truth.norm.vcf.gz proj/vcfs/bubbles.k31.vcf.gz -p truthisec
	@A=`$(VCFALLELES) truthisec/0000.vcf`; \
	 B=`$(VCFALLELES) truthisec/0001.vcf`; \
	 C=`$(VCFALLELES) truthisec/0002.vcf`; \
	MC=`$(VCFALLELES) proj/vcfs/bubbles.k31.vcf.gz`; \
	NT=`$(VCFALLELES) truth.norm.vcf.gz`; \
	awk 'BEGIN{printf("McCortex Missed: %4d / %4d (%5.2f%%)\n",'$$A','$$NT',100*'$$A'/'$$NT')}'; \
	awk 'BEGIN{printf("McCortex FP:     %4d / %4d (%5.2f%%)\n",'$$B','$$MC',100*'$$B'/'$$MC')}'; \
	awk 'BEGIN{printf("McCortex Found:  %4d / %4d (%5.2f%%)\n",'$$C','$$NT',100*'$$C'/'$$NT')}';

run-cortex: cortex/cortex.norm.vcf.gz cortex/cortex.norm.vcf.gz.csi
	$(BCFTOOLS) isec truth.norm.vcf.gz cortex/cortex.norm.vcf.gz -p truthisec
	@A=`$(VCFALLELES) truthisec2/0000.vcf`; \
	 B=`$(VCFALLELES) truthisec2/0001.vcf`; \
	 C=`$(VCFALLELES) truthisec2/0002.vcf`; \
	CT=`$(VCFALLELES) cortex/cortex.norm.vcf.gz`; \
	NT=`$(VCFALLELES) truth.norm.vcf.gz`; \
	awk 'BEGIN{printf("Cortex Missed: %4d / %4d (%5.2f%%)\n",'$$A','$$NT',100*'$$A'/'$$NT')}'; \
	awk 'BEGIN{printf("Cortex FP:     %4d / %4d (%5.2f%%)\n",'$$B','$$CT',100*'$$B'/'$$CT')}'; \
	awk 'BEGIN{printf("Cortex Found:  %4d / %4d (%5.2f%%)\n",'$$C','$$NT',100*'$$C'/'$$NT')}';

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt task.mk truth.*vcf* truthisec

.PHONY: all clean run-mccortex run-cortex
