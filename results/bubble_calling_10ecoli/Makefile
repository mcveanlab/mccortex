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
VCFISEC=$(CTXDIR)/scripts/make-isec.sh

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
	$(CTXPIPELINE) --stampy --stampy-base $(REF) -r $(REF) $(K) $(OUTDIR) $< > $@

run-mccortex: task.k$(K).mk $(TRUTH_FILES) $(READS)
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) JOINT_CALLING=yes USE_LINKS=yes
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) JOINT_CALLING=yes USE_LINKS=no
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) JOINT_CALLING=no  USE_LINKS=yes
	$(MAKE) -f $< CTXDIR=$(CTXDIR) MEM=$(MEM) JOINT_CALLING=no  USE_LINKS=no 
	@echo "McCortex Bubbles (joint+links) found:"
	$(VCFISEC) mc_isec_bub_joint_links truth.k$(K).norm.vcf.gz proj/vcfs/bubbles.joint.links.k$(K).vcf.gz
	$(VCFISEC) mc_isec_bub_joint_plain truth.k$(K).norm.vcf.gz proj/vcfs/bubbles.joint.plain.k$(K).vcf.gz
	$(VCFISEC) mc_isec_bub_1by1_links  truth.k$(K).norm.vcf.gz proj/vcfs/bubbles.1by1.links.k$(K).vcf.gz
	$(VCFISEC) mc_isec_bub_1by1_plain  truth.k$(K).norm.vcf.gz proj/vcfs/bubbles.1by1.plain.k$(K).vcf.gz
	$(VCFISEC) mc_isec_brk_joint_links truth.k$(K).norm.vcf.gz proj/vcfs/breakpoints.joint.links.k$(K).vcf.gz
	$(VCFISEC) mc_isec_brk_joint_plain truth.k$(K).norm.vcf.gz proj/vcfs/breakpoints.joint.plain.k$(K).vcf.gz
	$(VCFISEC) mc_isec_brk_1by1_links  truth.k$(K).norm.vcf.gz proj/vcfs/breakpoints.1by1.links.k$(K).vcf.gz
	$(VCFISEC) mc_isec_brk_1by1_plain  truth.k$(K).norm.vcf.gz proj/vcfs/breakpoints.1by1.plain.k$(K).vcf.gz

run-cortex:
	cd cortex && $(MAKE) K=$(K)
	@echo "Cortex Bubble Caller 1by1 found:"
	$(VCFISEC) ctx_isec_bub_1by1 truth.k$(K).norm.vcf.gz cortex/cortex.1by1.k$(K).norm.vcf.gz
	$(VCFISEC) ctx_isec_bub_joint truth.k$(K).norm.vcf.gz cortex/cortex.joint.k$(K).norm.vcf.gz

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) $(OUTDIR) samples.txt task.k$(K).mk truth.*vcf*
	rm -rf mc_isec_* ctx_isec_*
	cd cortex && $(MAKE) clean k=$(K)

.PHONY: all clean run-mccortex run-cortex
