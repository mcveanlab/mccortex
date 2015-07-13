
SHELL=/bin/bash -euo pipefail

# Need to install these separately yourself
CORTEXDIR=~/cortex/releases/CORTEX_release_v1.0.5.21
STAMPY=/apps/well/stampy/1.0.23-py2.6/stampy.py
VCFTOOLSDIR=~/bioinf/vcftools_0.1.12b/
VCFREF=~/c/vcf-hack/bin/vcfref
#

RUNCALLS=$(CORTEXDIR)/scripts/calling/run_calls.pl
CORTEX=$(CORTEXDIR)/bin/cortex_var_31_c1

CTXDIR=../../../
MCCORTEX=$(CTXDIR)/bin/mccortex31
VCFHEADER=$(CTXDIR)/libs/biogrok/vcf-header
VCFSORT=$(CTXDIR)/libs/biogrok/vcf-sort
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip

# Need absolute path to reference
REF=$(shell cd ../../data/ecoli/; echo `pwd`/ecoli.fa)

K=31
DIRS=ref lists

# Dependencies that we create. We delete them on clean
FALISTS=$(shell echo lists/reads{0..9}.1.falist lists/reads{0..9}.2.falist)

DEPS=ref/ecoli.sthash ref/ecoli.stidx \
     ref/ref.k$(K).ctx ref/ref.falist \
     samples.txt $(FALISTS)

CTXVCF_1BY1=cortex_1by1/vcfs/ecoli_union_BC_calls_k$(K).decomp.vcf
CTXVCF_JOINT=cortex_joint/vcfs/ecoli_union_BC_calls_k$(K).decomp.vcf
VCFS=cortex.k$(K).1by1.norm.vcf.gz cortex.k$(K).joint.norm.vcf.gz
CSIS=$(VCFS:.vcf.gz=.vcf.gz.csi)

TGTS=$(CTXVCF_1BY1) $(CTXVCF_JOINT) $(VCFS) $(CSIS)

all: $(TGTS) | checks

ref/ref.falist: | $(DIRS)
	echo $(REF) > $@

# Make stampy hash
ref/ecoli.stidx: $(REF) | $(DIRS)
	$(STAMPY) -G ref/ecoli $(REF)

ref/ecoli.sthash: ref/ecoli.stidx $(REF) | $(DIRS)
	$(STAMPY) -g ref/ecoli -H ref/ecoli

# Build reference graph file
ref/ref.k$(K).ctx: $(REF) | $(DIRS)
	$(MCCORTEX) build -k $(K) -s REF -1 $< $@ >& $@.log

lists/reads%.1.falist: lists
	echo `pwd`/../reads/chrom$*.50X.1.fa.gz > $@
lists/reads%.2.falist: lists
	echo `pwd`/../reads/chrom$*.50X.2.fa.gz > $@

samples.txt:
	for i in {0..9}; do \
		printf "Ecoli$$i\t.\t%s\t%s\n" lists/reads$$i.1.falist lists/reads$$i.2.falist; \
  done > $@

COMMON_ARGS= \
--first_kmer $(K) \
--last_kmer $(K) \
--kmer_step 2 \
--fastaq_index samples.txt \
--auto_cleaning yes \
--bc yes \
--pd no \
--outvcf ecoli \
--ploidy 2 \
--stampy_hash $(REF) \
--stampy_bin $(STAMPY) \
--list_ref_fasta ref/ref.falist \
--refbindir ref/ \
--genome_size 1000000 \
--qthresh 5 \
--mem_height 20 --mem_width 100 \
--vcftools_dir $(VCFTOOLSDIR) \
--do_union yes \
--ref CoordinatesAndInCalling

$(CTXVCF_1BY1): $(DEPS)
	$(RUNCALLS) $(COMMON_ARGS) \
--outdir cortex_1by1 \
--workflow independent \
--logfile runcalls.1by1.k$(K).log

$(CTXVCF_JOINT): $(DEPS)
	$(RUNCALLS) $(COMMON_ARGS) \
--outdir cortex_joint \
--workflow joint \
--logfile runcalls.joint.k$(K).log

# Add '##contig=<ID=ecoli,length=1000000,assembly=hg19>'
# to header, and fix an INFO field
new_header.k$(K).txt: cortex_1by1/vcfs/ecoli_union_BC_calls_k$(K).decomp.vcf
	( REFNAME=`$(DNACAT) -N $(REF) | awk '{print $$1}'`; \
		REFLEN=`$(DNACAT) -L $(REF) | awk '{print $2}'`; \
		$(VCFHEADER) $< | \
	  grep -v '^##contig' | \
	  grep -v '^#CHROM' | \
	  sed 's/, Description=/,Description=/g'; \
	  echo '##INFO=<ID=KMER,Number=1,Type=Integer,Description="Kmer used for calling">'; \
	  echo "##contig=<ID=$$REFNAME,length=$$REFLEN,assembly=$REFNAME>'; \
	  $(VCFHEADER) $< | grep '^#CHROM' ) > new_header.k$(K).txt

# Put new header on, filter ref mismatches, sort entries
cortex.k$(K).%.sort.vcf: cortex_%/vcfs/ecoli_union_BC_calls_k$(K).decomp.vcf new_header.k$(K).txt
	( cat new_header.k$(K).txt; \
	  $(VCFREF) -s $< $(REF) | grep -v '^#' | sort -k1,1d -k2,2n ) > $@

# Normalise indels and compress
cortex.k$(K).%.norm.vcf.gz: cortex.k$(K).%.sort.vcf
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref $(REF) --multiallelics +both $< > cortex.k$(K).$*.norm.vcf
	$(BGZIP) cortex.k$(K).$*.norm.vcf

# Index
%.norm.vcf.gz.csi: %.norm.vcf.gz
	$(BCFTOOLS) index $<

$(DIRS):
	mkdir -p $@

# Check there are executables for cortex
checks:
	@[ -x $(CORTEXDIR)/bin/cortex_var_31_c1 ] || ( echo 'Error: need cortex_var_31_c11' 1>&2 && false )
	@[ -x $(CORTEXDIR)/bin/cortex_var_31_c2 ] || ( echo 'Error: need cortex_var_31_c11' 1>&2 && false )
	@[ -x $(CORTEXDIR)/bin/cortex_var_31_c10 ] || ( echo 'Error: need cortex_var_31_c10' 1>&2 && false )
	@[ -x $(CORTEXDIR)/bin/cortex_var_31_c11 ] || ( echo 'Error: need cortex_var_31_c11' 1>&2 && false )

clean:
	rm -rf ref $(DEPS) $(DIRS) $(TGTS) $(TMP) new_header.k$(K).txt
	rm -rf runcalls.{1by1,joint}.k$(K).log cortex_{1by1,joint} cortex.k$(K).*

.PHONY: all clean checks
