
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
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
VCFHEADER=$(CTXDIR)/libs/biogrok/vcf-header
VCFSORT=$(CTXDIR)/libs/biogrok/vcf-sort
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip

REF=$(shell cd ../../data/chr22/; echo `pwd`/chr22_17M_18M.fa)

K=31
DIRS=ref

# Dependencies that we create. We delete them on clean
CTXVCF=cortex_run/vcfs/chr22_17M_18M_union_BC_calls_k$(K).decomp.vcf

DEPS=ref/chr22_17M_18M.sthash ref/chr22_17M_18M.stidx \
     ref/ref.k$(K).ctx ref/ref.falist \
     samples.txt reads.1.falist reads.2.falist

TMP=cortex.k$(K).sort.vcf cortex.k$(K).norm.vcf

TGTS=cortex.k$(K).norm.vcf.gz cortex.k$(K).norm.vcf.gz.csi

all: $(TGTS)

ref/ref.falist: | $(DIRS)
	echo $(REF) > $@

# Make stampy hash
ref/chr22_17M_18M.stidx: $(REF) | $(DIRS)
	$(STAMPY) -G ref/chr22_17M_18M $(REF)

ref/chr22_17M_18M.sthash: ref/chr22_17M_18M.stidx $(REF) | $(DIRS)
	$(STAMPY) -g ref/chr22_17M_18M -H ref/chr22_17M_18M

# Build reference graph file
ref/ref.k$(K).ctx: $(REF) | $(DIRS)
	$(MCCORTEX) build -k $(K) -s REF -1 $< $@ >& $@.log

reads.%.falist:
	(echo `pwd`/../reads/chrom0.30X.$*.fa.gz; \
	 echo `pwd`/../reads/chrom1.30X.$*.fa.gz) > reads.$*.falist

samples.txt:
	printf "MrSample\t.\t%s\t%s\n" reads.1.falist reads.2.falist > $@

$(CTXVCF): $(DEPS)
	$(RUNCALLS) \
--first_kmer $(K) \
--last_kmer $(K) \
--kmer_step 2 \
--fastaq_index samples.txt \
--auto_cleaning yes \
--bc yes \
--pd no \
--outdir cortex_run \
--outvcf chr22_17M_18M \
--ploidy 2 \
--stampy_hash ref/chr22_17M_18M \
--stampy_bin $(STAMPY) \
--list_ref_fasta ref/ref.falist \
--refbindir ref/ \
--genome_size 1000000 \
--qthresh 5 \
--mem_height 20 --mem_width 100 \
--vcftools_dir $(VCFTOOLSDIR) \
--do_union yes \
--ref CoordinatesAndInCalling \
--workflow independent \
--logfile runcalls.k$(K).log

# Add '##contig=<ID=chr22_17M_18M,length=1000000,assembly=hg19>'
# to header, and fix an INFO field
new_header.k$(K).txt: $(CTXVCF)
	( REFNAME=`$(DNACAT) -N $(REF) | awk '{print $$1}'`; \
		REFLEN=`$(DNACAT) -L $(REF) | awk '{print $2}'`; \
		$(VCFHEADER) $(CTXVCF) | \
	  grep -v '^##contig' | \
	  grep -v '^#CHROM' | \
	  sed 's/, Description=/,Description=/g'; \
	  echo '##INFO=<ID=KMER,Number=1,Type=Integer,Description="Kmer used for calling">'; \
	  echo "##contig=<ID=$$REFNAME,length=$$REFLEN,assembly=$$REFNAME>"; \
	  $(VCFHEADER) $(CTXVCF) | grep '^#CHROM' ) > new_header.k$(K).txt

# Put new header on, filter ref mismatches, sort entries
cortex.k$(K).sort.vcf: new_header.k$(K).txt $(CTXVCF)
	( cat new_header.k$(K).txt; \
	  $(VCFREF) -s $(CTXVCF) $(REF) | grep -v '^#' | sort -k1,1d -k2,2n ) > $@

# Normalise indels and compress
cortex.k$(K).norm.vcf.gz: cortex.k$(K).sort.vcf
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref $(REF) --multiallelics +both $< > cortex.k$(K).norm.vcf
	$(BGZIP) cortex.k$(K).norm.vcf

# Index
%.norm.vcf.gz.csi: %.norm.vcf.gz
	$(BCFTOOLS) index $<

$(DIRS):
	mkdir -p $@

clean:
	rm -rf ref $(DEPS) $(CTXVCF) $(TGTS) $(TMP) new_header.k$(K).txt runcalls.k$(K).log cortex_run

.PHONY: all clean
