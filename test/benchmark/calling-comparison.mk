
# Example input
# SEQ=chr1.fa
NUM_INDIVS=1
PLOIDY=2
KMER=31
# SNPS=
# INDELS=
# INV=
# INVLEN=
READLEN=100
MPSIZE=250
ALLELECOVG=30
# ERRPROF=
MEMWIDTH=20
MEMHEIGHT=20
# GENOMESIZE=

SHELL := /bin/bash

# current_dir = $(shell pwd)
current_dir := $(dir $(lastword $(MAKEFILE_LIST)))

CORTEX_PATH=$(HOME)/cortex/releases/CORTEX_release_v1.0.5.20
CTX_PATH=$(current_dir)/../../

# External tools
BCFTOOLS=$(HOME)/bioinf/bcftools/bcftools
STAMPY=$(HOME)/bioinf/stampy-1.0.20/stampy.py
SAMTOOLS=samtools
VCFTOOLSDIR=$(HOME)/bioinf/vcftools_0.1.11

UNAME=$(shell uname -s)
ifeq ($(UNAME),Darwin)
	STAMPY_BIN=python2.6 $(STAMPY)
else
	STAMPY_BIN=python $(STAMPY)
endif

$(shell (echo '#!/bin/bash'; echo '$(STAMPY_BIN) $$@';) > stampy.sh; chmod +x stampy.sh)
STAMPY_BIN=./stampy.sh

NUMCOLS=$(shell echo $$(($(NUM_INDIVS)+1)))
MEM=$(shell bc <<< '( $(MEMWIDTH) * 2^$(MEMHEIGHT) * (8+8+4)*8+$(NUMCOLS) ) / 8')

RELEASECTX=$(CORTEX_PATH)/bin/cortex_var_31_c$(NUMCOLS) --kmer_size $(KMER) --mem_height $(MEMHEIGHT) --mem_width $(MEMWIDTH)
BUILDCTX=$(CTX_PATH)/bin/ctx31 build
CLEANCTX=$(CTX_PATH)/bin/ctx31 clean
JOINCTX=$(CTX_PATH)/bin/ctx31 join
INFERCTX=$(CTX_PATH)/bin/ctx31 inferedges --pop
THREADCTX=$(CTX_PATH)/bin/ctx31 thread
CALLCTX=$(CTX_PATH)/bin/ctx31 call
PROCCTX=$(CTX_PATH)/bin/ctx31 unique
PLACECTX=$(CTX_PATH)/bin/ctx31 place
TRAVERSE=$(CTX_PATH)/bin/traversal31

RUNCALLS=time $(CORTEX_PATH)/scripts/calling/run_calls.pl

BIOINF=$(CTX_PATH)/libs/bioinf-perl
READSIM=$(CTX_PATH)/libs/readsim/readsim
SEQCAT=$(CTX_PATH)/libs/seq_file/bin/seqcat
HAPLEN=$(CTX_PATH)/scripts/longest-haplotype.sh
OLDCLEAN=$(CTX_PATH)/scripts/clean_bubbles.pl

# Measure genome size if not passed
ifndef GENOMESIZE
	GENOMESIZE=$(shell $(SEQCAT) $(SEQ) | tr -d '\n' | wc | grep -o '[0-9]*$$')
endif

# Calculate some numbers
NCHROMS=$(shell bc <<< '$(NUM_INDIVS) * $(PLOIDY)')
LASTCHROM=$(shell bc <<< '$(NCHROMS) - 1')
LASTINDIV=$(shell bc <<< '$(NUM_INDIVS) - 1')
NINDIVS_REF=$(shell echo $$(($(NUM_INDIVS) + 1)))

SEQPATH=$(realpath $(SEQ))

# Generate file names
GENOMES=$(shell echo genomes/genome{0..$(LASTCHROM)}.fa)
READS=$(shell echo reads/reads{0..$(LASTCHROM)}.{1..2}.fa.gz)
RAWGRAPHS=$(shell echo k$(KMER)/graphs/sample{0..$(LASTINDIV)}.raw.ctx)
CLEANGRAPHS=$(RAWGRAPHS:.raw.ctx=.clean.ctx)

DOOLD=$(shell if [ $(NUM_INDIVS) -gt 10 ]; then echo 'no'; else echo 'yes'; fi)

ifeq ($(DOOLD),yes)
	SET={oldbc,newbc,se,pe,sepe}
else
	# don't do oldbc if more than 10 individuals
	SET={newbc,se,pe,sepe}
endif

# Can't do noref if we only have one chrom
ifeq ($(NCHROMS),1)
	# Always need ref
	PATHS=$(shell echo k$(KMER)/graphs/pop.{se,pe,sepe}.ref.ctp)
	BUBBLES=$(shell echo `eval echo k$(KMER)/bubbles/samples.$(SET).ref.bubbles.gz`)
	BUBBLEVCFS=$(shell echo `eval echo k$(KMER)/vcfs/samples.$(SET).ref.bub.vcf`)
	TRUTHBUBBLES=k$(KMER)/vcfs/truth.ref.bub.vcf
	BUBOLDCMPRULES=compare-old-ref-bubbles
	BUBNEWCMPRULES=compare-new-ref-bubbles
	NORMCMPRULES=$(shell echo `eval echo compare-$(SET).ref-norm compare-runcalls-norm`)
else
	PATHS=$(shell echo k$(KMER)/graphs/pop.{se,pe,sepe}.{noref,ref}.ctp)
	BUBBLES=$(shell echo `eval echo k$(KMER)/bubbles/samples.$(SET).{noref,ref}.bubbles.gz`)
	BUBBLEVCFS=$(shell echo `eval echo k$(KMER)/vcfs/samples.$(SET).{noref,ref}.bub.vcf`)
	TRUTHBUBBLES=k$(KMER)/vcfs/truth.ref.bub.vcf k$(KMER)/vcfs/truth.noref.bub.vcf
	# don't do BUBOLDCMPRULES
	BUBNEWCMPRULES=$(shell echo compare-new-{noref,ref}-bubbles)
	NORMCMPRULES=$(shell echo `eval echo compare-$(SET).{noref,ref}-norm compare-runcalls-norm`)
endif

READLISTS=$(shell echo reads/reads{0..$(LASTINDIV)}.{1,2}.falist)
MGLIST=$(shell for i in {0..$(LASTCHROM)}; do echo -n " genomes/genome$$i.fa genomes/mask$$i.fa"; done)

se_list=$(shell for i in `seq 0 $(LASTINDIV)`; do \
	echo -n " --col $$i $$i"; \
	for j in `seq $$(($$i * $(PLOIDY))) $$(($$i * $(PLOIDY) + $(PLOIDY) - 1))`; do \
		echo -n " --seq reads/reads$$j.1.fa.gz --seq reads/reads$$j.2.fa.gz"; \
	done; \
done)

pe_list=$(shell for i in `seq 0 $(LASTINDIV)`; do \
	echo -n " --col $$i $$i"; \
	for j in `seq $$(($$i * $(PLOIDY))) $$(($$i * $(PLOIDY) + $(PLOIDY) - 1))`; do \
		echo -n " --seq2 reads/reads$$j.1.fa.gz reads/reads$$j.2.fa.gz"; \
	done; \
done)

sepe_list=$(shell for i in `seq 0 $(LASTINDIV)`; do \
	echo -n " --col $$i $$i"; \
	for j in `seq $$(($$i * $(PLOIDY))) $$(($$i * $(PLOIDY) + $(PLOIDY) - 1))`; do \
		echo -n " --seq reads/reads$$j.1.fa.gz --seq reads/reads$$j.2.fa.gz"; \
		echo -n " --seq2 reads/reads$$j.1.fa.gz reads/reads$$j.2.fa.gz"; \
	done; \
done)

PLACEVCFS=$(BUBBLEVCFS:.bub.vcf=.decomp.vcf)
PASSVCFS=$(BUBBLEVCFS:.bub.vcf=.pass.vcf)
NORMVCFS=$(BUBBLEVCFS:.bub.vcf=.norm.vcf)

FLANKFILES=$(BUBBLEVCFS:.vcf=.5pflanks.fa.gz)
SAMFILES=$(BUBBLEVCFS:.vcf=.5pflanks.sam)

ifdef ERRPROF
	GRAPHS_noref=$(CLEANGRAPHS)
else
	GRAPHS_noref=$(RAWGRAPHS)
endif

GRAPHS_ref=$(GRAPHS_noref) ref/ref.k$(KMER).ctx
REFARGS=--ref $(NUM_INDIVS)

GENOMES_noref=$(GENOMES)
GENOMES_ref=$(GENOMES) $(SEQ)

KEEP=$(GENOMES) $(READS) $(PATHS) $(BUBBLES) $(PASSVCFS) $(PLACEVCFS) $(NORMVCFS)

all: checkcmds $(KEEP) compare-bubbles compare-normvcf traverse

checkcmds:
	@if [ '$(SEQ)' == '' ]; then echo "You need to specify SEQ=.. Please and thank you."; exit -1; fi;

test:
	@echo DOOLD: $(DOOLD) $(NUM_INDIVS)
	@echo KEEP: $(KEEP)

compare-bubbles: $(BUBBLEVCFS) $(TRUTHBUBBLES) $(BUBOLDCMPRULES) $(BUBNEWCMPRULES)

# % is noref or ref
$(BUBOLDCMPRULES): compare-old-%-bubbles: k$(KMER)/vcfs/samples.oldbc.%.bub.vcf k$(KMER)/vcfs/truth.%.bub.vcf 
	@echo == Released Cortex $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.$*.bub.vcf k$(KMER)/vcfs/samples.oldbc.$*.bub.vcf k$(KMER)/vcfs/truth.oldbc.$*.vcf OLDBC k$(KMER)/vcfs/falsepos.oldbc.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.oldbc.$*.bub.vcf

# % is noref or ref
$(BUBNEWCMPRULES): $(BUBBLEVCFS)
$(BUBNEWCMPRULES): compare-new-%-bubbles: k$(KMER)/vcfs/truth.%.bub.vcf k$(KMER)/vcfs/samples.newbc.%.bub.vcf k$(KMER)/vcfs/samples.se.%.bub.vcf k$(KMER)/vcfs/samples.pe.%.bub.vcf k$(KMER)/vcfs/samples.sepe.%.bub.vcf
	@echo == New Bubble Caller $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.$*.bub.vcf k$(KMER)/vcfs/samples.newbc.$*.bub.vcf k$(KMER)/vcfs/truth.newbc.$*.vcf NEWBC k$(KMER)/vcfs/falsepos.newbc.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.newbc.$*.bub.vcf
	@echo == Paths se $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.newbc.$*.vcf k$(KMER)/vcfs/samples.se.$*.bub.vcf k$(KMER)/vcfs/truth.se.$*.vcf PAC k$(KMER)/vcfs/falsepos.se.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.se.$*.bub.vcf
	@echo == Paths pe $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.se.$*.vcf k$(KMER)/vcfs/samples.pe.$*.bub.vcf k$(KMER)/vcfs/truth.pe.$*.vcf PAC k$(KMER)/vcfs/falsepos.pe.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.pe.$*.bub.vcf
	@echo == Paths sepe $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.pe.$*.vcf k$(KMER)/vcfs/samples.sepe.$*.bub.vcf k$(KMER)/vcfs/truth.sepe.$*.vcf PAC k$(KMER)/vcfs/falsepos.sepe.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.sepe.$*.bub.vcf
	@echo == Truth ==
	$(HAPLEN) k$(KMER)/vcfs/truth.$*.bub.vcf

compare-normvcf: $(NORMVCFS) k$(KMER)/vcfs/truth.norm.vcf $(NORMCMPRULES)
$(NORMCMPRULES): compare-%-norm: k$(KMER)/vcfs/samples.%.norm.vcf
	@echo == $< ==
	$(BIOINF)/vcf_scripts/vcf_isec.pl k$(KMER)/vcfs/truth.norm.vcf $< > /dev/null

traverse: $(PATHS) k$(KMER)/graphs/pop.ref.ctx
	$(TRAVERSE)                                     --nsamples 10000 k$(KMER)/graphs/pop.ref.ctx
	$(TRAVERSE) -p k$(KMER)/graphs/pop.se.ref.ctp   --nsamples 10000 k$(KMER)/graphs/pop.ref.ctx
	$(TRAVERSE) -p k$(KMER)/graphs/pop.pe.ref.ctp   --nsamples 10000 k$(KMER)/graphs/pop.ref.ctx
	$(TRAVERSE) -p k$(KMER)/graphs/pop.sepe.ref.ctp --nsamples 10000 k$(KMER)/graphs/pop.ref.ctx

ref/stampy.stidx:
	mkdir -p ref
	$(STAMPY_BIN) -G ref/stampy $(SEQ)

ref/stampy.sthash:
	mkdir -p ref
	$(STAMPY_BIN) -g ref/stampy -H ref/stampy

$(SEQ).fai:
	samtools faidx $(SEQ)

DIRS=ref genomes reads k$(KMER) runcalls

clean:
	rm -rf $(DIRS) gap_sizes.*.csv mp_sizes.*.csv stampy.sh

#
# Patterns
#

$(GENOMES):
	mkdir -p genomes
	$(BIOINF)/sim_mutations/sim_mutations.pl --snps $(SNPS) --indels $(INDELS) --invs $(INV) --invlen $(INVLEN) genomes/ $(NCHROMS) $(SEQ)

$(READS): $(GENOMES)

$(RAWGRAPHS): $(READS)

$(CLEANGRAPHS): $(RAWGRAPHS)

# Reads
reads/reads%.1.fa.gz reads/reads%.2.fa.gz: genomes/genome%.fa
	mkdir -p reads
	cat genomes/genome$*.fa | tr -d '-' | $(READSIM) -r - -i $(MPSIZE) -v 0.2 -l $(READLEN) -d $(ALLELECOVG) $(USECALIB) reads/reads$*

reads/reads%.1.falist reads/reads%.2.falist:
	mkdir -p reads
	echo -n '' > reads/reads$*.1.falist; echo -n '' > reads/reads$*.2.falist;
	a=$$(($* * $(PLOIDY))); b=$$(($$a + $(PLOIDY) - 1)); \
	for i in `seq $$a $$b`; do \
		echo reads$$i.1.fa.gz >> reads/reads$*.1.falist; \
		echo reads$$i.2.fa.gz >> reads/reads$*.2.falist; \
	done

k$(KMER)/graphs/sample%.clean.ctx: k$(KMER)/graphs/sample%.raw.ctx
	$(CLEANCTX) $@ $<

k$(KMER)/graphs/sample%.raw.ctx: $(READS)
	mkdir -p k$(KMER)/graphs
	a=$$(($* * $(PLOIDY))); b=$$(($$a+$(PLOIDY)-1)); \
	files=$$(for k in `seq $$a $$b`; do echo -n " --seq2 reads/reads$$k.1.fa.gz reads/reads$$k.2.fa.gz"; done); \
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample Sample$* $$files k$(KMER)/graphs/sample$*.raw.ctx;

k$(KMER)/graphs/pop.noref.ctx: $(GRAPHS_noref)
k$(KMER)/graphs/pop.ref.ctx: $(GRAPHS_ref)
k$(KMER)/graphs/pop.%.ctx:
	$(JOINCTX) -m $(MEM) $@ $(GRAPHS_$*)
	$(INFERCTX) $@

# Paths
$(PATHS): k$(KMER)/graphs/pop.noref.ctx k$(KMER)/graphs/pop.ref.ctx

k$(KMER)/graphs/pop.%.noref.ctp: k$(KMER)/graphs/pop.noref.ctx
	$(THREADCTX) -t 1 $($*_list) $(NUM_INDIVS) $@ $< $<
	for f in *_sizes.*.csv; do mv $$f k$(KMER)/graphs/se.$$f; done

k$(KMER)/graphs/pop.%.ref.ctp: k$(KMER)/graphs/pop.ref.ctx
	$(THREADCTX) -t 1 $($*_list) --col $(NUM_INDIVS) $(NUM_INDIVS) --seq $(SEQ) $(NUMCOLS) $@ $< $<
	for f in *_sizes.*.csv; do mv $$f k$(KMER)/graphs/se.$$f; done

# Bubbles
$(BUBBLES): k$(KMER)/graphs/pop.noref.ctx k$(KMER)/graphs/pop.ref.ctx $(READLISTS)

k$(KMER)/bubbles/samples.oldbc.%.bubbles.gz: k$(KMER)/graphs/pop.%.ctx
	mkdir -p k$(KMER)/bubbles
	callargs=`if [ '$*' == 'ref' ]; then echo '--ref_colour $(NUM_INDIVS)'; fi`; \
	time $(RELEASECTX) --multicolour_bin $< $$callargs --detect_bubbles1 -1/-1 --output_bubbles1 k$(KMER)/bubbles/samples.oldbc.$*.bubbles --print_colour_coverages
	mv k$(KMER)/bubbles/samples.oldbc.$*.bubbles k$(KMER)/bubbles/samples.oldbc.$*.bubbles.dirty
	$(OLDCLEAN) $(KMER) k$(KMER)/bubbles/samples.oldbc.$*.bubbles.dirty | gzip -c > $@

k$(KMER)/bubbles/samples.newbc.%.bubbles.gz: k$(KMER)/graphs/pop.%.ctx
	mkdir -p k$(KMER)/bubbles
	callargs=`if [ $* == 'ref' ]; then echo '--ref $(NUM_INDIVS)'; fi`; \
	$(CALLCTX) -t 1 $$callargs $< $@

# % => {se,pe,sepe}.{ref.noref}
k$(KMER)/bubbles/samples.%.bubbles.gz: k$(KMER)/graphs/pop.%.ctp
	mkdir -p k$(KMER)/bubbles
	r=`echo $@ | grep -oE '(no)?ref'`; \
	$(CALLCTX) -t 1 -m $(MEM) $(REFARGS) -p $< k$(KMER)/graphs/pop.$$r.ctx $@

k$(KMER)/vcfs/truth.noref.bub.vcf: $(GENOMES)
	$(BIOINF)/sim_mutations/sim_bubble_vcf.pl $(KMER) $(MGLIST) > $@

k$(KMER)/vcfs/truth.ref.bub.vcf: $(GENOMES)
	$(BIOINF)/sim_mutations/sim_bubble_vcf.pl --ref $(SEQ) $(KMER) $(MGLIST) > $@

k$(KMER)/vcfs/truth.decomp.vcf: $(SEQ) $(GENOMES)
	$(BIOINF)/sim_mutations/sim_decomp_vcf.pl $(SEQ) $(MGLIST) > k$(KMER)/vcfs/truth.decomp.vcf

k$(KMER)/vcfs/samples.%.bub.vcf k$(KMER)/vcfs/samples.%.bub.5pflanks.fa.gz: k$(KMER)/bubbles/samples.%.bubbles.gz
	mkdir -p k$(KMER)/vcfs
	$(PROCCTX) k$(KMER)/bubbles/samples.$*.bubbles.gz k$(KMER)/vcfs/samples.$*.bub
	gzip -d -f k$(KMER)/vcfs/samples.$*.bub.vcf.gz

$(SAMFILES): ref/stampy.stidx ref/stampy.sthash
k$(KMER)/vcfs/samples.%.bub.5pflanks.sam: k$(KMER)/vcfs/samples.%.bub.5pflanks.fa.gz
	mkdir -p ref
	$(STAMPY_BIN) -g ref/stampy -h ref/stampy --inputformat=fasta -M $< > $@

k$(KMER)/vcfs/samples.%.decomp.vcf: k$(KMER)/vcfs/samples.%.bub.vcf k$(KMER)/vcfs/samples.%.bub.5pflanks.sam
	$(PLACECTX) k$(KMER)/vcfs/samples.$*.bub.vcf k$(KMER)/vcfs/samples.$*.bub.5pflanks.sam $(SEQ) > $@

k$(KMER)/vcfs/samples.%.pass.vcf: k$(KMER)/vcfs/samples.%.decomp.vcf
	cat $< | awk '$$1 ~ /^#/ || $$7 ~ /^(PASS|\.)$$/' | vcf-sort > $@

reads/reads.index: $(READLISTS)
	for i in {0..$(LASTINDIV)}; do echo -e Sample$$i"\t"."\t"reads/reads$$i.1.falist"\t"reads/reads$$i.2.falist; done > reads/reads.index

ref/ref.falist:
	mkdir -p ref
	echo $(SEQPATH) > ref/ref.falist

ref/ref.k$(KMER).ctx:
	mkdir -p ref
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample ref --seq $(SEQ) ref/ref.k$(KMER).ctx

k$(KMER)/vcfs/samples.runcalls.norm.vcf: reads/reads.index ref/ref.falist ref/ref.k$(KMER).ctx $(CORTEX_PATH)/bin/cortex_var_31_c2 $(CORTEX_PATH)/bin/cortex_var_31_c$(NINDIVS_REF)
	rm -rf runcalls/runcalls.log
	mkdir -p runcalls
	$(RUNCALLS) --first_kmer $(KMER) --last_kmer $(KMER) \
	            --fastaq_index reads/reads.index --auto_cleaning yes \
	            --mem_width $(MEMWIDTH) --mem_height $(MEMHEIGHT) \
						  --ploidy $(PLOIDY) --bc yes --pd no \
						  --outdir runcalls --outvcf samples.runcalls \
						  --stampy_hash ref/stampy --stampy_bin '$(STAMPY_BIN)' \
						  --list_ref_fasta ref/ref.falist --refbindir ref \
						  --genome_size $(GENOMESIZE) \
						  --qthresh 5 --vcftools_dir $(VCFTOOLSDIR) \
						  --do_union yes --ref CoordinatesAndInCalling \
						  --workflow independent --logfile runcalls/runcalls.log
						  # --apply_pop_classifier
	cp runcalls/vcfs/samples.runcalls_union_BC_calls_k31.decomp.vcf k$(KMER)/vcfs/samples.runcalls.norm.vcf

k$(KMER)/vcfs/truth.norm.vcf: k$(KMER)/vcfs/truth.decomp.vcf
	$(BCFTOOLS) norm --remove-duplicate -f $(SEQ) $< > $@

$(NORMVCFS): $(SEQ).fai

k$(KMER)/vcfs/samples.%.norm.vcf: k$(KMER)/vcfs/samples.%.pass.vcf
	$(BCFTOOLS) norm --remove-duplicate -f $(SEQ) $< > k$(KMER)/vcfs/samples.$*.norm.vcf


.PHONY: all clean test checkcmds
.PHONY: compare-bubbles compare-old-bubbles compare-new-bubbles
.PHONY: compare-normvcf $(NORMCMPRULES)
.PHONY: traverse
