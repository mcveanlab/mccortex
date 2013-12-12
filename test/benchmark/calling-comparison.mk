
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
# MAPARGS=
MINMAPQ=40
MAXALLELE=500
NTHREADS=4
# How many contigs to print
NCONTIGS=10000

SHELL := /bin/bash

# External tools not in repo
STAMPY=$(HOME)/bioinf/stampy-1.0.20/stampy.py
SAMTOOLS=samtools
VCFTOOLSDIR=$(HOME)/bioinf/vcftools_0.1.11
CORTEX_PATH=$(HOME)/cortex/releases/CORTEX_release_v1.0.5.20

# current_dir = $(shell pwd)
current_dir := $(dir $(lastword $(MAKEFILE_LIST)))
CTX_PATH=$(current_dir)/../../

UNAME=$(shell uname -s)
ifeq ($(UNAME),Darwin)
	STAMPY_BIN=python2.6 $(STAMPY)
else
	STAMPY_BIN=python $(STAMPY)
endif

$(shell (echo '#!/bin/bash'; echo '$(STAMPY_BIN) $$@';) > stampy.sh; chmod +x stampy.sh)
STAMPY_BIN=./stampy.sh

NUMCOLS=$(shell echo $$(($(NUM_INDIVS)+1)))
MEM=$(shell bc <<< '( $(MEMWIDTH) * 2^$(MEMHEIGHT) * (8+8+4)*8+$(NUMCOLS) ) / 4')

RELEASECTX=$(CORTEX_PATH)/bin/cortex_var_31_c$(NUMCOLS) --kmer_size $(KMER) --mem_height $(MEMHEIGHT) --mem_width $(MEMWIDTH)
BUILDCTX=$(CTX_PATH)/bin/ctx31 build
CLEANCTX=$(CTX_PATH)/bin/ctx31 clean
JOINCTX=$(CTX_PATH)/bin/ctx31 join
INFERCTX=$(CTX_PATH)/bin/ctx31 inferedges --pop
THREADCTX=$(CTX_PATH)/bin/ctx31 thread
CALLCTX=$(CTX_PATH)/bin/ctx31 call
PROCCTX=$(CTX_PATH)/bin/ctx31 unique
PLACECTX=$(CTX_PATH)/bin/ctx31 place
TRAVERSE=$(CTX_PATH)/bin/ctx31 contigs
CTXSTATS=$(CTX_PATH)/scripts/cortex_stats.pl

RUNCALLS=time $(CORTEX_PATH)/scripts/calling/run_calls.pl

BIOINF=$(CTX_PATH)/libs/bioinf-perl
READSIM=$(CTX_PATH)/libs/readsim/readsim
SEQCAT=$(CTX_PATH)/libs/seq_file/bin/seqcat
FACAT=$(CTX_PATH)/libs/seq_file/bin/facat
BCFTOOLS=$(CTX_PATH)/libs/bcftools/bcftools
VCFLIBALIGN=$(CTX_PATH)/libs/vcflib/bin/vcfleftalign
VCFLIBDECOMP=$(CTX_PATH)/libs/vcflib/bin/vcfallelicprimitives
HAPLEN=$(CTX_PATH)/scripts/longest-haplotype.sh
OLDCLEAN=$(CTX_PATH)/scripts/clean_bubbles.pl

# Measure genome size if not passed
ifndef GENOMESIZE
	GENOMESIZE=$(shell $(SEQCAT) $(SEQ) | tr -d '\n' | wc | grep -o '[0-9]*$$')
endif

# Calculate some numbers
NCHROMS=$(shell bc <<< '$(NUM_INDIVS) * $(PLOIDY)')
NINDIVS_REF=$(shell echo $$(($(NUM_INDIVS) + 1)))

# Generate file names
GENOMES=$(shell echo genomes/genome{1..$(NCHROMS)}.fa)
READS=$(shell echo reads/reads{1..$(NCHROMS)}.{1..2}.fa.gz)
RAWGRAPHS=$(shell echo k$(KMER)/graphs/sample{1..$(NUM_INDIVS)}.raw.ctx)
CLEANGRAPHS=$(RAWGRAPHS:.raw.ctx=.clean.ctx)
REFGENOME=ref/genome0.fa

DOOLD=$(shell if [ $(NUM_INDIVS) -gt 10 ]; then echo 'no'; else echo 'yes'; fi)

ifeq ($(DOOLD),yes)
	SET={oldbc,newbc,se,pe,sepe}
else
	# don't do oldbc if more than 10 individuals
	SET={newbc,se,pe,sepe}
endif

# ref is calling with the ref in the graph
# noref is calling without the ref

# Can't do noref if we only have one chrom
ifeq ($(NCHROMS),1)
	# Always need ref
	PATHS=$(shell echo k$(KMER)/graphs/pop.{se,pe,sepe}.ref.ctp)
	BUBBLES=$(shell echo `eval echo k$(KMER)/bubbles/samples.$(SET).ref.bubbles.gz`)
	TRUTHBUBBLES=k$(KMER)/vcfs/truth.ref.bub.vcf
	BUBBLESCMPRULES=compare-ref-bubbles
	NORMCMPRULES=$(shell echo `eval echo compare-$(SET).ref-norm compare-runcalls-norm`)
else
	ifeq ($(NCHROMS),2)
		TRUTHBUBBLES=k$(KMER)/vcfs/truth.noref.bub.vcf
		BUBBLESCMPRULES=compare-noref-bubbles
	endif

	PATHS=$(shell echo k$(KMER)/graphs/pop.{se,pe,sepe}.{noref,ref}.ctp)
	BUBBLES=$(shell echo `eval echo k$(KMER)/bubbles/samples.$(SET).{noref,ref}.bubbles.gz`)
	NORMCMPRULES=$(shell echo `eval echo compare-$(SET).{noref,ref}-norm compare-runcalls-norm`)
endif

READLISTS=$(shell echo reads/reads{1..$(NUM_INDIVS)}.{1,2}.falist)

MGLIST=$(shell for i in {1..$(NCHROMS)}; do echo -n " genomes/genome$$i.fa genomes/mask$$i.fa"; done)
MGLIST_DECOMP_noref=ref/genome0.fa ref/mask0.clean.fa $(MGLIST)
MGLIST_DECOMP_ref=ref/genome0.fa ref/mask0.fa $(MGLIST)

MGLIST_BUBBLES_noref=$(MGLIST)
MGLIST_BUBBLES_ref=ref/genome0.fa ref/mask0.fa $(MGLIST)

se_list=$(shell for i in `seq 1 $(NUM_INDIVS)`; do \
	j=$$(($$i-1)); echo -n " --col $$j"; \
	for k in `seq $$(($$j * $(PLOIDY) + 1)) $$(($$i * $(PLOIDY)))`; do \
		echo -n " --seq reads/reads$$k.1.fa.gz --seq reads/reads$$k.2.fa.gz"; \
	done; \
done)

pe_list=$(shell for i in `seq 1 $(NUM_INDIVS)`; do \
	j=$$(($$i-1)); echo -n " --col $$j"; \
	for k in `seq $$(($$j * $(PLOIDY) + 1)) $$(($$i * $(PLOIDY)))`; do \
		echo -n " --seq2 reads/reads$$k.1.fa.gz reads/reads$$k.2.fa.gz"; \
	done; \
done)

# sepe_list=$(shell for i in `seq 1 $(NUM_INDIVS)`; do \
# 	j=$$(($$i-1)); echo -n " --col $$j"; \
# 	for k in `seq $$(($$j * $(PLOIDY) + 1)) $$(($$i * $(PLOIDY)))`; do \
# 		echo -n " --p reads/reads$$k.1.fa.gz --seq reads/reads$$k.2.fa.gz"; \
# 		echo -n " --seq2 reads/reads$$k.1.fa.gz reads/reads$$k.2.fa.gz"; \
# 	done; \
# done)

BUBBLEVCFS=$(subst .bubbles.gz,.bub.vcf,$(subst /bubbles/,/vcfs/,$(BUBBLES)))
TRUTHDECOMP=$(TRUTHBUBBLES:.bub.vcf=.decomp.vcf)
TRUTHVCFS=$(TRUTHBUBBLES:.bub.vcf=.norm.vcf)

PLACEVCFS=$(BUBBLEVCFS:.bub.vcf=.decomp.vcf)
PASSVCFS=$(BUBBLEVCFS:.bub.vcf=.pass.vcf)
NORMVCFS=$(BUBBLEVCFS:.bub.vcf=.norm.vcf)

FLANKFILES=$(BUBBLEVCFS:.vcf=.5pflanks.fa.gz)
SAMFILES=$(BUBBLEVCFS:.vcf=.5pflanks.sam)

ifdef ERRPROF
	USECALIB=-p $(ERRPROF)
  GRAPHS_noref=$(CLEANGRAPHS)
else
  GRAPHS_noref=$(RAWGRAPHS)
endif

GRAPHS_ref=$(GRAPHS_noref) ref/ref.k$(KMER).ctx

GENOMES_noref=$(GENOMES)
GENOMES_ref=$(GENOMES) $(REFGENOME)

KEEP=$(GENOMES) $(READS) $(PATHS) $(TRUTHBUBBLES) $(TRUTHDECOMP) $(TRUTHVCFS) $(BUBBLES) $(PASSVCFS) $(PLACEVCFS) $(NORMVCFS)

#
# Phony commands
#

all: repo checkcmds $(KEEP) compare-bubbles compare-normvcf traverse

checkcmds:
	@if [ '$(SEQ)' == '' ]; then echo "You need to specify SEQ=.. Please and thank you."; exit -1; fi;

test:
	@echo KEEP: $(KEEP)

repo:
	@echo Last git commit: `git log -1 --format="%H: %s [%aD]"`

# % is noref or ref
compare-bubbles: $(BUBBLESCMPRULES)
$(BUBBLESCMPRULES): $(BUBBLEVCFS) $(TRUTHBUBBLES)
$(BUBBLESCMPRULES): compare-%-bubbles:
	@echo == Released Cortex $* ==
	$(BIOINF)/sim_mutations/sim_compare.pl k$(KMER)/vcfs/truth.$*.bub.vcf k$(KMER)/vcfs/samples.oldbc.$*.bub.vcf k$(KMER)/vcfs/truth.oldbc.$*.vcf OLDBC k$(KMER)/vcfs/falsepos.oldbc.$*.vcf $(GENOMES_$*)
	$(HAPLEN) k$(KMER)/vcfs/samples.oldbc.$*.bub.vcf
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

compare-normvcf: $(NORMCMPRULES)
$(NORMCMPRULES): $(NORMVCFS) k$(KMER)/vcfs/truth.ref.norm.vcf k$(KMER)/vcfs/truth.noref.norm.vcf
$(NORMCMPRULES): compare-%-norm: k$(KMER)/vcfs/samples.%.norm.vcf
	@echo == $< ==
	r=`echo $* | grep -oE '(no)*ref'`; if [[ $$r == '' ]]; then r='ref'; fi; \
	$(BIOINF)/vcf_scripts/vcf_isec.pl k$(KMER)/vcfs/truth.$$r.norm.vcf $< > /dev/null

# 1..PLOIDY chroms per colour
# e.g. genomes/genome{1,2}.fa loaded into colour 0 when PLOIDY=2
traverse: $(PATHS) k$(KMER)/graphs/pop.ref.ctx
	$(TRAVERSE)                                     --ncontigs $(NCONTIGS) --colour 0 --print k$(KMER)/graphs/pop.ref.ctx | $(BIOINF)/sim_mutations/sim_substrings.pl $(KMER) 0.1 - genomes/genome{1..$(PLOIDY)}.fa
	$(TRAVERSE) -p k$(KMER)/graphs/pop.se.ref.ctp   --ncontigs $(NCONTIGS) --colour 0 --print k$(KMER)/graphs/pop.ref.ctx | $(BIOINF)/sim_mutations/sim_substrings.pl $(KMER) 0.1 - genomes/genome{1..$(PLOIDY)}.fa
	$(TRAVERSE) -p k$(KMER)/graphs/pop.pe.ref.ctp   --ncontigs $(NCONTIGS) --colour 0 --print k$(KMER)/graphs/pop.ref.ctx | $(BIOINF)/sim_mutations/sim_substrings.pl $(KMER) 0.1 - genomes/genome{1..$(PLOIDY)}.fa
	$(TRAVERSE) -p k$(KMER)/graphs/pop.sepe.ref.ctp --ncontigs $(NCONTIGS) --colour 0 --print k$(KMER)/graphs/pop.ref.ctx | $(BIOINF)/sim_mutations/sim_substrings.pl $(KMER) 0.1 - genomes/genome{1..$(PLOIDY)}.fa
	$(CTXSTATS) k$(KMER)/graphs/pop.noref.ctx:0
	@echo == ref copy number ==
	$(CTX_PATH)/bin/ctx31 view --kmers ref/ref.k$(KMER).ctx | awk '{n[$$2]++} END {for (i in n) print i,n[i]}' | sort -n

clean:
	rm -rf ref genomes reads k$(KMER) runcalls gap_sizes.*.csv mp_sizes.*.csv stampy.sh

# .NOTPARALLEL: $(NORMCMPRULES) compare-bubbles compare-normvcf $(NORMCMPRULES)

#
# Patterns
#

ref/stampy.stidx: ref/ref.fa
	$(STAMPY_BIN) -G ref/stampy ref/ref.fa

ref/stampy.sthash: ref/stampy.stidx
	$(STAMPY_BIN) -g ref/stampy -H ref/stampy

ref/ref.fa.fai: ref/ref.fa
	samtools faidx ref/ref.fa

# Generate genomes
ref/mask0.clean.fa $(GENOMES): ref/ref.fa
ref/ref.fa:
	mkdir -p genomes ref
	$(BIOINF)/sim_mutations/sim_mutations.pl --snps $(SNPS) --indels $(INDELS) --invs $(INV) --invlen $(INVLEN) genomes/ $$(($(NCHROMS)+1)) $(SEQ)
	mv genomes/genome0.fa genomes/mask0.fa ref/
	awk 'BEGIN{print">mask";for(i=0;i<$(GENOMESIZE);i++) {printf "."}print""}' > ref/mask0.clean.fa
	cat ref/genome0.fa | tr -d '-' | $(FACAT) -w 50 > ref/ref.fa
	# To generate a mutation free reference genome:
	# cp $(SEQ) ref/ref.fa && cp $(SEQ) ref/genome0.fa && cp ref/mask0.clean.fa ref/mask0.fa

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
	b=$$(($* * $(PLOIDY))); a=$$(($$b-$(PLOIDY)+1)); \
	for i in `seq $$a $$b`; do \
		echo reads$$i.1.fa.gz >> reads/reads$*.1.falist; \
		echo reads$$i.2.fa.gz >> reads/reads$*.2.falist; \
	done

k$(KMER)/graphs/sample%.clean.ctx: k$(KMER)/graphs/sample%.raw.ctx
	$(CLEANCTX) --supernodes --tips 61 --covgs k$(KMER)/graphs/sample$*.contig_covg.csv $@ $<

k$(KMER)/graphs/sample%.raw.ctx: $(READS)
	mkdir -p k$(KMER)/graphs
	b=$$(($* * $(PLOIDY))); a=$$(($$b-$(PLOIDY)+1)); \
	files=$$(for k in `seq $$a $$b`; do echo -n " --seq2 reads/reads$$k.1.fa.gz reads/reads$$k.2.fa.gz"; done); \
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample Sample$* $$files k$(KMER)/graphs/sample$*.raw.ctx;

# Delete .raw.ctx and .clean.ctx paths after running
# .INTERMEDIATE: $(GRAPHS_noref)

k$(KMER)/graphs/pop.noref.ctx: $(GRAPHS_noref)
k$(KMER)/graphs/pop.ref.ctx: $(GRAPHS_ref)
k$(KMER)/graphs/pop.%.ctx:
	$(JOINCTX) --ncols $(NUMCOLS) -m $(MEM) $@ $(GRAPHS_$*)
	$(INFERCTX) $@

# Paths
$(PATHS): k$(KMER)/graphs/pop.noref.ctx k$(KMER)/graphs/pop.ref.ctx

k$(KMER)/graphs/pop.sepe.%.ctp: k$(KMER)/graphs/pop.%.ctx k$(KMER)/graphs/pop.se.%.ctp
	$(THREADCTX) -m $(MEM) -t $(NTHREADS) -p k$(KMER)/graphs/pop.se.$*.ctp $(pe_list) $@ k$(KMER)/graphs/pop.$*.ctx
	for f in *_sizes.*.csv; do mv $$f k$(KMER)/graphs/se.$$f; done

k$(KMER)/graphs/pop.%.noref.ctp: k$(KMER)/graphs/pop.noref.ctx
	$(THREADCTX) -m $(MEM) -t $(NTHREADS) $($*_list) $@ $<
	for f in *_sizes.*.csv; do mv $$f k$(KMER)/graphs/se.$$f; done

k$(KMER)/graphs/pop.%.ref.ctp: k$(KMER)/graphs/pop.ref.ctx ref/ref.fa
	$(THREADCTX) -m $(MEM) -t $(NTHREADS) $($*_list) --col $(NUM_INDIVS) --seq ref/ref.fa $@ $<
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
	callargs=`if [ '$*' == 'ref' ]; then echo '--ref $(NUM_INDIVS)'; fi`; \
	$(CALLCTX) -t $(NTHREADS) -m $(MEM) --maxallele $(MAXALLELE) $$callargs $< $@

# % => {se,pe,sepe}.{ref.noref}
k$(KMER)/bubbles/samples.%.bubbles.gz: k$(KMER)/graphs/pop.%.ctp
	mkdir -p k$(KMER)/bubbles
	r=`echo $@ | grep -oE '(no)?ref'`; \
	callargs=`if [ '$*' == 'ref' ]; then echo '--ref $(NUM_INDIVS)'; fi`; \
	$(CALLCTX) -t $(NTHREADS) -m $(MEM) --maxallele $(MAXALLELE) $$callargs -p $< k$(KMER)/graphs/pop.$$r.ctx $@

k$(KMER)/vcfs/truth.%.bub.vcf: ref/ref.fa $(GENOMES)
	mkdir -p k$(KMER)/vcfs
	$(BIOINF)/sim_mutations/sim_bubble_vcf.pl $(KMER) $(MGLIST_BUBBLES_$*) > $@

k$(KMER)/vcfs/truth.%.decomp.vcf: ref/ref.fa $(GENOMES)
	mkdir -p k$(KMER)/vcfs
	$(BIOINF)/sim_mutations/sim_decomp_vcf.pl $(MGLIST_DECOMP_$*) > k$(KMER)/vcfs/truth.$*.decomp.vcf

k$(KMER)/vcfs/samples.%.bub.vcf k$(KMER)/vcfs/samples.%.bub.5pflanks.fa.gz: k$(KMER)/bubbles/samples.%.bubbles.gz
	mkdir -p k$(KMER)/vcfs
	$(PROCCTX) k$(KMER)/bubbles/samples.$*.bubbles.gz k$(KMER)/vcfs/samples.$*.bub
	gzip -d -f k$(KMER)/vcfs/samples.$*.bub.vcf.gz

$(SAMFILES): ref/stampy.stidx ref/stampy.sthash
k$(KMER)/vcfs/samples.%.bub.5pflanks.sam: k$(KMER)/vcfs/samples.%.bub.5pflanks.fa.gz
	mkdir -p ref
	$(STAMPY_BIN) -g ref/stampy -h ref/stampy $(MAPARGS) --inputformat=fasta -M $< > $@

k$(KMER)/vcfs/samples.%.decomp.vcf: k$(KMER)/vcfs/samples.%.bub.vcf k$(KMER)/vcfs/samples.%.bub.5pflanks.sam ref/ref.fa
	$(PLACECTX) --minmapq $(MINMAPQ) k$(KMER)/vcfs/samples.$*.bub.vcf k$(KMER)/vcfs/samples.$*.bub.5pflanks.sam ref/ref.fa > $@

k$(KMER)/vcfs/samples.%.pass.vcf: k$(KMER)/vcfs/samples.%.decomp.vcf
	cat $< | awk '$$1 ~ /^#/ || $$7 ~ /^(PASS|\.)$$/' | vcf-sort > $@

reads/reads.index: $(READLISTS)
	for i in {1..$(NUM_INDIVS)}; do echo -e Sample$$i"\t"."\t"reads/reads$$i.1.falist"\t"reads/reads$$i.2.falist; done > reads/reads.index

REFPATH=$(realpath ref/ref.fa)

ref/ref.falist: ref/ref.fa
	echo $(REFPATH) > ref/ref.falist

ref/ref.k$(KMER).ctx: ref/ref.fa
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample ref --seq ref/ref.fa ref/ref.k$(KMER).ctx

# Left align with vcflib or bcftools:
# LEFTALIGN=$(VCFLIBALIGN) --reference ref/ref.fa
LEFTALIGN=$(BCFTOOLS) norm --remove-duplicate -f ref/ref.fa -

k$(KMER)/vcfs/samples.runcalls.norm.vcf: ref/ref.fa.fai reads/reads.index ref/ref.falist ref/ref.k$(KMER).ctx $(CORTEX_PATH)/bin/cortex_var_31_c2 $(CORTEX_PATH)/bin/cortex_var_31_c$(NINDIVS_REF)
	@echo == Run Calls ==
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
						  --max_var_len $(MAXALLELE) \
						  --workflow independent --logfile runcalls/runcalls.log
						  # --apply_pop_classifier
	# 1) Add KMER and contig lines
  # 2) Remove ref mismatches
	# 3) Left align indels
  # 4) Remove duplicates
	$(BIOINF)/vcf_scripts/vcf_header.pl --entries \
	  runcalls/vcfs/samples.runcalls_union_BC_calls_k31.decomp.vcf \
	  '+tag:INFO,KMER,.,Integer,"Kmer called at"' \
	  '+meta:contig=<ID=ref,length=1000>' | \
	$(BIOINF)/vcf_scripts/vcf_header_add_contigs.pl - ref/ref.fa | \
	$(BIOINF)/vcf_scripts/vcf_filter_by_ref.pl - ref/ref.fa | \
	$(VCFLIBDECOMP) | $(LEFTALIGN) | \
	$(BIOINF)/vcf_scripts/vcf_remove_dupes.pl > $@

$(NORMVCFS): ref/ref.fa.fai

$(TRUTHVCFS): $(TRUTHBUBBLES) ref/ref.fa.fai

# % is ref or noref
k$(KMER)/vcfs/truth.%.norm.vcf: k$(KMER)/vcfs/truth.%.decomp.vcf
	cat $< | \
	$(VCFLIBDECOMP) | $(LEFTALIGN) | \
	$(BIOINF)/vcf_scripts/vcf_remove_dupes.pl > $@

k$(KMER)/vcfs/samples.%.norm.vcf: k$(KMER)/vcfs/samples.%.pass.vcf
	cat $< | \
	$(VCFLIBDECOMP) | $(LEFTALIGN) | \
	$(BIOINF)/vcf_scripts/vcf_remove_dupes.pl > k$(KMER)/vcfs/samples.$*.norm.vcf

.PHONY: all clean test repo checkcmds
.PHONY: compare-bubbles compare-normvcf $(NORMCMPRULES)
.PHONY: traverse
.FORCE: repo
