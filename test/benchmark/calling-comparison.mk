
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

# USEREF=1
# DEV: add ref to calling

SHELL := /bin/bash

CORTEX_PATH=$(HOME)/cortex/releases/CORTEX_release_v1.0.5.20
SHADE_PATH=$(HOME)/cortex/versions/current
CTX_PATH=$(HOME)/cortex/versions/compact_hash

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

ifdef USEREF
	NUMCOLS=$(shell echo $$(($(NUM_INDIVS)+1)))
else
	NUMCOLS=$(NUM_INDIVS)
endif

RELEASECTX=$(CORTEX_PATH)/bin/cortex_var_31_c$(NUMCOLS) --kmer_size $(KMER) --mem_height $(MEMHEIGHT) --mem_width $(MEMWIDTH)
SHADECTX=$(SHADE_PATH)/bin/cortex_var_31_c$(NUMCOLS)_s8 --kmer_size $(KMER) --mem_height $(MEMHEIGHT) --mem_width $(MEMWIDTH)
BUILDCTX=$(CTX_PATH)/bin/ctx31 build
CLEANCTX=$(CTX_PATH)/bin/ctx31 clean
JOINCTX=$(CTX_PATH)/bin/ctx31 join
INFERCTX=$(CTX_PATH)/bin/ctx31 inferedges --pop
THREADCTX=$(CTX_PATH)/bin/ctx31 thread
CALLCTX=$(CTX_PATH)/bin/ctx31 call
PROCCTX=$(CTX_PATH)/bin/ctx31 unique
PLACECTX=$(CTX_PATH)/bin/ctx31 place
TRAVERSE=$(CTX_PATH)/bin/traversal31

RUNCALLS=$(CORTEX_PATH)/scripts/calling/run_calls.pl

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

MEM=$(shell bc <<< '(($(MEMWIDTH) * 2^$(MEMHEIGHT)) * (8+1+4)*8+1) / 8')

SEQPATH=$(realpath $(SEQ))

# Generate file names
GENOMES=$(shell echo genomes/genome{0..$(LASTCHROM)}.fa)
READS=$(shell echo reads/reads{0..$(LASTCHROM)}.{1..2}.fa.gz)
RAWGRAPHS=$(shell echo graphs/sample{0..$(LASTINDIV)}.raw.ctx)
CLEANGRAPHS=$(RAWGRAPHS:.raw.ctx=.clean.ctx)
PATHS=$(shell echo graphs/pop.{se,pe,sepe}.ctp)
BUBBLES=$(shell echo bubbles/samples.{oldbc,newbc,shaded,se,pe,sepe}.bubbles.gz)

READ1LISTS=$(shell echo reads/reads{0..$(LASTINDIV)}.1.falist)
READ2LISTS=$(shell echo reads/reads{0..$(LASTINDIV)}.2.falist)
SHADEDLIST=$(shell for i in {0..$(LASTINDIV)}; do echo -n " --pe_list reads/reads$$i.1.falist,reads/reads$$i.2.falist"; done)
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

BUBBLEVCFS=$(shell echo vcfs/samples.{oldbc,newbc,shaded,se,pe,sepe}.bub.vcf)
PLACEVCFS=$(BUBBLEVCFS:.bub.vcf=.decomp.vcf)
PASSVCFS=$(BUBBLEVCFS:.bub.vcf=.pass.vcf)
NORMVCFS=$(BUBBLEVCFS:.bub.vcf=.norm.vcf)

FLANKFILES=$(shell echo vcfs/samples.{oldbc,newbc,shaded,se,pe,sepe}.bub.5pflanks.fa.gz)
SAMFILES=$(FLANKFILES:.fa.gz=.sam)

ifdef ERRPROF
	SAMPLEGRAPHS=$(CLEANGRAPHS)
else
	SAMPLEGRAPHS=$(RAWGRAPHS)
endif

ifdef USEREF
	GRAPHS=$(SAMPLEGRAPHS) ref/ref.k$(KMER).ctx
	REFARGS=--ref $(NUM_INDIVS)
	ALLCHROMS=$(GENOMES) $(SEQ)
else
	GRAPHS=$(SAMPLEGRAPHS)
	ALLCHROMS=$(GENOMES)
endif

NORMCMPRULES=$(shell echo compare-{oldbc,newbc,se,pe,sepe,runcalls}-norm)

KEEP=$(GENOMES) $(READS) $(SAMPLEGRAPHS) $(PATHS) $(BUBBLES) $(PASSVCFS) $(PLACEVCFS) $(NORMVCFS)

# Directories
DIRS=ref genomes reads graphs bubbles vcfs runcalls

all: check-cmds $(DIRS) $(KEEP) compare-bubbles compare-normvcf traverse

check-cmds:
	@if [ '$(SEQ)' == '' ]; then echo "You need to specify SEQ=.. Please and thank you."; exit -1; fi;

test:
	echo $(NUMCOLS)

$(READS): $(GENOMES)

compare-bubbles: $(BUBBLEVCFS) vcfs/truth.bub.vcf
	@echo == Released Cortex ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.bub.vcf vcfs/samples.oldbc.bub.vcf vcfs/truth.oldbc.vcf OLDBC vcfs/falsepos.oldbc.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.oldbc.bub.vcf
	@echo == New Bubble Caller ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.oldbc.vcf vcfs/samples.newbc.bub.vcf vcfs/truth.newbc.vcf NEWBC vcfs/falsepos.newbc.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.newbc.bub.vcf
	@echo == Shaded Bubble Caller ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.newbc.vcf vcfs/samples.shaded.bub.vcf vcfs/truth.shaded.vcf SHADED vcfs/falsepos.shaded.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.shaded.bub.vcf
	@echo == Paths se ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.shaded.vcf vcfs/samples.se.bub.vcf vcfs/truth.se.vcf PAC vcfs/falsepos.se.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.se.bub.vcf
	@echo == Paths pe ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.se.vcf vcfs/samples.pe.bub.vcf vcfs/truth.pe.vcf PAC vcfs/falsepos.pe.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.pe.bub.vcf
	@echo == Paths sepe ==
	$(BIOINF)/sim_mutations/sim_compare.pl vcfs/truth.pe.vcf vcfs/samples.sepe.bub.vcf vcfs/truth.sepe.vcf PAC vcfs/falsepos.sepe.vcf $(ALLCHROMS)
	$(HAPLEN) vcfs/samples.sepe.bub.vcf
	@echo == Truth ==
	$(HAPLEN) vcfs/truth.bub.vcf

compare-normvcf: $(NORMCMPRULES)
$(NORMCMPRULES): $(NORMVCFS) vcfs/truth.norm.vcf
$(NORMCMPRULES): compare-%-norm: vcfs/samples.%.norm.vcf
	@echo == $< ==
	$(BIOINF)/vcf_scripts/vcf_isec.pl vcfs/truth.norm.vcf $< > /dev/null

traverse: $(PATHS) graphs/pop.ctx
	$(TRAVERSE)                        --nsamples 10000 graphs/pop.ctx
	$(TRAVERSE) -p graphs/pop.se.ctp   --nsamples 10000 graphs/pop.ctx
	$(TRAVERSE) -p graphs/pop.pe.ctp   --nsamples 10000 graphs/pop.ctx
	$(TRAVERSE) -p graphs/pop.sepe.ctp --nsamples 10000 graphs/pop.ctx

ref/stampy.stidx:
	$(STAMPY_BIN) -G ref/stampy $(SEQ)

ref/stampy.sthash:
	$(STAMPY_BIN) -g ref/stampy -H ref/stampy

$(SEQ).fai:
	samtools faidx $(SEQ)

$(DIRS):
	mkdir -p $(DIRS)

clean:
	rm -rf $(DIRS) gap_sizes.*.csv mp_sizes.*.csv stampy.sh

#
# Patterns
#

$(GENOMES):
	zcat -f $(SEQ) | $(BIOINF)/sim_mutations/sim_mutations.pl --snps $(SNPS) --indels $(INDELS) --invs $(INV) --invlen $(INVLEN) genomes/ $(NCHROMS) -

# Reads
reads/reads%.1.fa.gz reads/reads%.2.fa.gz: genomes/genome%.fa
	cat genomes/genome$*.fa | tr -d '-' | $(READSIM) -r - -i $(MPSIZE) -v 0.2 -l $(READLEN) -d $(ALLELECOVG) $(USECALIB) reads/reads$*

reads/reads%.1.falist reads/reads%.2.falist:
	echo reads$$(($* * 2)).1.fa.gz >> reads/reads$*.1.falist
	echo reads$$(($* * 2+1)).1.fa.gz >> reads/reads$*.1.falist
	echo reads$$(($* * 2)).2.fa.gz >> reads/reads$*.2.falist
	echo reads$$(($* * 2+1)).2.fa.gz >> reads/reads$*.2.falist

$(SAMPLEGRAPHS): $(READS)

$(CLEANGRAPHS): $(RAWGRAPHS)

graphs/sample%.clean.ctx: graphs/sample%.raw.ctx
	$(CLEANCTX) $@ $<

graphs/sample%.raw.ctx:
	a=$$(($* * $(PLOIDY))); b=$$(($$first+$(PLOIDY)-1)); \
	files=$$(for k in `seq $$a $$b`; do echo -n " --seq2 reads/reads$$k.1.fa.gz reads/reads$$k.2.fa.gz"; done); \
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample Sample$* $$files graphs/sample$*.raw.ctx;

graphs/pop.ctx: $(GRAPHS)
	$(JOINCTX) -m $(MEM) graphs/pop.ctx $(GRAPHS)
	$(INFERCTX) graphs/pop.ctx

# Paths
graphs/pop.%.ctp: graphs/pop.ctx
	$(THREADCTX) -t 1 $($*_list) $(NUMCOLS) $@ graphs/pop.ctx graphs/pop.ctx
	for f in *_sizes.*.csv; do mv $$f graphs/se.$$f; done

# Bubbles
$(BUBBLES): graphs/pop.ctx

bubbles/samples.oldbc.bubbles.gz: graphs/pop.ctx
	time $(RELEASECTX) --multicolour_bin $< --detect_bubbles1 -1/-1 --output_bubbles1 bubbles/samples.oldbc.bubbles --print_colour_coverages
	mv bubbles/samples.oldbc.bubbles bubbles/samples.oldbc.bubbles.dirty
	$(OLDCLEAN) $(KMER) bubbles/samples.oldbc.bubbles.dirty | gzip -c > $@

bubbles/samples.newbc.bubbles.gz: graphs/pop.ctx
	$(CALLCTX) -t 1 $(REFARGS) $< $@

bubbles/samples.shaded.bubbles.gz: $(READ1LISTS) $(READ2LISTS)
bubbles/samples.shaded.bubbles.gz: graphs/pop.ctx
	time $(SHADECTX) --load_binary $< --add_shades $(SHADEDLIST) --paths_caller $@ --paths_caller_cols -1

bubbles/samples.%.bubbles.gz: graphs/pop.%.ctp
	$(CALLCTX) -t 1 -m $(MEM) $(REFARGS) -p $< graphs/pop.ctx $@

vcfs/truth.bub.vcf: $(GENOMES)
	$(BIOINF)/sim_mutations/sim_vcf.pl $(KMER) $(MGLIST) > vcfs/truth.bub.vcf

vcfs/truth.decomp.vcf: $(SEQ) $(GENOMES)
	zcat -f $(SEQ) | $(BIOINF)/sim_mutations/sim_decomp_vcf.pl - $(GENOMES) > vcfs/truth.decomp.vcf

vcfs/samples.%.bub.vcf vcfs/samples.%.bub.5pflanks.fa.gz: bubbles/samples.%.bubbles.gz
	$(PROCCTX) bubbles/samples.$*.bubbles.gz vcfs/samples.$*.bub
	gzip -d -f vcfs/samples.$*.bub.vcf.gz

$(SAMFILES): ref/stampy.stidx ref/stampy.sthash
vcfs/samples.%.bub.5pflanks.sam: vcfs/samples.%.bub.5pflanks.fa.gz
	$(STAMPY_BIN) -g ref/stampy -h ref/stampy --inputformat=fasta -M $< > $@

vcfs/samples.%.decomp.vcf: vcfs/samples.%.bub.vcf vcfs/samples.%.bub.5pflanks.sam
	$(PLACECTX) vcfs/samples.$*.bub.vcf vcfs/samples.$*.bub.5pflanks.sam $(SEQ) > $@

vcfs/samples.%.pass.vcf: vcfs/samples.%.decomp.vcf
	cat $< | awk '$$1 ~ /^#/ || $$7 ~ /^(PASS|\.)$$/' | vcf-sort > $@

reads/reads.index:
	for i in {0..$(LASTINDIV)}; do echo -e Sample$$i"\t"."\t"reads/reads$$i.1.falist"\t"reads/reads$$i.1.falist; done > reads/reads.index

ref/ref.falist:
	echo $(SEQPATH) > ref/ref.falist

ref/ref.k$(KMER).ctx:
	$(BUILDCTX) -k $(KMER) -m $(MEM) --sample ref --seq $(SEQ) ref/ref.k$(KMER).ctx

vcfs/samples.runcalls.norm.vcf runcalls.log: reads/reads.index ref/ref.falist ref/ref.k$(KMER).ctx $(CORTEX_PATH)/bin/cortex_var_31_c2 $(CORTEX_PATH)/bin/cortex_var_31_c$(NINDIVS_REF)
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
	cp runcalls/vcfs/samples.runcalls_union_BC_calls_k31.decomp.vcf vcfs/samples.runcalls.norm.vcf

vcfs/truth.norm.vcf: vcfs/truth.decomp.vcf
	$(BCFTOOLS) norm --remove-duplicate -f $(SEQ) vcfs/truth.decomp.vcf > vcfs/truth.norm.vcf

$(NORMVCFS): $(SEQ).fai

vcfs/samples.%.norm.vcf: vcfs/samples.%.pass.vcf
	$(BCFTOOLS) norm --remove-duplicate -f $(SEQ) $< > vcfs/samples.$*.norm.vcf


.PHONY: all clean compare-bubbles compare-normvcf $(NORMCMPRULES) traverse
