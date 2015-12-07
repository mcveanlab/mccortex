SHELL=/bin/bash -euo pipefail

# build1: test --intersect and --graph arguments 

CTXDIR=../../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex31
CONTIGSTATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl
K=11

SEQA=TCCCTGGTATCAACTTGTCTCTGGTCGGCC
SEQB=TACCAATCGGGACAACACGGGTCACTACGA
SEQA_SUB=GGTATCAACTTGTCTC
SEQB_SUB=CGGGACAACACG
EXTRA_KMERS=GTTATTAGAATGTTTAATTAATCATAGAAGCAACCTGGGGACACGTCCTCGTGACTGGAGCAGTACAACCCATCAATTAACTCTATTACTATACGGGAAC

all: merge.k$(K).ctx
	@echo "Testing we get exactly 8 kmers..."
	@[[ `$(MCCORTEX) view -q -k $< | awk '{x+=1}END{print x}'` -eq 8 ]] || \
	  (echo "Merged graph does not contain 8 kmers" && false)
	@echo "All looks good."

clean:
	rm -rf {seq,isec0,isec1,small,extra}.fa
	rm -rf isec{0,1}.k$(K).ctx small.k$(K).ctx merge.k$(K).ctx

extra.fa:
	echo $(EXTRA_KMERS) | $(DNACAT) -F - > $@

seq.fa: extra.fa
	echo $(SEQA_SUB) | $(DNACAT) -F - $< > $@

small.fa: extra.fa
	echo $(SEQB_SUB) | $(DNACAT) -F - $< > $@

isec0.fa:
	echo $(SEQA) | $(DNACAT) -F -n 9 - > $@

isec1.fa:
	echo $(SEQB) | $(DNACAT) -F -n 9 - > $@

merge.k$(K).ctx: seq.fa isec0.k$(K).ctx isec1.k$(K).ctx small.k$(K).ctx
	$(MCCORTEX) build -q -m 1M -k $(K) \
	                  --intersect isec1.k$(K).ctx \
	                  --intersect isec0.k$(K).ctx \
	                  --graph 1:small.k$(K).ctx \
	                  --sample Spiderman --seq seq.fa $@
	$(MCCORTEX) check -q $@


isec0.k$(K).ctx: isec0.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample isec0 --seq $< $@

isec1.k$(K).ctx: isec1.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample isec1 --seq $< $@

small.k$(K).ctx: small.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample small --seq $< $@

.PHONY: all clean
