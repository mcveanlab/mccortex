# Calling variants when the alleles contain duplicated kmers

SHELL:=/bin/bash

K=9
CTX=../../bin/ctx31
STAMPY=python2.6 ~/bioinf/stampy-1.0.23/stampy.py

SEQ=seq1.fa seq2.fa
GRAPHS=seq.k$(K).ctx
PATHS=seq.k$(K).ctp
PLOTS=seq.k$(K).pdf seq.k$(K).kmers.pdf
REF=ref.fa ref.stidx ref.sthash

CALLS0=calls.k$(K).bubbles.gz
CALLS1=calls.k$(K).5pflanks.fa.gz calls.k$(K).vcf.gz calls.k$(K).5pflanks.sam
CALLS2=calls.k$(K).placed.vcf

CALLS=$(CALLS0) $(CALLS1) $(CALLS2)

KEEP=$(SEQ) $(REF) $(GRAPHS) $(PATHS) $(PLOTS) $(CALLS)

all: $(KEEP)

clean:
	rm -rf $(KEEP)
	rm -rf gap_sizes.*.csv mp_sizes.*.csv

plots: $(PLOTS)

# repeat of CCTCGTCTCATACGCTGCACTTACG
# GATATAGCTCCGCACATCGAGCCTGC-CCTCGTCTCATACGCTGCACTTACGCCTCGTCTCATACGCTGCACTTACG-AACGAGGAGAGATTAAGCCGCCCAG
# GATATAGCTCCGCACATCGAGCCTGC-GGA-AACGAGGAGAGATTAAGCCGCCCAG

# GATATAGCTCCGCACATCGAGCCTGCGGAAACGAGGAGAGATTAAGCCGCCCAG
# GATATAGCTCCGCACATCGAGCCTGCCCTCGTCTCATACGCTGCACTTACGCCTCGTCTCATACGCTGCACTTACGAACGAGGAGAGATTAAGCCGCCCAG

seq1.fa:
	echo -e ">seq1\nGATATAGCTCCGCACATCGAGCCTGCCCTCGTCTCATACGCTGCACTTACGCCTCGTCTCATACGCTGCACTTACGAACGAGGAGAGATTAAGCCGCCCAG" > seq1.fa

seq2.fa:
	echo -e ">seq2\nGATATAGCTCCGCACATCGAGCCTGCGGAAACGAGGAGAGATTAAGCCGCCCAG" > seq2.fa

ref.fa: seq1.fa seq2.fa
	cp seq2.fa ref.fa

ref.stidx: ref.fa
	$(STAMPY) -G ref ref.fa

ref.sthash: ref.stidx
	$(STAMPY) -g ref -H ref

seq.k$(K).ctx: seq1.fa seq2.fa
	$(CTX) build -k $(K) --sample seq1 --seq seq1.fa --sample seq2 --seq seq2.fa seq.k$(K).ctx

seq.k$(K).ctp: seq1.fa seq2.fa
	$(CTX) thread --seq seq1.fa -o seq.0.k$(K).ctp seq.k$(K).ctx:0
	$(CTX) thread --seq seq2.fa -o seq.1.k$(K).ctp seq.k$(K).ctx:1
	$(CTX) pjoin -o $@ seq.0.k$(K).ctp seq.1.k$(K).ctp
	rm seq.0.k$(K).ctp seq.1.k$(K).ctp

seq.k$(K).pdf: seq.k$(K).ctx
	 $(CTX) unitigs --graphviz seq.k$(K).ctx | dot -Tpdf > seq.k$(K).pdf

seq.k$(K).kmers.pdf: seq.k$(K).ctx
	../../scripts/cortex_to_graphviz.pl seq.k$(K).ctx | dot -Tpdf > seq.k$(K).kmers.pdf

calls.k$(K).bubbles.gz:
	$(CTX) bubbles -p seq.k$(K).ctp -o calls.k$(K).bubbles.gz seq.k$(K).ctx

calls.k$(K).5pflanks.fa.gz calls.k$(K).vcf.gz:
	$(CTX) unique calls.k$(K).bubbles.gz calls.k$(K)

calls.k$(K).5pflanks.sam: calls.k$(K).5pflanks.fa.gz ref.stidx ref.sthash
	$(STAMPY) -g ref -h ref --inputformat=fasta -M calls.k$(K).5pflanks.fa.gz > $@

calls.k$(K).placed.vcf: calls.k$(K).vcf.gz calls.k$(K).5pflanks.sam
	$(CTX) place --out calls.k$(K).placed.vcf calls.k$(K).vcf.gz calls.k$(K).5pflanks.sam ref.fa

.PHONY: all clean plot
