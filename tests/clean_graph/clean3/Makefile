SHELL:=/bin/bash -euo pipefail

#
# Three sample cleaning, one colour in the graph.
# Compares output from 'mccortex build --sort', 'mccortex sort' and
# 'mccortex clean --sort'.
# Works with kmers 11..39 inclusive.
#

K=31
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
VIEWKMERS=$(MCCORTEX) view -qk - | cut -d' ' -f 1

INPUT_FASTA=$(shell echo sample.{1,2,3}.fa)
TRUTH_FASTA=$(shell echo truth.{1,2,3}.fa)
SAMPLE_GRAPHS=$(shell echo sample.{1,2,3}.k$(K).raw.ctx)

TGTS=$(INPUT_FASTA) $(TRUTH_FASTA) $(SAMPLE_GRAPHS) pop.k$(K).clean.ctx truth.k$(K).ctx

all: $(TGTS) check

# Each sample has a private mutation that should be cleaned off
# Each sample also has a region of private sequence that will be kept
#                  *                                       *                                          *
sample.1.fa:
	echo CAGGGTCGGCTCgCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCGCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGAAGAGG > $@
	echo TTGAACGTTTTTTACTGGTTACTCTCTACATGATTGGTTAAGTTTACCTC >> $@
	echo TTGAACGTTTTTTACTGGTTACTCTCTACATGATTGGTTAAGTTTACCTC >> $@
sample.2.fa:
	echo CAGGGTCGGCTCCCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCtCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGAAGAGG > $@
	echo AGGCTGTGTTATTGTCCAGTCATAGCATAACTGATTAATACTTTAAGTCT >> $@
	echo AGGCTGTGTTATTGTCCAGTCATAGCATAACTGATTAATACTTTAAGTCT >> $@
sample.3.fa:
	echo CAGGGTCGGCTCCCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCGCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGcAGAGG > $@
	echo CAGGACTACTCAAGTGAGGTATCGGGCCCAAGACACTGCTATTTAGCGTT >> $@
	echo CAGGACTACTCAAGTGAGGTATCGGGCCCAAGACACTGCTATTTAGCGTT >> $@

truth.1.fa:
	echo CAGGGTCGGCTCCCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCGCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGAAGAGG > $@
	echo TTGAACGTTTTTTACTGGTTACTCTCTACATGATTGGTTAAGTTTACCTC >> $@

truth.2.fa:
	echo CAGGGTCGGCTCCCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCGCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGAAGAGG > $@
	echo AGGCTGTGTTATTGTCCAGTCATAGCATAACTGATTAATACTTTAAGTCT >> $@

truth.3.fa:
	echo CAGGGTCGGCTCCCCCGAATCGCAGGGGTCAGCGTTAGACAAGAAACAAGCCGCTTGCGTCACAGTGGGACGAAACCTATTTTCAGTCTCGGTGAAGAGG > $@
	echo CAGGACTACTCAAGTGAGGTATCGGGCCCAAGACACTGCTATTTAGCGTT >> $@

sample.%.k$(K).raw.ctx: sample.%.fa
	$(MCCORTEX) build -q -k $(K) --sort --sample sample$* --seq $< $@

pop.k$(K).clean.ctx: $(SAMPLE_GRAPHS)
	$(MCCORTEX) clean -q -m 32K --sort --ncols 1 --unitigs=2 --tips=0 -o $@ $(SAMPLE_GRAPHS)

truth.k$(K).ctx:
	$(MCCORTEX) build -q -k $(K) --sample sample1 --seq truth.1.fa \
															 --sample sample2 --seq truth.2.fa \
															 --sample sample3 --seq truth.3.fa $@
	$(MCCORTEX) sort -q $@

check: pop.k$(K).clean.ctx truth.k$(K).ctx
	diff -q <($(MCCORTEX) view -qk pop.k$(K).clean.ctx | cut -d' ' -f 1) <($(MCCORTEX) view -qk truth.k$(K).ctx | cut -d' ' -f 1)
	# check each sample is properly sorted
	diff -q <(cat sample.1.k$(K).raw.ctx | $(VIEWKMERS)) <(cat sample.1.k$(K).raw.ctx | $(VIEWKMERS) | sort)
	diff -q <(cat sample.2.k$(K).raw.ctx | $(VIEWKMERS)) <(cat sample.2.k$(K).raw.ctx | $(VIEWKMERS) | sort)
	diff -q <(cat sample.3.k$(K).raw.ctx | $(VIEWKMERS)) <(cat sample.3.k$(K).raw.ctx | $(VIEWKMERS) | sort)

clean:
	rm -rf $(TGTS)

.PHONY: all clean check
