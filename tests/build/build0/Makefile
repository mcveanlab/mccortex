SHELL:=/bin/bash -euo pipefail

#
# Generate 60 bases of sequence, build graph, save, sort and index it, then
# re-assemble the sequence. We should get back a single contig of length 60
#

K=21
CTXDIR=../../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CONTIGSTATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl

TGTS=seq.fa seq.k$(K).ctx sort.k$(K).ctx sort.k$(K).ctx.idx

all: $(TGTS) test_assemble

clean:
	rm -rf $(TGTS) contig.stats.txt

seq.fa:
	$(DNACAT) -F -n 60 > $@

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -q -m 1M -k $(K) \
	                  --sample Wallace \
	                  --sample Gromit --seq seq.fa \
	                  --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(MCCORTEX) check -q $@
	$(MCCORTEX) view -q $@
	$(MCCORTEX) view -q -k 1,3,5:$@:2,1,0 > /dev/null

sort.k$(K).ctx: seq.k$(K).ctx
	cp $< $@
	$(MCCORTEX) view -q -k $< > /dev/null
	$(MCCORTEX) sort -q $@
	$(MCCORTEX) check -q $@
	$(MCCORTEX) view -k $@ > /dev/null

%.ctx.idx: %.ctx
	$(MCCORTEX) index -q $< > $@
	cat $@

contig.stats.txt: sort.k$(K).ctx
	$(MCCORTEX) contigs -q -o - $<:1 | \
	  $(MCCORTEX) rmsubstr -q -n 1M -k $(K) - | \
	  $(CONTIGSTATS) - > $@

test_assemble: contig.stats.txt
	grep -q 'contigs: .* 1$$' $<
	grep -q 'length: .* 60$$' $<

.PHONY: all clean test_assemble
