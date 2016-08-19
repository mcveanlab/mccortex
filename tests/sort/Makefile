SHELL:=/bin/bash -euo pipefail

# Any kmer (K=) should work

K=51
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

GRAPHS=seq.fa graph.k$(K).ctx build.then.sort.k$(K).ctx build.and.sort.k$(K).ctx
MISC=kmers.sorted.k$(K).txt build.then.sort.k$(K).ctx.idx
LOGS=$(addsuffix .log,$(GRAPHS) $(MISC))

all: title $(GRAPHS) $(MISC) check

title:
	@echo "-- Testing sort k=$(K) --"

clean:
	rm -rf $(GRAPHS) $(MISC) $(LOGS)

seq.fa:
	$(DNACAT) -F -n 100 > $@

graph.k$(K).ctx: seq.fa
	$(MCCORTEX) build -k $(K) --sample Jimmy --seq $< $@ >& $@.log
	$(MCCORTEX) check -q $@

build.then.sort.k$(K).ctx: graph.k$(K).ctx
	$(MCCORTEX) sort -o $@ $< >& $@.log
	$(MCCORTEX) check -q $@

build.and.sort.k$(K).ctx: seq.fa
	$(MCCORTEX) build -k $(K) --sort --sample Jimmy --seq $< $@ >& $@.log
	$(MCCORTEX) check -q $@

%.ctx.idx: %.ctx
	$(MCCORTEX) index --out $@ --block-kmers 11 $< >& $@.log

kmers.sorted.k$(K).txt: graph.k$(K).ctx
	$(MCCORTEX) view -q --kmers $< | sort > $@

check: kmers.sorted.k$(K).txt build.then.sort.k$(K).ctx build.and.sort.k$(K).ctx
	diff -q $< <($(MCCORTEX) view -q -k build.then.sort.k$(K).ctx)
	diff -q $< <($(MCCORTEX) view -q -k build.and.sort.k$(K).ctx)

.PHONY: all clean check title
