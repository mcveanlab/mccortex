SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex63
K=51

GRAPHS=seq.fa graph.k$(K).ctx build.then.sort.k$(K).ctx build.and.sort.k$(K).ctx
MISC=kmers.sorted.txt build.then.sort.k$(K).ctx.idx
LOGS=$(addsuffix .log,$(GRAPHS)) $(addsuffix .log,$(MISC))

all: $(GRAPHS) $(MISC) check

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

kmers.sorted.txt: graph.k$(K).ctx
	$(MCCORTEX) view -q --kmers $< | sort > $@

check: kmers.sorted.txt build.then.sort.k$(K).ctx build.and.sort.k$(K).ctx
	diff -q kmers.sorted.txt <($(MCCORTEX) view -q -k build.then.sort.k$(K).ctx)
	diff -q kmers.sorted.txt <($(MCCORTEX) view -q -k build.and.sort.k$(K).ctx)

.PHONY: all clean
