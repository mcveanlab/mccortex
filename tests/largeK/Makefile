SHELL:=/bin/bash -euo pipefail

# Run build and clean commands followed by the check command
# to test we can build valid graphs for random sequence
# Works for any K value (e.g. K=11, K=39, K=61)

K=39
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

GRAPHS=rnd.k$(K).ctx rnd.clean.k$(K).ctx
LOGS=$(addsuffix .log,$(GRAPHS))
TGTS=rnd.fa $(GRAPHS)

all: $(TGTS)

clean:
	rm -rf $(TGTS) $(LOGS)

rnd.fa:
	$(DNACAT) -F -n 200 > $@

rnd.k$(K).ctx: rnd.fa
	$(MCCORTEX) build -k $(K) --sample rnd --seq rnd.fa $@ >& $@.log
	$(MCCORTEX) check -q $@

rnd.clean.k$(K).ctx: rnd.k$(K).ctx
	$(MCCORTEX) clean -o $@ --unitigs=2 $< >& $@.log
	$(MCCORTEX) check -q $@

.PHONY: all clean
