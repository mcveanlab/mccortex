SHELL:=/bin/bash -euo pipefail

TARGETS=rnd.fa rnd.k39.ctx rnd.clean.k39.ctx

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex63

all: $(TARGETS)

clean:
	rm -rf $(TARGETS)

rnd.fa:
	$(DNACAT) -F -n 200 > $@

rnd.k39.ctx: rnd.fa
	$(MCCORTEX) build -k 39 --sample rnd --seq rnd.fa $@

rnd.clean.k39.ctx: rnd.k39.ctx
	$(MCCORTEX) clean -o $@ --unitigs=2 $<
	$(MCCORTEX) check $@

.PHONY: all clean
