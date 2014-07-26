SHELL:=/bin/bash

TARGETS=rnd.fa rnd.k39.ctx rnd.clean.k39.ctx

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/ctx63

all: $(TARGETS)

clean:
	rm -rf $(TARGETS)

rnd.fa:
	$(DNACAT) -F -n 200 > $@

rnd.k39.ctx: rnd.fa
	$(CTX) build -k 39 --sample rnd --seq rnd.fa $@

rnd.clean.k39.ctx: rnd.k39.ctx
	$(CTX) clean -o $@ $<
	$(CTX) check $@

.PHONY: all clean
