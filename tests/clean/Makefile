SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31
MKDOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
PLOTCOVG=$(CTXDIR)/scripts/plot-covg-hist.R
PLOTLENS=$(CTXDIR)/scripts/plot-length-hist.R

SEQ=seq.fa
GRAPHS=seq.k9.raw.ctx seq.k9.clean.ctx
STATS=lens.before.csv lens.after.csv covgs.before.csv covgs.after.csv
KEEP=$(SEQ) $(GRAPHS) $(STATS)
PLOTS=$(GRAPHS:.ctx=.pdf) $(STATS:.csv=.pdf)

all: $(KEEP)

plots: $(PLOTS)

seq.fa: Makefile
	echo ACACAGAGAGTCCCT > seq.fa
	echo ACACAGAGAGTCACTCCCC >> seq.fa
	echo ACACAGAGAGTCACTCCCC >> seq.fa
	echo ACACAGAGACTCACTCCCC >> seq.fa
	echo ACACAGAGACTCACTCCCC >> seq.fa

seq.k9.raw.ctx: seq.fa
	$(CTX) build -m 10M -k 9 --sample SeqJr --seq $< $@

covgs.before.csv: seq.k9.clean.ctx
covgs.after.csv: seq.k9.clean.ctx
lens.before.csv: seq.k9.clean.ctx
lens.after.csv: seq.k9.clean.ctx
seq.k9.clean.ctx: seq.k9.raw.ctx
	$(CTX) clean --covg-before covgs.before.csv --covg-after covgs.after.csv \
	             --len-before lens.before.csv --len-after lens.after.csv \
	             --supernodes=2 --tips=62 --out $@ $<

seq.k9.%.dot: seq.k9.%.ctx
	$(MKDOT) $< > $@

seq.k9.%.pdf: seq.k9.%.dot
	dot -Tpdf $< > $@

covgs.%.pdf: covgs.%.csv
	R --vanilla --file=$(PLOTCOVG) --args $< $@

lens.%.pdf: lens.%.csv
	R --vanilla --file=$(PLOTLENS) --args $< $@

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all clean plots
