SHELL:=/bin/bash -euo pipefail

K=9
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
MKDOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
PLOTCOVG=$(CTXDIR)/scripts/plot-covg-hist.R
PLOTLENS=$(CTXDIR)/scripts/plot-length-hist.R

SEQ=seq.fa
GRAPHS=seq.k$(K).raw.ctx seq.k$(K).clean.ctx
STATS=lens.before.k$(K).csv lens.after.k$(K).csv covgs.before.k$(K).csv covgs.after.k$(K).csv
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

seq.k$(K).raw.ctx: seq.fa
	$(MCCORTEX) build -q -m 10M -k $(K) --sample SeqJr --seq $< $@

covgs.before.k$(K).csv: seq.k$(K).clean.ctx
covgs.after.k$(K).csv: seq.k$(K).clean.ctx
lens.before.k$(K).csv: seq.k$(K).clean.ctx
lens.after.k$(K).csv: seq.k$(K).clean.ctx
seq.k$(K).clean.ctx: seq.k$(K).raw.ctx
	$(MCCORTEX) clean -q --covg-before covgs.before.k$(K).csv \
	                     --covg-after covgs.after.k$(K).csv \
	                     --len-before lens.before.k$(K).csv \
	                     --len-after lens.after.k$(K).csv \
	                  --unitigs=2 --tips=62 --out $@ $<

seq.k$(K).%.dot: seq.k$(K).%.ctx
	$(MKDOT) $< > $@

seq.k$(K).%.pdf: seq.k$(K).%.dot
	dot -Tpdf $< > $@

covgs.%.pdf: covgs.%.csv
	R --vanilla --file=$(PLOTCOVG) --args $< $@

lens.%.pdf: lens.%.csv
	R --vanilla --file=$(PLOTLENS) --args $< $@

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all clean plots
