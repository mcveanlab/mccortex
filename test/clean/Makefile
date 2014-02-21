
CTXDIR=../../
CTX=$(CTXDIR)/bin/ctx31
MKDOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
PLOTCOVG=$(CTXDIR)/scripts/plot-covg-hist.R
PLOTLENS=$(CTXDIR)/scripts/plot-length-hist.R

KEEP=seq.fa seq.k9.raw.ctx seq.k9.clean.ctx len.before.csv len.after.csv covgs.csv
PLOTS=seq.k9.clean.pdf seq.k9.raw.pdf len.before.pdf len.after.pdf covgs.pdf

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

covgs.csv: seq.k9.clean.ctx
len.before.csv: seq.k9.clean.ctx
len.after.csv: seq.k9.clean.ctx
seq.k9.clean.ctx: seq.k9.raw.ctx
	$(CTX) clean --covgs covgs.csv --supernodes --threshold 2 --tips 62 -o $@ --len-before len.before.csv --len-after len.after.csv $<

seq.k9.%.dot: seq.k9.%.ctx
	$(MKDOT) $< > $@

seq.k9.%.pdf: seq.k9.%.dot
	dot -Tpdf $< > $@

covgs.pdf: covgs.csv
	R --vanilla --file=$(PLOTCOVG) --args $< $@

len.%.pdf: len.%.csv
	R --vanilla --file=$(PLOTLENS) --args $< $@

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all clean plots
