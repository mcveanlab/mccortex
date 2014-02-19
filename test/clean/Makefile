
CTXDIR=../../
CTX=$(CTXDIR)/bin/ctx31
MKDOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

KEEP=seq.fa seq.k9.raw.ctx seq.k9.clean.ctx len.before.csv len.after.csv
PLOTS=seq.k9.clean.pdf seq.k9.raw.pdf

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

len.before.csv: seq.k9.clean.ctx
len.after.csv: seq.k9.clean.ctx
seq.k9.clean.ctx: seq.k9.raw.ctx
	$(CTX) clean --covgs covg.csv --supernodes --threshold 2 --tips 62 -o $@ --len-before len.before.csv --len-after len.after.csv $<

seq.k9.%.dot: seq.k9.%.ctx
	$(MKDOT) $< > $@

seq.k9.%.pdf: seq.k9.%.dot
	dot -Tpdf $< > $@

clean:
	rm -rf $(KEEP) $(PLOTS)

.PHONY: all clean plots
