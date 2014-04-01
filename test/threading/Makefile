
all: rnd.ctx rnd.ctp

plots: rnd.dot rnd.pdf

CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
CTX2GV=$(CTXPATH)/scripts/cortex_to_graphviz.pl
KMER=31

rnd.ctx: rnd.fa
	$(CTX) build -k $(KMER) -m 100MB --sample MsSample --seq rnd.fa rnd.ctx

rnd.dot: rnd.ctx
	$(CTX2GV) --simplify $< > $@

rnd.pdf: rnd.dot
	dot -Tpdf $< > $@

rnd.ctp: rnd.ctx
	$(CTX) thread --col 0 --seq rnd.fa -o $@ $<

clean:
	rm -rf *.pdf *.dot *.ctx *.ctp gap_sizes.*.csv mp_sizes.*.csv

.PHONY: all plots clean
