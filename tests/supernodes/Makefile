
CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
SEQRND=$(CTXPATH)/libs/seq_file/bin/seqrnd

K=7
PLOTS=genome.k$(K).ctx.pdf genome.k$(K).perl.pdf
KEEP=genome.fa genome.k$(K).ctx $(PLOTS:.pdf=.dot)

all: $(KEEP)

clean:
	rm -rf $(KEEP) $(PLOTS) $(PLOTS:.pdf=.dot)

# Sample random genome
genome.fa:
	$(SEQRND) 200 > genome.fa

genome.k$(K).ctx: genome.fa
	$(CTX) build -m 1M -k $(K) --sample MssrGenome --seq $< $@

genome.k$(K).ctx.dot: genome.k$(K).ctx
	$(CTX) supernodes -m 1M --graphviz --points $< > $@

genome.k$(K).perl.dot: genome.k$(K).ctx
	../../scripts/cortex_to_graphviz.pl --simplify $< > $@

%.pdf: %.dot
	dot -Tpdf $< > $@

plots: $(PLOTS)

.PHONY: all clean plots
