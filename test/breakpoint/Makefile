

#
# s0: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
#
# s1: GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG
#
# mx: CCCGTAGGTAAGGGCGTTAGCGGGTTGGAGTTGGCCAAAGAAGTTCAACG
#     11111111111111111111222222222222222222222222222222
#

K=11
CTXDIR=../../
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl

SEQS=seq0.fa seq1.fa mix.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
TGTS=breakpoints.txt $(GRAPHS)
# join.k$(K).ctx

all: $(TGTS)

seq0.fa:
	echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA > seq0.fa

seq1.fa:
	echo GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG > seq1.fa

mix.fa:
	echo CCCGTAGGTAAGGGCGTTAGCGGGTTGGAGTTGGCCAAAGAAGTTCAACG > mix.fa

%.k$(K).ctx: %.fa
	$(CTX) build -m 10M -k $(K) --sample $< --seq $< $@

breakpoints.txt.gz: $(GRAPHS) $(SEQS)
	$(CTX) breakpoints -t 1 -m 10M --seq seq0.fa --seq seq1.fa --out $@ mix.k$(K).ctx

breakpoints.txt: breakpoints.txt.gz
	gzip -d $<
	cat breakpoints.txt

join.k$(K).ctx: $(GRAPHS)
	$(CTX) join --flatten $@ $(GRAPHS)

seq.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: seq.k$(K).pdf

clean:
	rm -rf $(TGTS) $(SEQS) breakpoints.txt.gz seq.k$(K).pdf join.k$(K).ctx

.PHONY: all clean plots
