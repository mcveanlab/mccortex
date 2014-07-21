
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31

all: seq.k15.ctx seq.k15.dot

plots: seq.k15.pdf

seq.fa:
	echo TACATTCCCCCATAGTCGTAGGCGTTAAATACA     > seq.fa
	echo TACATTCCCCCATAGTAGTAGGCGTTAAATACA    >> seq.fa
	echo GGCGTTCGCTTATCCGGATAAGCGAACGCC       >> seq.fa
	echo ATATATATATATATATATATATATATATATATATAT >> seq.fa

seq.k15.ctx: seq.fa
	$(CTX) build -k 15 --sample MasterSeq --seq seq.fa seq.k15.ctx

seq.k15.dot: seq.k15.ctx
	$(CTX) supernodes --dot $< > $@

seq.k15.pdf: seq.k15.dot
	cat $< | dot -Tpdf > $@

clean:
	rm -rf seq.fa seq.k15.ctx seq.k15.dot seq.k15.pdf

.PHONY: all clean plots
