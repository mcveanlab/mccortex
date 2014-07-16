
TGTS=ref.fa ref.ctx ref.ctp.gz
PLOTS=ref.pdf ref.dot

all: $(TGTS)

plots: ref.dot ref.pdf

CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31
CTX2GV=$(CTXPATH)/scripts/cortex_to_graphviz.pl
KMER=31

ref.fa:
	echo CGATTGAATTCCACCGATAATGCAGATGTGAGCCTCAGCATCTACTGCTTCCTCGTCGTCGGGGACTTTTGTTGACCTACCACATGATACATGCGGCCAT >  $@
	echo ACAAGCTAAAGAAGCTAGCCAGTGCAGGCTCCCTTCAGCATCTACTGCTTCCTCGTCGTCGGGGACTAGAAACGTGACCATCGGCCACCGAAAGATAAGG >> $@
	echo ACAAGCTAAAGAAGCTAGCCAGTGCAGGCTCCCTTCAGCATCTACTGCTGGCTCGTCGTCGGGGACTAGAAACGTGACCATCGGCCACCGAAAGATAAGG >> $@

ref.ctx: ref.fa
	$(CTX) build -m 1M -k $(KMER) --sample MsSample --seq ref.fa ref.ctx

ref.ctp.gz: ref.ctx
	$(CTX) thread -m 1M --seq ref.fa -o $@ $<

ref.dot: ref.ctx
	$(CTX2GV) --simplify $< > $@

ref.pdf: ref.dot
	dot -Tpdf $< > $@

clean:
	rm -rf $(TGTS) $(PLOTS)

.PHONY: all plots clean
