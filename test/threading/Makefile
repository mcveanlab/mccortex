
all: rnd.ctx rnd.ctp rnd.dot rnd.pdf

CTX=~/cortex/versions/compact_hash/bin/ctx31
CTX2GV=~/cortex/versions/compact_hash/scripts/cortex_to_graphviz.pl --simplify
KMER=31

rnd.ctx: rnd.fa
	$(CTX) build -k $(KMER) -m 100MB --sample MsSample --seq rnd.fa rnd.ctx

rnd.dot: rnd.ctx
	$(CTX2GV) rnd.ctx > rnd.dot

rnd.pdf: rnd.dot
	dot -Tpdf rnd.dot > rnd.pdf

rnd.ctp: rnd.ctx
	$(CTX) thread --col 0 --seq rnd.fa rnd.ctp rnd.ctx

clean:
	rm -rf *.pdf *.dot *.ctx *.ctp gap_sizes.*.csv mp_sizes.*.csv

.PHONY: all clean
