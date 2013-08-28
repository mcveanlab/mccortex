
all: rnd.pdf rnd.bubbles.gz rnd.clean.pdf

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
	$(CTX) thread --col 0 0 --seq rnd.fa 1 rnd.ctp rnd.ctx rnd.ctx

rnd.bubbles.gz: rnd.ctx rnd.ctp
	$(CTX) call -t 1 -m 100MB -p rnd.ctp rnd.ctx rnd.bubbles.gz

rnd.clean.ctx: rnd.ctx
	$(CTX) clean --tips 40 --supernodes --threshold 2 --covgs covgs.csv rnd.clean.ctx rnd.ctx

rnd.clean.dot: rnd.clean.ctx
	$(CTX2GV) rnd.clean.ctx > rnd.clean.dot

rnd.clean.pdf: rnd.clean.dot
	dot -Tpdf rnd.clean.dot > rnd.clean.pdf

clean:
	rm -rf *.pdf *.dot *.ctx *.ctp gap_sizes.*.csv mp_sizes.*.csv covgs.csv rnd.bubbles.gz

.PHONY: all clean
