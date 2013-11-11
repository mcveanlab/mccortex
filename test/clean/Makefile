
KEEP=seq.fa seq.k9.raw.ctx seq.k9.clean.ctx seq.k9.clean.pdf seq.k9.raw.pdf

all: $(KEEP)

seq.fa: Makefile
	echo ACACAGAGAGTCCCT > seq.fa
	echo ACACAGAGAGTCACTCCCC >> seq.fa
	echo ACACAGAGAGTCACTCCCC >> seq.fa
	echo ACACAGAGACTCACTCCCC >> seq.fa
	# echo ACACAGAGACTCACTCCCC >> seq.fa

seq.k9.raw.ctx: seq.fa
	ctx31 build -m 10M -k 9 --sample SeqJr --seq $< $@

seq.k9.clean.ctx: seq.k9.raw.ctx
	ctx31 clean --supernodes --threshold 2 --tips 62 $@ $<

seq.k9.%.dot: seq.k9.%.ctx
	~/cortex/versions/compact_hash/scripts/cortex_to_graphviz.pl $< > $@

seq.k9.%.pdf: seq.k9.%.dot
	dot -Tpdf $< > $@

clean:
	rm -rf $(KEEP)

.PHONY: all clean
