CTX=../../bin/ctx31

all: pop.ctx noedges.fixed.ctx noedges.nofix.ctx

edges.fa:
	echo AATATGATGGCACGTGCCAACC > edges.fa
noedges.fa:
	echo AATATGATGGCACGTGCCAAC > noedges.fa
	echo ATATGATGGCACGTGCCAACC >> noedges.fa
edges.ctx: edges.fa
	$(CTX) build -k 21 --sample MrSamples --seq edges.fa edges.ctx
noedges.ctx: noedges.fa
	$(CTX) build -k 21 --sample MrsSamples --seq noedges.fa noedges.ctx

noedges.nofix.ctx: noedges.ctx
	cp noedges.ctx noedges.nofix.ctx
	$(CTX) inferedges --pop noedges.nofix.ctx
	$(CTX) view --kmers noedges.ctx
	$(CTX) view --kmers noedges.nofix.ctx
noedges.fixed.ctx: noedges.ctx
	cp noedges.ctx noedges.fixed.ctx
	$(CTX) inferedges --all noedges.fixed.ctx
	$(CTX) view --kmers noedges.ctx
	$(CTX) view --kmers noedges.fixed.ctx

pair.ctx: edges.ctx noedges.ctx
	$(CTX) join pair.ctx edges.ctx noedges.ctx
pop.ctx: pair.ctx
	cp pair.ctx pop.ctx
	$(CTX) inferedges pop.ctx
	$(CTX) view --kmers pair.ctx
	$(CTX) view --kmers pop.ctx

clean:
	rm -rf edges.fa noedges.fa edges.ctx noedges.ctx
	rm -rf noedges.fixed.ctx noedges.nofix.ctx pair.ctx pop.ctx

.PHONY: all clean
