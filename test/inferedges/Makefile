
all: pop.ctx noedges.fixed.ctx noedges.nofix.ctx

edges.fa:
	echo AATATGATGGCACGTGCCAACC > edges.fa
noedges.fa:
	echo AATATGATGGCACGTGCCAAC > noedges.fa
	echo ATATGATGGCACGTGCCAACC >> noedges.fa
edges.ctx: edges.fa
	ctx31 build -k 21 --sample MrSamples --seq edges.fa edges.ctx
noedges.ctx: noedges.fa
	ctx31 build -k 21 --sample MrsSamples --seq noedges.fa noedges.ctx

noedges.nofix.ctx: noedges.ctx
	cp noedges.ctx noedges.nofix.ctx
	ctx31 inferedges --pop noedges.nofix.ctx
	ctx31 view --kmers noedges.ctx
	ctx31 view --kmers noedges.nofix.ctx
noedges.fixed.ctx: noedges.ctx
	cp noedges.ctx noedges.fixed.ctx
	ctx31 inferedges --all noedges.fixed.ctx
	ctx31 view --kmers noedges.ctx
	ctx31 view --kmers noedges.fixed.ctx

pair.ctx: edges.ctx noedges.ctx
	ctx31 join pair.ctx edges.ctx noedges.ctx
pop.ctx: pair.ctx
	cp pair.ctx pop.ctx
	ctx31 inferedges pop.ctx
	ctx31 view --kmers pair.ctx
	ctx31 view --kmers pop.ctx

clean:
	rm -rf edges.fa noedges.fa edges.ctx noedges.ctx
	rm -rf noedges.fixed.ctx noedges.nofix.ctx pair.ctx pop.ctx

.PHONY: all clean
