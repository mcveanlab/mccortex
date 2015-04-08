SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31

all: subgraph0.ctx subgraph1.ctx subgraph10.ctx

seed.fa:
	echo ACAATGCAGCATT > seed.fa

seq.fa:
	echo AGGGGCAGAAAATGCAGCAT > seq.fa

graph.ctx: seq.fa
	$(CTX) build -k 9 --sample MsGraph --seq seq.fa graph.ctx

subgraph%.ctx: graph.ctx seed.fa
	$(CTX) subgraph --seed seed.fa --dist $* -o subgraph$*.ctx graph.ctx

clean:
	rm -rf subgraph*.ctx graph.ctx seed.fa seq.fa

.PHONY: all clean
