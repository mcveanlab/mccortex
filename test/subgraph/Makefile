
CTXPATH=../../
CTX=$(CTXPATH)/bin/ctx31

all: subgraph0.ctx subgraph1.ctx subgraph10.ctx

seed.fa:
	echo ACAATGCAGCATT > seed.fa

seq.fa:
	echo AGGGGCAGAAAATGCAGCAT > seq.fa

graph.ctx: seq.fa
	$(CTX) build -k 9 --sample MsGraph --seq seq.fa graph.ctx

subgraph0.ctx: graph.ctx seed.fa
	$(CTX) subgraph --seed seed.fa -o subgraph0.ctx 0 graph.ctx

subgraph1.ctx: graph.ctx seed.fa
	$(CTX) subgraph --seed seed.fa -o subgraph1.ctx 1 graph.ctx

subgraph10.ctx: graph.ctx seed.fa
	$(CTX) subgraph --seed seed.fa -o subgraph10.ctx 10 graph.ctx

clean:
	rm -rf subgraph*.ctx graph.ctx seed.fa seq.fa

.PHONY: all clean
