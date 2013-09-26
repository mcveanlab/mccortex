

all: subgraph0.ctx subgraph1.ctx subgraph10.ctx

seed.fa:
	echo ACAATGCAGCATT > seed.fa
seq.fa:
	echo AGGGGCAGAAAATGCAGCAT > seq.fa

graph.ctx: seq.fa
	ctx31 build -k 9 --sample MsGraph --seq seq.fa graph.ctx

subgraph0.ctx: graph.ctx
	ctx31 subgraph --seed seed.fa subgraph.ctx 0 graph.ctx
subgraph1.ctx: graph.ctx
	ctx31 subgraph --seed seed.fa subgraph.ctx 1 graph.ctx
subgraph10.ctx: graph.ctx
	ctx31 subgraph --seed seed.fa subgraph.ctx 10 graph.ctx

clean:
	rm -rf subgraph*.ctx graph.ctx

.PHONY: all clean
