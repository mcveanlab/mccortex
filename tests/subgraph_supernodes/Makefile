
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31

all: subgraph0.ctx subgraph1.ctx

seed.fa:
	echo TGCCTAGAAGG > seed.fa

seq.fa:
	echo ATGGTGCCTAGAAGGTA > seq.fa
	echo cTGGTGCCTAGAAGGTg >> seq.fa

graph.ctx: seq.fa
	$(CTX) build -m 1M -k 11 --sample MsGraph --seq seq.fa graph.ctx

subgraph%.ctx: graph.ctx seed.fa
	$(CTX) subgraph --seed seed.fa --supernodes --dist $* -o subgraph$*.ctx graph.ctx

clean:
	rm -rf subgraph*.ctx graph.ctx seed.fa seq.fa

.PHONY: all clean
