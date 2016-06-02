SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex31

all: subgraph0.ctx subgraph1.ctx

seed.fa:
	echo TGCCTAGAAGG > seed.fa

seq.fa:
	echo ATGGTGCCTAGAAGGTA > seq.fa
	echo cTGGTGCCTAGAAGGTg >> seq.fa

graph.ctx: seq.fa
	$(MCCORTEX) build -m 1M -k 11 --sample MsGraph --seq seq.fa graph.ctx

subgraph%.ctx: graph.ctx seed.fa
	$(MCCORTEX) subgraph --seed seed.fa --unitigs --dist $* -o subgraph$*.ctx graph.ctx

clean:
	rm -rf subgraph*.ctx graph.ctx seed.fa seq.fa

.PHONY: all clean
