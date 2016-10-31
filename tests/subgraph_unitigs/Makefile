SHELL:=/bin/bash -euo pipefail

K=11
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])

all: subgraph0.k$(K).ctx subgraph1.k$(K).ctx

seed.fa:
	echo TGCCTAGAAGG > seed.fa

seq.fa:
	echo ATGGTGCCTAGAAGGTA > seq.fa
	echo cTGGTGCCTAGAAGGTg >> seq.fa

graph.k$(K).ctx: seq.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample MsGraph --seq seq.fa graph.k$(K).ctx

subgraph%.k$(K).ctx: graph.k$(K).ctx seed.fa
	$(MCCORTEX) subgraph -q --seed seed.fa --unitigs --dist $* -o subgraph$*.k$(K).ctx graph.k$(K).ctx

clean:
	rm -rf subgraph*.k$(K).ctx graph.k$(K).ctx seed.fa seq.fa

.PHONY: all clean
