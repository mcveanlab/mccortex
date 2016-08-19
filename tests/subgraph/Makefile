SHELL:=/bin/bash -euo pipefail

K=9
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])

GRAPHS=graph.one.k$(K).ctx graph.many.k$(K).ctx
SUBGRAPHS=subgraph.0.one.k$(K).ctx subgraph.0.many.k$(K).ctx \
          subgraph.1.one.k$(K).ctx subgraph.1.many.k$(K).ctx \
          subgraph.10.one.k$(K).ctx subgraph.10.many.k$(K).ctx

all: check

seed.fa:
	echo ACAATGCAGCATT > seed.fa

seq.fa:
	echo AGGGGCAGAAAATGCAGCAT > seq.fa

graph.one.k$(K).ctx: seq.fa
	$(MCCORTEX) build -q -k $(K) --sample MsGraph --seq $< $@

graph.many.k$(K).ctx: graph.one.k$(K).ctx
	$(MCCORTEX) join -q -o $@ 0:$< 2:$<

subgraph.%.one.k$(K).ctx: graph.one.k$(K).ctx seed.fa
	$(MCCORTEX) subgraph -q --seed seed.fa --dist $* -o subgraph.$*.one.k$(K).ctx $<

subgraph.%.many.k$(K).ctx: graph.many.k$(K).ctx seed.fa
	$(MCCORTEX) subgraph -q --seed seed.fa --dist $* -o subgraph.$*.many.k$(K).ctx $<

check: $(GRAPHS) $(SUBGRAPHS)
	@[ `$(MCCORTEX) view -q -k subgraph.0.one.k$(K).ctx   | awk 'END{print NR}'` -eq  2 ]
	@[ `$(MCCORTEX) view -q -k subgraph.0.many.k$(K).ctx  | awk 'END{print NR}'` -eq  2 ]
	@[ `$(MCCORTEX) view -q -k subgraph.1.one.k$(K).ctx   | awk 'END{print NR}'` -eq  3 ]
	@[ `$(MCCORTEX) view -q -k subgraph.1.many.k$(K).ctx  | awk 'END{print NR}'` -eq  3 ]
	@[ `$(MCCORTEX) view -q -k subgraph.10.one.k$(K).ctx  | awk 'END{print NR}'` -eq 12 ]
	@[ `$(MCCORTEX) view -q -k subgraph.10.many.k$(K).ctx | awk 'END{print NR}'` -eq 12 ]
	@echo "Looks good."

clean:
	rm -rf subgraph*.k$(K).ctx graph*.k$(K).ctx seed.fa seq.fa

.PHONY: all clean
