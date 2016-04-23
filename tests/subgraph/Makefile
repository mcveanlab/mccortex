SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex31

GRAPHS=graph.one.ctx graph.many.ctx
SUBGRAPHS=subgraph.0.one.ctx subgraph.0.many.ctx \
          subgraph.1.one.ctx subgraph.1.many.ctx \
          subgraph.10.one.ctx subgraph.10.many.ctx

all: check

seed.fa:
	echo ACAATGCAGCATT > seed.fa

seq.fa:
	echo AGGGGCAGAAAATGCAGCAT > seq.fa

graph.one.ctx: seq.fa
	$(MCCORTEX) build -q -k 9 --sample MsGraph --seq $< $@

graph.many.ctx: graph.one.ctx
	$(MCCORTEX) join -q -o $@ 0:$< 2:$<

subgraph.%.one.ctx: graph.one.ctx seed.fa
	$(MCCORTEX) subgraph -q --seed seed.fa --dist $* -o subgraph.$*.one.ctx $<

subgraph.%.many.ctx: graph.many.ctx seed.fa
	$(MCCORTEX) subgraph -q --seed seed.fa --dist $* -o subgraph.$*.many.ctx $<

check: $(GRAPHS) $(SUBGRAPHS)
	@[ `$(MCCORTEX) view -q -k subgraph.0.one.ctx   | awk 'END{print NR}'` -eq  2 ]
	@[ `$(MCCORTEX) view -q -k subgraph.0.many.ctx  | awk 'END{print NR}'` -eq  2 ]
	@[ `$(MCCORTEX) view -q -k subgraph.1.one.ctx   | awk 'END{print NR}'` -eq  3 ]
	@[ `$(MCCORTEX) view -q -k subgraph.1.many.ctx  | awk 'END{print NR}'` -eq  3 ]
	@[ `$(MCCORTEX) view -q -k subgraph.10.one.ctx  | awk 'END{print NR}'` -eq 12 ]
	@[ `$(MCCORTEX) view -q -k subgraph.10.many.ctx | awk 'END{print NR}'` -eq 12 ]
	@echo "Looks good."

clean:
	rm -rf subgraph*.ctx graph*.ctx seed.fa seq.fa

.PHONY: all clean
