#
# We assume fragment length is 400bp for PE threading
#

REQFIELDS=INPUT NAME K REF

ifndef INPUT
  $(error "Error: you need to pass the input file INPUT=$(INPUT) ($(REQFIELDS))")
endif

ifndef NAME
  $(error "Error: you need to pass the project dir NAME= ($(REQFIELDS))")
endif

ifndef K
  $(error "Error: you need to pass kmer size K= ($(REQFIELDS))")
endif

ifndef REF
  $(error "Error: you need to pass the ref file REF= ($(REQFIELDS))")
endif

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)
MCCORTEX31=$(CTXDIR)/bin/mccortex 31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
PYSTATS=python $(CTXDIR)/scripts/python/break-contigs-vs-truth.py

DIR=$(NAME)/k$(K)
MEM=1G

RAWGRAPH=$(DIR)/graph.k$(K).raw.ctx
CLEANGRAPH=$(DIR)/graph.k$(K).clean.ctx

GRAPH=$(RAWGRAPH)
SELINKS=$(DIR)/graph.k$(K).se.raw.ctp.gz
PELINKS=$(DIR)/graph.k$(K).pe.raw.ctp.gz
LINKSTATS=$(DIR)/graph.k$(K).linkstats.txt
PERFECTGRAPH=perfect_cov/k$(K)/graph.k$(K).raw.ctx
TGTS=$(DIR)/stats.plain.txt $(DIR)/stats.links.txt $(DIR)/stats.pe.txt

ifdef CLEAN
  GRAPH=$(CLEANGRAPH)
  SELINKS=$(DIR)/graph.k$(K).se.clean.ctp.gz
  PELINKS=$(DIR)/graph.k$(K).pe.clean.ctp.gz
  TGTS := $(TGTS) $(DIR)/graph.k$(K).dist.txt
endif

# Keep all files
.SECONDARY:

all: $(TGTS)

clean:
	rm -rf $(DIR)

$(DIR)/graph.k$(K).raw.ctx: $(INPUT) | $(DIR)
	$(MCCORTEX) build -m $(MEM) -k $(K) -s KmerExperiment -1 $(INPUT) $@ >& $@.log

$(DIR)/graph.k$(K).clean.ctx: $(DIR)/graph.k$(K).raw.ctx
	$(MCCORTEX) clean -m $(MEM) --fallback 3 -o $@ $< >& $@.log

$(DIR)/graph.k$(K).se.raw.ctp.gz: $(GRAPH) $(INPUT)
	$(MCCORTEX) thread -m $(MEM) -o $@ -1 $(INPUT) $(GRAPH) >& $@.log

$(DIR)/graph.k$(K).pe.raw.ctp.gz: $(GRAPH) $(INPUT) $(SELINKS)
	$(MCCORTEX) thread -m $(MEM) -p $(SELINKS) -0 -l 350 -L 450 -o $@ -i $(INPUT) $(GRAPH) >& $@.log

$(DIR)/graph.k$(K).linkstats.txt: $(DIR)/graph.k$(K).se.raw.ctp.gz
	$(MCCORTEX) links -T $@ -L 1000 $< 2> $@.log

$(DIR)/graph.k$(K).%.clean.ctp.gz: $(DIR)/graph.k$(K).%.raw.ctp.gz $(LINKSTATS)
	( LINK_THRESH=`grep 'suggested_cutoff=' $(LINKSTATS) | grep -oE '[0-9,]+$$'`; \
	  $(MCCORTEX) links --clean $$LINK_THRESH -o $@ $< >& $@.log )

$(DIR)/graph.k$(K).dist.txt: $(RAWGRAPH) $(CLEANGRAPH) $(PERFECTGRAPH)
	$(MCCORTEX) dist -m $(MEM) -o $@ $(RAWGRAPH) $(CLEANGRAPH) $(PERFECTGRAPH) >& $@.log

$(DIR)/contigs.plain.fa: $(GRAPH)
	$(MCCORTEX) contigs -m $(MEM) -o $@ $< >& $@.log

$(DIR)/contigs.links.fa: $(GRAPH) $(SELINKS)
	$(MCCORTEX) contigs -m $(MEM) -p $(SELINKS) -o $@ $< >& $@.log

$(DIR)/contigs.pe.fa: $(GRAPH) $(PELINKS)
	$(MCCORTEX) contigs -m $(MEM) -p $(PELINKS) -o $@ $< >& $@.log

$(DIR)/contigs.%.rmdup.fa: $(DIR)/contigs.%.fa
	$(MCCORTEX31) rmsubstr -m $(MEM) -n 50M -k 21 -o $@ $< >& $@.log

$(DIR)/stats.%.txt: $(DIR)/contigs.%.rmdup.fa
	$(DNACAT) -P $(REF) $< | $(PYSTATS) 21 2> $@ 1> $(DIR)/stats.$*.out

$(DIR):
	mkdir -p $(DIR)

.PHONY: all clean
