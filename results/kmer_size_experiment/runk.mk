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

GRAPH=$(DIR)/graph.k$(K).raw.ctx
LINKS=$(DIR)/graph.k$(K).raw.ctp.gz
LINKSTATS=$(DIR)/graph.k$(K).linkstats.txt

ifdef CLEAN:
  GRAPH=$(DIR)/graph.k$(K).clean.ctx
  LINKS=$(DIR)/graph.k$(K).clean.ctp.gz
endif

all: $(DIR)/stats.plain.txt $(DIR)/stats.links.txt

clean:
	rm -rf $(DIR)

$(DIR)/graph.k$(K).raw.ctx: $(INPUT) | $(DIR)
	$(MCCORTEX) build -m $(MEM) -k $(K) -s KmerExperiment -1 $< $@ >& $@.log

$(DIR)/graph.k$(K).clean.ctx: $(DIR)/graph.k$(K).raw.ctx
	$(MCCORTEX) clean -m $(MEM) -o $@ $< >& $@.log

$(DIR)/graph.k$(K).raw.ctp.gz: $(GRAPH) $(INPUT)
	$(MCCORTEX) thread -m $(MEM) -o $@ -1 $(INPUT) $(GRAPH) >& $@.log

$(DIR)/graph.k$(K).linkstats.txt: $(DIR)/graph.k$(K).raw.ctp.gz
	$(MCCORTEX) links -T 0.001 -L 1000 $(LINKS) 2> $@

$(DIR)/graph.k$(K).clean.ctp.gz: $(DIR)/graph.k$(K).raw.ctp.gz $(LINKSTATS)
	LINK_THRESH=`grep 'suggested_cutoffs=' $(LINKSTATS) | grep -oE '[0-9,]+$'`
	$(MCCORTEX) links -m $(MEM) --clean $LINK_THRESH -o $@ $< >& $@.log

$(DIR)/contigs.plain.fa: $(GRAPH)
	$(MCCORTEX) contigs -m $(MEM) -o $@ $< >& $@.log

$(DIR)/contigs.links.fa: $(GRAPH) $(LINKS)
	$(MCCORTEX) contigs -m $(MEM) -p $(LINKS) -o $@ $< >& $@.log

$(DIR)/contigs.%.rmdup.fa: $(DIR)/contigs.%.fa
	$(MCCORTEX31) rmsubstr -m 2G -o $@ $< >& $@.log

$(DIR)/stats.%.txt: $(DIR)/contigs.%.rmdup.fa
	$(DNACAT) -P $(REF) $< | $(PYSTATS) $(K) 2> $@ 1> $(DIR)/stats.$*.out

$(DIR):
	mkdir -p $(DIR)

.PHONY: all clean
