SHELL=/bin/bash -euo pipefail

# 1. Sample seq from ref with error
# 2. Build + clean graph
# 3. Build + clean links
# 4. Assemble contigs
# 5. check contigs match seq

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat


K=17
# Coverage must be at least one, bash {1..N} gives '1 0' (i.e. 2) if N is 0
REFCOV=20
ERRCOV=2

SEQ=ref.fa err.fa reads.fa
GRAPHS=graph.raw.k$(K).ctx graph.clean.k$(K).ctx
LINKS=graph.raw.k$(K).ctp.gz graph.clean.k$(K).ctp.gz
CONTIGS=contigs.raw.fa contigs.fa
LOGS=$(GRAPHS:.ctx=.ctx.log) $(LINKS:.ctp.gz=.ctp.gz.log) $(CONTIGS:.fa=.fa.log)
DOTS=$(GRAPHS:.ctx=.dot)
PDFS=$(DOTS:.dot=.pdf)

FILES=$(SEQ) $(GRAPHS) $(LINKS) $(CONTIGS) $(LOGS)

all: test

ref.fa:
	(echo "CACAA""TAACATTAGGATGGGTGTCC""TAGGA"; \
	 echo "CTCCC""TAACATTAGGATGGGTGTCC""AGCAT"; ) > $@

err.fa:
	(echo "CACAA""TAACATTAGGATGGGTGTCC""AGCAT"; \
	 echo "CTCCC""TAACATTAGGATGGGTGTCC""TAGGA"; ) > $@

reads.fa: ref.fa err.fa
	( for i in {1..$(REFCOV)}; do cat ref.fa; done; \
		for i in {1..$(ERRCOV)}; do cat err.fa; done; ) > $@

graph.raw.k$(K).ctx: reads.fa
	$(MCCORTEX) build -m 10M -k $(K) --sample SeqUniq --seq $< $@ >& $@.log

# Cleaning should not remove any kmers
graph.clean.k$(K).ctx: graph.raw.k$(K).ctx
	$(MCCORTEX) clean --unitigs=5 --tips=0 --out $@ $< >& $@.log

graph.raw.k$(K).ctp.gz: graph.clean.k$(K).ctx reads.fa
	$(MCCORTEX) thread --seq reads.fa --out $@ graph.clean.k$(K).ctx >& $@.log

graph.clean.k$(K).ctp.gz: graph.raw.k$(K).ctp.gz
	$(MCCORTEX) links --clean 5 --out $@ $< >& $@.log

contigs.raw.fa: graph.clean.k$(K).ctx graph.clean.k$(K).ctp.gz
	$(MCCORTEX) contigs -q --no-missing-check -o $@ -p graph.clean.k$(K).ctp.gz graph.clean.k$(K).ctx

contigs.fa: contigs.raw.fa
	$(MCCORTEX) rmsubstr -q -n 1M -k $(K) $< > $@

test: contigs.fa
	@echo Checking if regenerated file matches original...
	diff -q <($(DNACAT) -r -k -P ref.fa | sort) <($(DNACAT) -r -k -P contigs.fa | sort)
	@echo "All looks good."

%.dot: %.ctx
	$(MCCORTEX) unitigs --dot --out $< $@

%.pdf: %.dot
	dot -Tpdf $< > $@

plots: $(PDFS)

clean:
	rm -rf $(FILES) $(DOTS) $(PDFS)

.PHONY: all test clean
