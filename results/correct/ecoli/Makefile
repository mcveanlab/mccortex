SHELL:=/bin/bash -euo pipefail

CTXDIR=../../..
CTX=$(CTXDIR)/bin/ctx$(K)
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
READSIM=$(CTXDIR)/libs/readsim/readsim
SIMCHECK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_corrected.pl
K=21

REF=$(CTXDIR)/results/data/ecoli/ecoli.fa

DEPTH=100
READLEN=250

READS=reads/ecoli.$(DEPTH)X.fa.gz
CORREADS=correct/ecoli.$(DEPTH)X.fa.gz
CORREADS_PATHS=correct/ecoli.$(DEPTH)X.paths.fa.gz

TGTS=$(READS) $(CORREADS) k$(K)/graph.k$(K).ctx k$(K)/cleaned.k$(K).ctx

DIRS=k$(K) correct reads logs

all: $(TGTS) test

# Simulate SE reads
reads/ecoli.$(DEPTH)X.fa.gz: $(REF) | $(DIRS)
	$(READSIM) -l $(READLEN) -r $< -d $(DEPTH) -e 0.01 -s reads/ecoli.$(DEPTH)X

# Cortex build k=$(K)
k$(K)/graph.k$(K).ctx: $(READS) | $(DIRS)
	$(CTX) build -m 2G -k $(K) --sample ecoliK12 --seq reads/ecoli.$(DEPTH)X.fa.gz $@

# Cortex clean
k$(K)/cleaned.k$(K).ctx: k$(K)/graph.k$(K).ctx
	$(CTX) clean -m 2G --covg-before logs/cleaning-covg-before.csv --covg-after logs/cleaning-covg-after.csv --len-before logs/cleaning-len-before.csv --len-after logs/cleaning-len-after.csv --out $@ $<

# Cortex thread reads
k$(K)/cleaned.k$(K).ctp.gz: k$(K)/cleaned.k$(K).ctx
	$(CTX) thread -m 2G --seq reads/ecoli.$(DEPTH)X.fa.gz -o $@ $<

# Cortex correct
correct/ecoli.$(DEPTH)X.fa.gz: k$(K)/cleaned.k$(K).ctx $(READS) | $(DIRS)
	$(CTX) correct -m 2G --format FASTA --seq reads/ecoli.$(DEPTH)X.fa.gz:correct/ecoli.$(DEPTH)X k$(K)/cleaned.k$(K).ctx

correct/ecoli.$(DEPTH)X.paths.fa.gz: k$(K)/cleaned.k$(K).ctx k$(K)/cleaned.k$(K).ctp.gz $(READS) | $(DIRS)
	$(CTX) correct -m 2G --format FASTA --seq reads/ecoli.$(DEPTH)X.fa.gz:correct/ecoli.$(DEPTH)X -p k$(K)/cleaned.k$(K).ctp.gz k$(K)/cleaned.k$(K).ctx

$(DIRS):
	mkdir -p $@

test: $(READS) $(CORREADS)
	gzip -dc $(READS)          | $(SIMCHECK) - $(REF)
	gzip -dc $(CORREADS)       | $(SIMCHECK) - $(REF)
	gzip -dc $(CORREADS_PATHS) | $(SIMCHECK) - $(REF)

clean:
	rm -rf $(TGTS) $(DIRS)

.PHONY: all clean
