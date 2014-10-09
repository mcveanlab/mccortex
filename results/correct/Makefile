SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx$(K)
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
READSIM=$(CTXDIR)/libs/readsim/readsim
SIMCHECK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_corrected.pl
K=31

REF=$(CTXDIR)/results/data/chr22/chr22_17M_18M.fa
PHIX=$(CTXDIR)/results/data/PhiX/PhiX.1.fq.gz
DEPTH=50
READLEN=250

READS=reads/chrom0.$(DEPTH)X.1.fa.gz reads/chrom1.$(DEPTH)X.1.fa.gz
CORREADS=correct/chrom0.$(DEPTH)X.1.fa.gz correct/chrom1.$(DEPTH)X.1.fa.gz

TGTS=diploid/genome0.fa diploid/genome1.fa \
     diploid/chrom0.fa diploid/chrom1.fa \
     $(READS) $(CORREADS) \
     k$(K)/graph.k$(K).ctx k$(K)/cleaned.k$(K).ctx

DIRS=k$(K) correct reads diploid logs

all: $(TGTS) test

# Generate a diploid genome
diploid/genome1.fa: diploid/genome0.fa
diploid/genome0.fa: $(REF) | $(DIRS)
	$(SIMMUT) --snps 1000 --indels 10 --invs 0 diploid 2 $<

# Remove deletion marks (-) and convert to uppercase
diploid/chrom%.fa: diploid/genome%.fa
	cat $< | tr -d '-' | $(DNACAT) -u -F - > $@

# Simulate PE reads of each chrom each 50X
reads/chrom%.$(DEPTH)X.1.fa.gz: diploid/chrom%.fa | $(DIRS)
	$(READSIM) -l $(READLEN) -r $< -d $(DEPTH) -p $(PHIX) reads/chrom$*.$(DEPTH)X

# Take only the left hand read for 50X diploid coverage
# Cortex build k=$(K)
k$(K)/graph.k$(K).ctx: $(READS) | $(DIRS)
	$(CTX) build -m 100M -k $(K) --sample chr22_17M_18M --seq reads/chrom0.$(DEPTH)X.1.fa.gz --seq reads/chrom1.$(DEPTH)X.1.fa.gz $@

# Cortex clean
k$(K)/cleaned.k$(K).ctx: k$(K)/graph.k$(K).ctx
	$(CTX) clean -m 100M --covg-before logs/cleaning-covg-before.csv --covg-after logs/cleaning-covg-after.csv --len-before logs/cleaning-len-before.csv --len-after logs/cleaning-len-after.csv --out $@ $<

# Cortex correct
correct/chrom1.$(DEPTH)X.1.fa.gz: correct/chrom0.$(DEPTH)X.1.fa.gz
correct/chrom0.$(DEPTH)X.1.fa.gz: k$(K)/cleaned.k$(K).ctx $(READS) | $(DIRS)
	$(CTX) correct -m 100M --format FASTA --seq reads/chrom0.$(DEPTH)X.1.fa.gz:correct/chrom0.$(DEPTH)X.1 --seq reads/chrom1.$(DEPTH)X.1.fa.gz:correct/chrom1.$(DEPTH)X.1 k$(K)/cleaned.k$(K).ctx

$(DIRS):
	mkdir -p $@

test: $(READS) $(CORREADS)
	zcat $(READS)    | $(SIMCHECK) - diploid/chrom{0,1}.fa
	zcat $(CORREADS) | $(SIMCHECK) - diploid/chrom{0,1}.fa

clean:
	rm -rf $(TGTS) $(DIRS)

.PHONY: all clean
