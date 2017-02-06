SHELL:=/bin/bash -euo pipefail

K=31

CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
SIMMUT=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_mutations.pl
READSIM=$(CTXDIR)/libs/readsim/readsim
SIMCHECK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_corrected.pl

REF=$(CTXDIR)/results/data/chr22/chr22_28M_29M.fa
PHIX=$(CTXDIR)/results/data/PhiX/PhiX.1.fq.gz

# simluating a diploid, actual depth 2*$(HAPDEPTH)
HAPDEPTH=25
READLEN=250

READS=reads/chrom0.$(HAPDEPTH)X.fa.gz reads/chrom1.$(HAPDEPTH)X.fa.gz
CORREADS=correct/chrom0.$(HAPDEPTH)X.fa.gz correct/chrom1.$(HAPDEPTH)X.fa.gz
CORREADS_LINKS=correct/chrom0.$(HAPDEPTH)X.links.fa.gz correct/chrom1.$(HAPDEPTH)X.links.fa.gz

TGTS=diploid/genome0.fa diploid/genome1.fa \
     diploid/chrom0.fa diploid/chrom1.fa \
     $(READS) $(CORREADS) $(CORREADS_LINKS) \
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
reads/chrom%.$(HAPDEPTH)X.fa.gz: diploid/chrom%.fa | $(DIRS)
	$(READSIM) -l $(READLEN) -r $< -d $(HAPDEPTH) -p $(PHIX) -s reads/chrom$*.$(HAPDEPTH)X

# McCortex build k=$(K)
k$(K)/graph.k$(K).ctx: $(READS) | $(DIRS)
	$(MCCORTEX) build -m 100M -k $(K) --sample chr22_17M_18M --seq reads/chrom0.$(HAPDEPTH)X.fa.gz --seq reads/chrom1.$(HAPDEPTH)X.fa.gz $@

# McCortex clean
k$(K)/cleaned.k$(K).ctx: k$(K)/graph.k$(K).ctx
	$(MCCORTEX) clean -m 100M --covg-before logs/cleaning-covg-before.csv --covg-after logs/cleaning-covg-after.csv --len-before logs/cleaning-len-before.csv --len-after logs/cleaning-len-after.csv --out $@ $<

# McCortex thread reads
k$(K)/cleaned.k$(K).ctp.gz: k$(K)/cleaned.k$(K).ctx $(READS) | $(DIRS)
	$(MCCORTEX) thread -m 100M --seq reads/chrom0.$(HAPDEPTH)X.fa.gz --seq reads/chrom1.$(HAPDEPTH)X.fa.gz -o $@ $<

# McCortex correct
correct/chrom1.$(HAPDEPTH)X.fa.gz: correct/chrom0.$(HAPDEPTH)X.fa.gz
correct/chrom0.$(HAPDEPTH)X.fa.gz: k$(K)/cleaned.k$(K).ctx $(READS) | $(DIRS)
	$(MCCORTEX) correct -m 100M --format FASTA --seq reads/chrom0.$(HAPDEPTH)X.fa.gz:correct/chrom0.$(HAPDEPTH)X --seq reads/chrom1.$(HAPDEPTH)X.fa.gz:correct/chrom1.$(HAPDEPTH)X k$(K)/cleaned.k$(K).ctx

correct/chrom1.$(HAPDEPTH)X.links.fa.gz: correct/chrom0.$(HAPDEPTH)X.fa.gz
correct/chrom0.$(HAPDEPTH)X.links.fa.gz: k$(K)/cleaned.k$(K).ctx k$(K)/cleaned.k$(K).ctp.gz $(READS) | $(DIRS)
	$(MCCORTEX) correct -m 100M --format FASTA --seq reads/chrom0.$(HAPDEPTH)X.fa.gz:correct/chrom0.$(HAPDEPTH)X.links --seq reads/chrom1.$(HAPDEPTH)X.fa.gz:correct/chrom1.$(HAPDEPTH)X.links -p k$(K)/cleaned.k$(K).ctp.gz k$(K)/cleaned.k$(K).ctx

$(DIRS):
	mkdir -p $@

test: $(READS) $(CORREADS) $(CORREADS_LINKS)
	gzip -dc $(READS)          | $(SIMCHECK) - diploid/chrom{0,1}.fa
	gzip -dc $(CORREADS)       | $(SIMCHECK) - diploid/chrom{0,1}.fa
	gzip -dc $(CORREADS_LINKS) | $(SIMCHECK) - diploid/chrom{0,1}.fa

clean:
	rm -rf $(TGTS) $(DIRS)

.PHONY: all clean
