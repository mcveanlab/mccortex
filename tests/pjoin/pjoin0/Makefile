SHELL:=/bin/bash -euo pipefail

CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex31

PATHS=paths.0.ctp.gz paths.1.ctp.gz
SEQ=genome.0.fa genome.1.fa
GRAPHS=$(SEQ:.fa=.ctx)
MERGED=genomes.ctx genomes.ctp.gz

TGTS=$(SEQ) $(GRAPHS) $(PATHS) $(MERGED)

# non-default target: genome.k9.pdf

all: $(TGTS)

clean:
	rm -rf $(TGTS)

genome.0.fa:
	echo TGGTGTCGCCTACA > $@
	echo TtGTGTCGCCTAgA >> $@

genome.1.fa:
	echo TtGTGTCGCCTACA > $@
	echo TGGTGTCGCCTAgA >> $@

genome.%.ctx: genome.%.fa
	$(MCCORTEX) build -q -m 1M -k 9 --sample Gnome$* --seq genome.$*.fa $@

paths.%.ctp.gz: genome.%.ctx genome.%.fa
	$(MCCORTEX) thread -q -m 1M --seq genome.$*.fa -o $@ genome.$*.ctx
	gunzip -c $@

genomes.ctx: $(PATHS)
	$(MCCORTEX) join -q -o $@ $(GRAPHS)

genomes.ctp.gz: $(PATHS)
	$(MCCORTEX) pjoin -q -o $@ $(PATHS)
	gunzip -c $@

.PHONY: all plots clean
