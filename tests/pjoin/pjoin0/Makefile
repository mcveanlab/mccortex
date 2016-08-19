SHELL:=/bin/bash -euo pipefail

K=9
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])

LINKS=paths.0.k$(K).ctp.gz paths.1.k$(K).ctp.gz
SEQ=genome.0.fa genome.1.fa
GRAPHS=$(SEQ:.fa=.k$(K).ctx)
MERGED=genomes.k$(K).ctx genomes.k$(K).ctp.gz

TGTS=$(SEQ) $(GRAPHS) $(LINKS) $(MERGED)

# non-default target: genome.k$(K).pdf

all: $(TGTS)

clean:
	rm -rf $(TGTS)

genome.0.fa:
	echo TGGTGTCGCCTACA > $@
	echo TtGTGTCGCCTAgA >> $@

genome.1.fa:
	echo TtGTGTCGCCTACA > $@
	echo TGGTGTCGCCTAgA >> $@

genome.%.k$(K).ctx: genome.%.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample Gnome$* --seq genome.$*.fa $@

paths.%.k$(K).ctp.gz: genome.%.k$(K).ctx genome.%.fa
	$(MCCORTEX) thread -q -m 1M --seq genome.$*.fa -o $@ genome.$*.k$(K).ctx
	gunzip -c $@

genomes.k$(K).ctx: $(LINKS)
	$(MCCORTEX) join -q -o $@ $(GRAPHS)

genomes.k$(K).ctp.gz: $(LINKS)
	$(MCCORTEX) pjoin -q -o $@ $(LINKS)
	gunzip -c $@

.PHONY: all plots clean
