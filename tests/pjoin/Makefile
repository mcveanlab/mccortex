SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31

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
	$(CTX) build -m 1M -k 9 --sample Gnome$* --seq genome.$*.fa $@

paths.%.ctp.gz: genome.%.ctx genome.%.fa
	$(CTX) thread -m 1M --seq genome.$*.fa -o $@ genome.$*.ctx
	zcat $@

genomes.ctx: $(PATHS)
	$(CTX) join -o $@ $(GRAPHS)

genomes.ctp.gz: $(PATHS)
	$(CTX) pjoin -o $@ $(PATHS)
	zcat $@

.PHONY: all plots clean
