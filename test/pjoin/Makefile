
CTX=../../bin/ctx31
PATHS=paths.0.ctp paths.1.ctp
SEQ=genome.0.fa genome.1.fa
GRAPHS=$(SEQ:.fa=.ctx)
MERGED=genomes.ctx genomes.ctp

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

paths.%.ctp: genome.%.ctx genome.%.fa
	$(CTX) thread -m 1M --seq genome.$*.fa -o $@ genome.$*.ctx
	$(CTX) pview $@

genomes.ctx: $(PATHS)
	$(CTX) join -o $@ $(GRAPHS)

genomes.ctp: $(PATHS)
	$(CTX) pjoin -o $@ $(PATHS)
	$(CTX) pview $@

.PHONY: all plots clean
