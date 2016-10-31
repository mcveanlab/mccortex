SHELL=/bin/bash -euo pipefail

K=31
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

# Mac doesn't have:
# - `sort -R` (sort with hash of key)
# - `shuf` shuffle
#  So have to use perl
SHUFFLE=perl -MList::Util=shuffle -e 'print shuffle<STDIN>'

all: truth.tsv dist.tsv
	diff -q truth.tsv dist.tsv
	@echo "Success."

tmp.fa:
	for i in {1..200}; do $(DNACAT) -n $(K); done | sort | uniq | $(SHUFFLE) > $@

beauty.fa: tmp.fa
	head -n 160 tmp.fa | $(SHUFFLE) > $@

beast.fa:
	tail -n 140 tmp.fa | $(SHUFFLE) > $@

truth.tsv: beauty.fa beast.fa
	perl -MList::Util=min,max -W -e 'use strict; use warnings; '\
'my @files = ("tmp.fa","beauty.fa","beast.fa"); '\
'my ($$a,$$b,$$c) = map {my $$x = `wc -l $$_`; $$x} @files; '\
'my ($$N,$$H,$$T) = map {my ($$x) = ($$_ =~ /(\d+)/); $$x} ($$a,$$b,$$c); '\
'print ".\tcol0\tcol1\n"; '\
'print "col0\t$$H\t".max(min($$H,$$N)-max($$N-$$T,0),0)."\n"; '\
'print "col1\t.\t$$T\n";' > $@

dist.tsv: beauty.ctx beast.ctx
	$(MCCORTEX) dist -q --out $@ beauty.ctx beast.ctx

%.ctx: %.fa
	$(MCCORTEX) build -q -m 1M -k $(K) --sample $* --seq $< $@

$(DIRS):
	mkdir -p $@

clean:
	rm -rf beauty.fa beast.fa tmp.fa
	rm -rf beauty.ctx beast.ctx
	rm -rf truth.tsv dist.tsv

.PHONY: all clean sams bams
