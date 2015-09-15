SHELL=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex31

K=31

# Mac doesn't have:
# - `sort -R` (sort with hash of key)
# - `shuf` shuffle
#  So have to use perl

all: truth.tsv dist.tsv
	diff -q truth.tsv dist.tsv
	@echo "Success."

tmp.fa:
	for i in {1..200}; do $(DNACAT) -n $(K); done | sort | uniq | \
	perl -MList::Util=shuffle -e 'print shuffle<STDIN>' > $@

beauty.fa: tmp.fa
	head -n 160 tmp.fa | \
	perl -MList::Util=shuffle -e 'print shuffle<STDIN>' > $@

beast.fa:
	tail -n 140 tmp.fa | \
	perl -MList::Util=shuffle -e 'print shuffle<STDIN>' > $@

truth.tsv: beauty.fa beast.fa
	perl -v
	perl -MList::Util=min,max -e ' \
($$a,$$b,$$c) = map {$$x = `wc -l $$_`; $$x} qw(tmp.fa beauty.fa beast.fa); \
($$N,$$H,$$T) = map {($$x) = ($$_ =~ /(\d+)/); $$x} ($$a,$$b,$$c); \
print ".\tcol0\tcol1\n"; \
print "col0\t$$H\t".max(min($$H,$$N)-max($$N-$$T,0),0)."\n"; \
print "col1\t0\t$$T\n";' > $@

dist.tsv: beauty.ctx beast.ctx
	$(MCCORTEX) dist --out $@ beauty.ctx beast.ctx

%.ctx: %.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample $* --seq $< $@

$(DIRS):
	mkdir -p $@

clean:
	rm -rf beauty.fa beast.fa tmp.fa
	rm -rf beauty.ctx beast.ctx
	rm -rf truth.tsv dist.tsv

.PHONY: all clean sams bams
