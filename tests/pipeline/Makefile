SHELL=/bin/bash -euo pipefail


K=11
CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
MKPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl

SEQ_FILES=ben.fa jerry.fa

all: job.mk $(SEQ_FILES)
	$(MAKE) -f job.mk CTXDIR=$(CTXDIR) CLEANING_ARGS='--supernodes=1'

%.fa:
	$(DNACAT) -n 200 > $@.tmp
	$(DNACAT) $@.tmp $@.tmp $@.tmp | $(DNACAT) -F - > $@
	rm $@.tmp

samples.txt:
	( echo "ben ben.fa . ."; echo "jerry jerry.fa" ) > $@

job.mk: samples.txt
	$(MKPIPELINE) 31 proj $< > $@

clean:
	rm -rf $(SEQ_FILES) samples.txt job.mk proj

.PHONY: all clean
