#
# Sanity check that merging matching link files gives correct counts
#

SHELL:=/bin/bash -euo pipefail

K=7
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

REFLEN=5000

TGTS=genome0.fa genome0.k$(K).ctx genome0.k$(K).ctp.gz \
     genome1.fa genome1.k$(K).ctx genome1.k$(K).ctp.gz \
     joint.k$(K).ctp.gz

all: joint.k$(K).ctp.gz

genome%.fa:
	$(DNACAT) -n $(REFLEN) -M <(echo ref) -F > $@

genome%.k$(K).ctx: genome%.fa
	$(MCCORTEX) build -q -k $(K) --sample Genome0 -1 $< $@

genome%.k$(K).ctp.gz: genome%.k$(K).ctx genome%.fa
	$(MCCORTEX) thread -q -o $@ -1 genome$*.fa genome$*.k$(K).ctx

joint.k$(K).ctp.gz: genome0.k$(K).ctp.gz genome1.k$(K).ctp.gz
	$(MCCORTEX) pjoin -q -n 1M -o $@ genome0.k$(K).ctp.gz genome0.k$(K).ctp.gz genome1.k$(K).ctp.gz genome0.k$(K).ctp.gz genome1.k$(K).ctp.gz

clean:
	rm -rf $(TGTS)

.PHONY: all clean
