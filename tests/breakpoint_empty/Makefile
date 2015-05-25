#
# Sanity check that calling breakpoints with only a ref and reads from the ref
# gives no calls
#

SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTXPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
READSIM=$(CTXDIR)/libs/readsim/readsim
VCFCOUNT=$(CTXDIR)/libs/biogrok/vcf-count

REFLEN=1000
K=31
SEQDEPTH=30
READLEN=100
OUTDIR=proj

all: run

ref.fa:
	$(DNACAT) -n $(REFLEN) -M <(echo ref) -F > $@

reads/reads.fa.gz: ref.fa
	mkdir -p reads
	$(READSIM) -r ref.fa -s -d $(SEQDEPTH) reads/reads

task.k$(K).mk:
	echo "RefReads reads/reads.fa.gz" | $(CTXPIPELINE) -r ref.fa $(K) proj - > $@

run: task.k$(K).mk reads/reads.fa.gz ref.fa
	$(MAKE) -f $< CTXDIR=$(CTXDIR) breakpoints-vcf
	@# Check no VCF entries
	(( `$(VCFCOUNT) proj/vcfs/breakpoints.k$(K).vcf.gz` == 0 )) || false
	@# Check no breakpoint call entries
	(( `grep -c '>brkpnt' proj/k$(K)/breakpoints/breakpoints.txt.gz` == 0 )) || false
	@echo 'Success - no breakpoint calls or VCF entries!'

clean:
	rm -rf ref.fa* reads proj task.k$(K).mk

.PHONY: all run clean
