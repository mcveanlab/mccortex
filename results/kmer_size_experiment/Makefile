SHELL=/bin/bash -euo pipefail
#
# Get N50 for assemblies with perfect coverage, stochastic coverage and
# stochastic coverage + sequencing error for k=21,31,41,51,63,75,99
#
# Fetch data before running:
#   cd mccortex/results/data && ./download.sh
#
# Generate plots with:
#   cd results && ./make-csvs-and-plots.sh
#
# Isaac Turner 2016-10-28

CTXDIR=../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
DNADEMUX=$(CTXDIR)/libs/seq_file/bin/dnademux
BFC=$(CTXDIR)/libs/bfc/bfc
GENREADS=$(CTXDIR)/scripts/python/generate-reads.py
COUNT_BAD_EDGES=python $(CTXDIR)/scripts/python/count-bad-edges.py
TIME=/usr/bin/time -l

# REF=$(CTXDIR)/results/data/chr22/chr22_17M_18M.fa
REF=$(CTXDIR)/results/data/chr22/chr22_28M_29M.fa
PERFECT=data/perfect_cov
STOCH=data/stoch_cov
STOCHERR=data/stocherr_cov
CORR=data/stocherr_corr

READS_PERFECT=$(PERFECT).1.fq.gz $(PERFECT).2.fq.gz
READS_STOCH=$(STOCH).1.fq.gz $(STOCH).2.fq.gz
READS_STOCHERR=$(STOCHERR).1.fq.gz $(STOCHERR).2.fq.gz
READS_CORR=$(CORR).1.fq.gz $(CORR).2.fq.gz

MKFILE=runk.mk
KMERS=21 31 41 51 61 71 81 91 99
FRAGLEN=400
READLEN=100
DEPTH=100
ERRRATE=0.005
SEED=2380999655

# {perfect_cov,stoch_cov,stocherr_cov}/k{21,31,41,51,61,71,81,91,99}/stats.links.txt
NAMES=perfect_cov stoch_cov stocherr_cov stocherr_corr
PLAIN_STATS=$(shell for d in $(NAMES); do for k in $(KMERS); do echo $$d/k$$k/stats.plain.txt; done; done)
LINKS_STATS=$(PLAIN_STATS:plain.txt=links.txt)
DIRS=data $(NAMES)

# Keep all files
.SECONDARY:

JOBS_PERFECT=$(foreach K,$(KMERS),perfect_k$(K))
JOBS_STOCH=$(foreach K,$(KMERS),stoch_k$(K))
JOBS_STOCHERR=$(foreach K,$(KMERS),stocherr_k$(K))
JOBS_ERRCORR=$(foreach K,$(KMERS),stocherr_corr_k$(K))

all: $(JOBS_PERFECT) $(JOBS_STOCH) $(JOBS_STOCHERR) $(JOBS_ERRCORR) results/latest/bad.edges.csv

$(REF):
	cd ../data && ./download.sh

reads: $(READS_PERFECT) $(READS_STOCH) $(READS_STOCHERR) $(READS_CORR)

$(PERFECT).2.fq.gz: $(PERFECT).1.fq.gz
$(PERFECT).1.fq.gz: $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | \
	$(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) | \
	$(DNADEMUX) -z -Q $(PERFECT).1.fq.gz $(PERFECT).2.fq.gz

$(STOCH).2.fq.gz: $(STOCH).1.fq.gz
$(STOCH).1.fq.gz: $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | \
	$(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) -d $(DEPTH) | \
	$(DNADEMUX) -z -Q $(STOCH).1.fq.gz $(STOCH).2.fq.gz

$(STOCHERR).2.fq.gz: $(STOCHERR).1.fq.gz
$(STOCHERR).1.fq.gz: $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | \
	$(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) -d $(DEPTH) -e $(ERRRATE) | \
	$(DNADEMUX) -z -Q $(STOCHERR).1.fq.gz $(STOCHERR).2.fq.gz

$(CORR).2.fq.gz: $(CORR).1.fq.gz
$(CORR).1.fq.gz: $(READS_STOCHERR) | $(DIRS)
	$(TIME) $(BFC) -s 1m <($(DNACAT) -i -F $(READS_STOCHERR)) <($(DNACAT) -i -F $(READS_STOCHERR)) 2> $@.log | \
	$(DNADEMUX) -z -Q $(CORR).1.fq.gz $(CORR).2.fq.gz

# perfect_cov/k%/stats.plain.txt perfect_cov/k%/stats.links.txt: $(READS_PERFECT)
perfect_k%: $(READS_PERFECT)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=perfect_cov IN1=$(PERFECT).1.fq.gz IN2=$(PERFECT).2.fq.gz

# stoch_cov/k%/stats.plain.txt stoch_cov/k%/stats.links.txt: $(READS_STOCH)
stoch_k%: $(READS_STOCH)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stoch_cov IN1=$(STOCH).1.fq.gz IN2=$(STOCH).2.fq.gz

# stocherr_cov/k%/stats.plain.txt stocherr_cov/k%/stats.links.txt: $(READS_STOCHERR)
stocherr_k%: $(READS_STOCHERR)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stocherr_cov IN1=$(STOCHERR).1.fq.gz IN2=$(STOCHERR).2.fq.gz CLEAN=1

stocherr_corr_k%: $(READS_CORR)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stocherr_corr IN1=$(CORR).1.fq.gz IN2=$(CORR).2.fq.gz CLEAN=1

# Find the number of sequencing errors that would add a new edges between two
# existing kmers
results/latest/bad.edges.csv: $(REF)
	mkdir -p results/latest
	$(DNACAT) -P $(REF) | $(COUNT_BAD_EDGES) $(DEPTH) $(ERRRATE) > $@

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) results/latest/bad.edges.csv

.PHONY: all clean
