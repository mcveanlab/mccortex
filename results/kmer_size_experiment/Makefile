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

KMERS=21 31 41 51 61 71 81 91 99
CTXDIR=../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
GENREADS=$(CTXDIR)/scripts/python/generate-reads.py
COUNT_BAD_EDGES=python $(CTXDIR)/scripts/python/count-bad-edges.py
REF=$(CTXDIR)/results/data/chr22/chr22_17M_18M.fa
READS_PERFECT=data/perfect_cov.fa.gz
READS_STOCH=data/stoch_cov.fa.gz
READS_STOCHERR=data/stocherr_cov.fa.gz
MKFILE=runk.mk
FRAGLEN=400
READLEN=100
DEPTH=100
ERRRATE=0.005
SEED=2380999655

# {perfect_cov,stoch_cov,stocherr_cov}/k{21,31,41,51,61,71,81,91,99}/stats.links.txt
NAMES=perfect_cov stoch_cov stocherr_cov
PLAIN_STATS=$(shell for d in $(NAMES); do for k in $(KMERS); do echo $$d/k$$k/stats.plain.txt; done; done)
LINKS_STATS=$(PLAIN_STATS:plain.txt=links.txt)
DIRS=data $(NAMES)

# Keep all files
.SECONDARY:

PERFECT=$(foreach K,$(KMERS),perfect_k$(K))
STOCH=$(foreach K,$(KMERS),stoch_k$(K))
STOCHERR=$(foreach K,$(KMERS),stocherr_k$(K))

# all: $(PLAIN_STATS) $(LINKS_STATS) bad.edges.csv
all: $(PERFECT) $(STOCH) $(STOCHERR) bad.edges.csv

$(REF):
	cd ../data && ./download.sh

$(READS_PERFECT): $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | $(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) | gzip -c > $@

$(READS_STOCH): $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | $(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) -d $(DEPTH) | gzip -c > $@

$(READS_STOCHERR): $(REF) | $(DIRS)
	$(DNACAT) -P $(REF) | $(GENREADS) -s $(SEED) -p $(FRAGLEN) -r $(READLEN) -d $(DEPTH) -e $(ERRRATE) | gzip -c > $@

# perfect_cov/k%/stats.plain.txt perfect_cov/k%/stats.links.txt: $(READS_PERFECT)
perfect_k%: $(READS_PERFECT)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=perfect_cov INPUT=$(READS_PERFECT)

# stoch_cov/k%/stats.plain.txt stoch_cov/k%/stats.links.txt: $(READS_STOCH)
stoch_k%: $(READS_STOCH)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stoch_cov INPUT=$(READS_STOCH)

# stocherr_cov/k%/stats.plain.txt stocherr_cov/k%/stats.links.txt: $(READS_STOCHERR)
stocherr_k%: $(READS_STOCHERR)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stocherr_cov INPUT=$(READS_STOCHERR) CLEAN=1

# Find the number of sequencing errors that would add a new edges between two
# existing kmers
bad.edges.csv: $(REF)
	$(DNACAT) -P $(REF) | $(COUNT_BAD_EDGES) $(DEPTH) $(ERRRATE) > $@

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) bad.edges.csv

.PHONY: all clean
