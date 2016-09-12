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
# Isaac Turner 2016-10-12

KMERS=21 31 41 51 61 71 81 91 99
CTXDIR=../../
READSIM=$(CTXDIR)/libs/readsim/readsim
PERFECTSIM=$(CTXDIR)/libs/seq_file/scripts/perfect_covg.sh
REF=$(CTXDIR)/results/data/chr22/chr22_17M_18M.fa
READS_PERFECT=data/perfect_cov.fa.gz
READS_STOCH=data/stoch_cov.fa.gz
READS_STOCHERR=data/stocherr_cov.fa.gz
MKFILE=runk.mk
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

all: $(PLAIN_STATS) $(LINKS_STATS)

$(REF):
	cd ../data && ./download.sh

$(READS_PERFECT): $(REF) | $(DIRS)
	$(PERFECTSIM) $(READLEN) $< | gzip -c > $@

$(READS_STOCH): $(REF) | $(DIRS)
	$(READSIM) -g $(SEED) -l $(READLEN) -r $< -d $(DEPTH) -e 0 -s data/stoch_cov

$(READS_STOCHERR): $(REF) | $(DIRS)
	$(READSIM) -g $(SEED) -l $(READLEN) -r $< -d $(DEPTH) -e $(ERRRATE) -s data/stocherr_cov

perfect_cov/k%/stats.plain.txt perfect_cov/k%/stats.links.txt: $(READS_PERFECT)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=perfect_cov INPUT=$(READS_PERFECT)

stoch_cov/k%/stats.plain.txt stoch_cov/k%/stats.links.txt: $(READS_STOCH)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stoch_cov INPUT=$(READS_STOCH)

stocherr_cov/k%/stats.plain.txt stocherr_cov/k%/stats.links.txt: $(READS_STOCHERR)
	$(MAKE) -f $(MKFILE) K=$* REF=$(REF) NAME=stocherr_cov INPUT=$(READS_STOCHERR) CLEAN=1

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS)

.PHONY: all clean
