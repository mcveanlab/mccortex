SHELL:=/bin/bash -euo pipefail
# Makefile
# Isaac Turner
# 2014-10-01
# Using 1Mb of chr22 and error free reads to measure effect of kmer-size

KMERS=15 21 31 41 51 63 75 99

# Expand these for all k values
DIRS=reads logs $(foreach K,$(KMERS),k$(K)) plots

all:
	./run.sh

plots:
	mkdir -p plots
	R --vanilla -f plot-results.R --args {perf,stoch,stocherr}.{links,plain}.join.csv
	cd report && $(MAKE)

clean:
	rm -rf $(DIRS) {perf,stoch,stocherr}.{plain,links}.join.csv

.PHONY: all clean plots
