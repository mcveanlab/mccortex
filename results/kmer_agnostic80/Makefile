SHELL:=/bin/bash -euo pipefail
# Makefile
# Isaac Turner
# 2014-10-01
# Using 1Mb of chr22 and error free reads to measure effect of kmer-size

KMERS=15 21 31 41 51 63 75 99

# Expand these for all k values
DIRS=reads $(foreach K,$(KMERS),k$(K)) plots

all:
	./run.sh
plots:
	mkdir -p plots
	R --vanilla -f plot-results.R --args gP.l0.plain.stats.csv \
	                                     gS.l0.plain.stats.csv \
	                                     gE.l0.plain.stats.csv \
	                                     gP.lP.raw.links.stats.csv \
	                                     gP.lP.raw.string.stats.csv \
	                                     gS.lS.raw.links.stats.csv \
	                                     gS.lS.raw.string.stats.csv \
	                                     gS.lS.clean.links.stats.csv \
	                                     gS.lS.clean.string.stats.csv \
	                                     gE.lE.clean.links.stats.csv \
	                                     gE.lE.clean.string.stats.csv \
	                                     gS.lE.clean.links.stats.csv \
	                                     gS.lE.clean.string.stats.csv \
	                                     gE.lS.clean.links.stats.csv \
	                                     gE.lS.clean.string.stats.csv
	cd report && $(MAKE)

clean:
	rm -rf $(DIRS) {perf,stoch,stocherr}.{plain,links,string}.stats.csv

.force:

.PHONY: all clean plots .force
