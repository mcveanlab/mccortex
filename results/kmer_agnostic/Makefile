SHELL:=/bin/bash -euo pipefail
# Makefile
# Isaac Turner
# 2014-10-01
# Using 1Mb of chr22 and error free reads to measure effect of kmer-size

# Expand these for all k values
DIRS=plots readlen80 readlen100 readlen125

LINK_100_CSVS=readlen100/k15/link_cleaning/link_cleaning.csv \
	            readlen100/k21/link_cleaning/link_cleaning.csv \
	            readlen100/k31/link_cleaning/link_cleaning.csv \
	            readlen100/k41/link_cleaning/link_cleaning.csv \
	            readlen100/k51/link_cleaning/link_cleaning.csv \
	            readlen100/k63/link_cleaning/link_cleaning.csv \
	            readlen100/k75/link_cleaning/link_cleaning.csv

KMER_100_CSVS=readlen100/k15/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k21/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k31/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k41/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k51/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k63/kmer_cleaning/kmer_cleaning.csv \
	            readlen100/k75/kmer_cleaning/kmer_cleaning.csv

RESULT_100_CSVS=readlen100/gP.l0.plain.stats.csv \
	              readlen100/gS.l0.plain.stats.csv \
	              readlen100/gE.l0.plain.stats.csv \
	              readlen100/gP.lP.raw.links.stats.csv \
	              readlen100/gP.lP.raw.string.stats.csv \
	              readlen100/gS.lS.raw.links.stats.csv \
	              readlen100/gS.lS.raw.string.stats.csv \
	              readlen100/gS.lS.clean.links.stats.csv \
	              readlen100/gS.lS.clean.string.stats.csv \
	              readlen100/gE.lE.clean.links.stats.csv \
	              readlen100/gE.lE.clean.string.stats.csv \
	              readlen100/gS.lE.clean.links.stats.csv \
	              readlen100/gS.lE.clean.string.stats.csv \
	              readlen100/gE.lS.clean.links.stats.csv \
	              readlen100/gE.lS.clean.string.stats.csv

RESULT_80_CSVS =$(subst readlen100,readlen80, $(RESULT_100_CSVS))
RESULT_125_CSVS=$(subst readlen100,readlen125,$(RESULT_100_CSVS))

PLOTS_100=plots/readlen100.pdf plots/link_100.pdf plots/kmer_100.pdf
PLOTS=$(substr 100,80,$(PLOTS_100)) $(PLOTS_100) $(substr 100,125,$(PLOTS_100))

all: $(PLOTS)

$(RESULT_80_CSVS):
	mkdir -p readlen80  && cd readlen80  && ../run.sh 80

$(RESULT_100_CSVS):
	mkdir -p readlen100 && cd readlen100 && ../run.sh 100

$(RESULT_125_CSVS):
	mkdir -p readlen125 && cd readlen125 && ../run.sh 125

plots/readlen80.pdf: $(RESULT_80_CSVS) | $(DIRS)
plots/readlen100.pdf: $(RESULT_100_CSVS) | $(DIRS)
plots/readlen125.pdf: $(RESULT_125_CSVS) | $(DIRS)

plots/kmer_100.pdf: $(KMER_100_CSVS) | $(DIRS)
	R --vanilla -f plot-link-err.R --args $@ $(KMER_100_CSVS)

plots/link_100.pdf: $(LINK_100_CSVS) | $(DIRS)
	R --vanilla -f plot-kmer-err.R --args $@ $(LINK_100_CSVS)

plots/readlen%.pdf:
	R --vanilla -f plot-results.R --args $@ \
	                                     readlen$*/gP.l0.plain.stats.csv \
	                                     readlen$*/gS.l0.plain.stats.csv \
	                                     readlen$*/gE.l0.plain.stats.csv \
	                                     readlen$*/gP.lP.raw.links.stats.csv \
	                                     readlen$*/gP.lP.raw.string.stats.csv \
	                                     readlen$*/gS.lS.raw.links.stats.csv \
	                                     readlen$*/gS.lS.raw.string.stats.csv \
	                                     readlen$*/gS.lS.clean.links.stats.csv \
	                                     readlen$*/gS.lS.clean.string.stats.csv \
	                                     readlen$*/gE.lE.clean.links.stats.csv \
	                                     readlen$*/gE.lE.clean.string.stats.csv \
	                                     readlen$*/gS.lE.clean.links.stats.csv \
	                                     readlen$*/gS.lE.clean.string.stats.csv \
	                                     readlen$*/gE.lS.clean.links.stats.csv \
	                                     readlen$*/gE.lS.clean.string.stats.csv

$(DIRS):
	mkdir -p $@

clean:
	rm -rf $(DIRS) {perf,stoch,stocherr}.{plain,links,string}.stats.csv

.force:

.PHONY: all clean plots .force
