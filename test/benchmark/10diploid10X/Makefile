
# Usage: ../run-sample-pop.sh <genome.fa> <stampyhsh> <indivs> <ploidy> <kmer>
#                             <snps> <indels> <invs> <invlen> <readlen> <mpsize>
#                             <allelecovg> <error|noerror> <run|print>
# Set stampy to the command needed
STAMPY_BIN=~/bioinf/stampy-1.0.20/stampy.py

UNAME=$(shell uname -s)
ifeq ($(UNAME),Darwin)
	STAMPY=python2.6 $(STAMPY_BIN)
else
	STAMPY=python $(STAMPY_BIN)
endif

all: smaller.stidx smaller.sthash smaller.fa.fai
	../run-sample-pop.sh smaller.fa smaller 10 2 31 10 10 10 10 100 200 30 ../PhiX.100K.1.fq.gz print 10 10 > run.sh
	chmod +x run.sh
	./run.sh

smaller.stidx: smaller.fa
	$(STAMPY) -G smaller smaller.fa
smaller.sthash: smaller.stidx
	$(STAMPY) -g smaller -H smaller
smaller.fa.fai: smaller.fa
	samtools faidx smaller.fa

clean:
	rm -rf genomes reads graphs bubbles vcfs gap_sizes* mp_sizes*

.PHONY: all clean
