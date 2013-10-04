
# Usage: ../run-sample-pop.sh <genome.fa> <stampyhsh> <indivs> <ploidy> <kmer>
#                             <snps> <indels> <invs> <invlen> <readlen> <mpsize>
#                             <allelecovg> <error|noerror> <run|print>

all: smaller.stidx smaller.sthash smaller.fa.fai
	../run-sample-pop.sh smaller.fa smaller 10 2 31 10 10 10 10 100 200 30 error print > run.sh
	chmod +x run.sh
	./run.sh

smaller.stidx: smaller.fa
	stampy.py -G smaller smaller.fa
smaller.sthash: smaller.stidx
	stampy.py -g smaller -H smaller
smaller.fa.fai: smaller.fa
	samtools faidx smaller.fa

clean:
	rm -rf genomes reads graphs bubbles vcfs gap_sizes* mp_sizes*

.PHONY: all clean
