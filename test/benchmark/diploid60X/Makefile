
# first 1Mbp of chr1
# 10000 SNPs (1%)
# 5000 indels (0.5%)
# 500 inversions of size 1-100bp (50 * 500 = 0.025%)

INPUT_SEQ=../chr21.1Mb.fa
STAMPY_HSH=../chr21
NINDIVS=1
PLOIDY=2
KMER=31
NSNPS=10000
NINDELS=5000
NINV=500
INVLEN=100
READLEN=100
MPSIZE=250
ALLELECOVG=30
ERROR=error
#error or noerror

# Set stampy to the command needed
STAMPY_BIN=~/bioinf/stampy-1.0.20/stampy.py

UNAME=$(shell uname -s)
ifeq ($(UNAME),Darwin)
	STAMPY=python2.6 $(STAMPY_BIN)
else
	STAMPY=python $(STAMPY_BIN)
endif

all: ../chr21.1Mb.fa ../chr21.1Mb.fa.fai ../chr21.stidx ../chr21.sthash
	../run-sample-pop.sh $(INPUT_SEQ) $(STAMPY_HSH) $(NINDIVS) $(PLOIDY) $(KMER) $(NSNPS) $(NINDELS) $(NINV) $(INVLEN) $(READLEN) $(MPSIZE) $(ALLELECOVG) $(ERROR) print > run.sh
	chmod +x run.sh
	./run.sh

../chr21.1Mb.fa:
	gzip -c -d ../chr21.1Mb.fa.gz > ../chr21.1Mb.fa
../chr21.1Mb.fa.fai:
	samtools faidx ../chr21.1Mb.fa

../chr21.stidx: ../chr21.1Mb.fa
	$(STAMPY) -G ../chr1 ../chr21.1Mb.fa
../chr21.sthash: ../chr21.stidx
	$(STAMPY) -g ../chr1 -H ../chr21

clean:
	rm -rf genomes reads graphs bubbles vcfs gap_sizes* mp_sizes*

.PHONY: all clean
