
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

all: $(INPUT_SEQ)
	../run-sample-pop.sh $(INPUT_SEQ) $(STAMPY_HSH) $(NINDIVS) $(PLOIDY) $(KMER) $(NSNPS) $(NINDELS) $(NINV) $(INVLEN) $(READLEN) $(MPSIZE) $(ALLELECOVG) $(ERROR) print > run.sh
	chmod +x run.sh
	./run.sh

clean:
	rm -rf genomes reads graphs bubbles vcfs gap_sizes* mp_sizes*

.PHONY: all clean
