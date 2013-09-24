
# first 1Mbp of chr1
# 10000 SNPs (1%)
# 5000 indels (0.5%)
# 500 inversions of size 1-100bp (50 * 500 = 0.025%)

INPUT_SEQ=../chr21.1Mb.fa.gz
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

CTX=~/cortex/versions/current/bin/cortex_var_31_c1_s8
CTX2=../../../bin/ctx31

PROGS=$(CTX) $(CTX2)

all: $(INPUT_SEQ) $(PROGS)
	../run-sample-pop.sh $(INPUT_SEQ) $(NINDIVS) $(PLOIDY) $(KMER) $(NSNPS) $(NINDELS) $(NINV) $(INVLEN) $(READLEN) $(MPSIZE) $(ALLELECOVG) $(ERROR)

$(CTX):
	cd ../../../../current; make NUM_SHADES=8

$(CTX2):
	cd ../../../; make

clean:
	rm -rf diploid* falsepos* reads* mask* genome* truth* gap_sizes* mp_sizes* pop.*

.PHONY: all clean
