
# first 1Mbp of chr1
# 10000 SNPs (1%)
# 5000 indels (0.5%)
# 500 inversions of size 1-100bp (50 * 500 = 0.025%)

INPUT_SEQ=../chr21.1Mb.fa.gz
SNPS=10000
INDELS=5000
INV=500
INVLEN=100
READLEN=100
MPSIZE=450
ALLELECOVG=30

CTX=~/cortex/versions/current/bin/cortex_var_31_c1_s8
CTX2=../../../bin/ctx31

PROGS=$(CTX) $(CTX2)

all: $(INPUT_SEQ) $(PROGS)
	../sample_diploid.sh $(INPUT_SEQ) $(SNPS) $(INDELS) $(INV) $(INVLEN) $(READLEN) $(MPSIZE) $(ALLELECOVG)

$(CTX):
	cd ../../../../current; make NUM_SHADES=8

$(CTX2):
	cd ../../../; make

clean:
	rm -rf *.ctx *.ctp *.fa *.falist
	rm -rf diploid.* falsepos.k31*.vcf truth.k31*.vcf
	rm -rf gap_sizes.*.csv mp_sizes.*.csv

.PHONY: all clean
