
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

BINARIES=diploid.k31.ctx diploid.shaded.k31.ctx
BUBBLES=diploid.oldbc.bubbles diploid.newbc.bubbles.gz diploid.newbc.shaded.bubbles.gz
FLANKS=diploid.oldbc.5pflanks.fa.gz diploid.newbc.5pflanks.fa.gz diploid.newbc.shaded.5pflanks.fa.gz
READFILES=reads0.0.fa reads0.1.fa reads1.0.fa reads1.1.fa reads.0.falist reads.1.falist
GENOMEFILES=genome0.fa genome1.fa mask0.fa mask1.fa
VCFS=truth.k31.vcf diploid.oldbc.vcf diploid.newbc.vcf diploid.newbc.shaded.vcf

CTX=~/cortex/versions/current/bin/cortex_var_31_c1_s8
CTX2=../../../bin/ctx_call
UNIQ=../../../bin/ctx_unique
PLACE=../../../bin/ctx_place

PROGS=$(CTX) $(CTX2) $(UNIQ) $(PLACE)

all: $(INPUT_SEQ) $(PROGS)
	../sample_diploid.sh $(INPUT_SEQ) $(SNPS) $(INDELS) $(INV) $(INVLEN) $(READLEN) $(MPSIZE) $(ALLELECOVG)

$(CTX):
	cd ../../../../current; make NUM_SHADES=8

$(CTX2):
	cd ../../../; make ctx_call

$(UNIQ):
	cd ../../../; make ctx_unique

$(PLACE):
	cd ../../../; make ctx_place

clean:
	rm -rf $(BINARIES) $(BUBBLES) $(FLANKS) $(READFILES) $(GENOMEFILES) $(VCFS)
	rm -rf gap_sizes.*.csv mp_sizes.*.csv
	rm -rf falsepos.k31.newbc.vcf falsepos.k31.oldbc.vcf falsepos.k31.shaded.vcf
	rm -rf truth.k31.newbc.vcf truth.k31.oldbc.vcf truth.k31.shaded.vcf

.PHONY: all clean
