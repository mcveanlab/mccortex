HGCHR1=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
STAPH=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=FN433596&format=fasta&style=raw
FACAT=../../../libs/seq_file/bin/facat

# PhiX.1.fq.gz needs to be downloaded
PROFILE_READS=PhiX.1.fq.gz

# Set stampy to the command needed
STAMPY_BIN=~/bioinf/stampy-1.0.20/stampy.py

UNAME=$(shell uname -s)
ifeq ($(UNAME),Darwin)
	STAMPY=python2.6 $(STAMPY_BIN)
else
	STAMPY=python $(STAMPY_BIN)
endif

all: chr1/chr1.fa.gz staph/staph.fa.gz run-files $(PROFILE_READS)

chr1/chr1.orig.fa.gz:
	mkdir -p chr1
	wget -O chr1/chr1.orig.fa.gz '$(HGCHR1)'

staph/staph.orig.fa:
	mkdir -p staph
	wget -O staph/staph.orig.fa '$(STAPH)'

chr1/chr1.fa.gz: chr1/chr1.orig.fa.gz
	zcat chr1/chr1.orig.fa.gz | tr -d 'N' | $(FACAT) -w 50 > chr1/chr1.fa
	gzip -c chr1/chr1.fa > chr1/chr1.fa.gz
	samtools faidx chr1/chr1.fa
	samtools faidx chr1/chr1.fa.gz
	$(STAMPY) -G chr1/chr1 chr1/chr1.fa.gz
	$(STAMPY) -g chr1/chr1 -H chr1/chr1

staph/staph.fa.gz: staph/staph.orig.fa
	cat staph/staph.orig.fa | tr -d 'N' | $(FACAT) -w 50 > staph/staph.fa
	gzip -c staph/staph.fa > staph/staph.fa.gz
	samtools faidx staph/staph.fa
	samtools faidx staph/staph.fa.gz
	$(STAMPY) -G staph/staph staph/staph.fa.gz
	$(STAMPY) -g staph/staph -H staph/staph

# Number of events [hs][SI][OM]
# h=human s=staph S=SNPs I=indels O=one M=many
hSO=225000
hSM=798242
hIO=25000
hIM=88694

sSO=27000
sSM=202081
sIO=3000
SIM=22453

run-files:
	mkdir -p OneHumanErrorReads100   OneStaphErrorReads100   OneHumanErrorReads250   OneStaphErrorReads250
	mkdir -p OneHumanNoErrorReads100 OneStaphNoErrorReads100 OneHumanNoErrorReads250 OneStaphNoErrorReads250
	mkdir -p ManyHumanErrorReads100  ManyStaphErrorReads100  ManyHumanErrorReads250  ManyStaphErrorReads250
	#
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    1 2 31 $(hSO) $(hIO) 0 0 100 250 30  ../PhiX.1.fq.gz print 40 21 > OneHumanErrorReads100/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 2 1 31 $(sSO) $(sIO) 0 0 100 250 100 ../PhiX.1.fq.gz print 40 21 > OneStaphErrorReads100/run.sh
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    1 2 31 $(hSO) $(hIO) 0 0 250 250 30  ../PhiX.1.fq.gz print 40 21 > OneHumanErrorReads250/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 2 1 31 $(sSO) $(sIO) 0 0 250 250 100 ../PhiX.1.fq.gz print 40 21 > OneStaphErrorReads250/run.sh
	#
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    1 2 31 $(hSO) $(hSO) 0 0 100 250 30  noerror print 40 21 > OneHumanNoErrorReads100/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 2 1 31 $(sSO) $(sSO) 0 0 100 250 100 noerror print 40 21 > OneStaphNoErrorReads100/run.sh
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    1 2 31 $(hSO) $(hSO) 0 0 250 250 30  noerror print 40 21 > OneHumanNoErrorReads250/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 2 1 31 $(sSO) $(sSO) 0 0 250 250 100 noerror print 40 21 > OneStaphNoErrorReads250/run.sh
	#
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    10   2 31 $(hSM) $(hSM) 0 0 100 250 30  ../PhiX.1.fq.gz print 40 21 > ManyHumanErrorReads100/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 1000 1 31 $(sSM) $(sSM) 0 0 100 250 100 ../PhiX.1.fq.gz print 40 21 > ManyStaphErrorReads100/run.sh
	../run-sample-pop.sh ../chr1/chr1.fa  ../chr1/chr1    10   2 31 $(hSM) $(hSM) 0 0 250 250 30  ../PhiX.1.fq.gz print 40 21 > ManyHumanErrorReads250/run.sh
	../run-sample-pop.sh ../staph/staph.fa ../staph/staph 1000 1 31 $(sSM) $(sSM) 0 0 250 250 100 ../PhiX.1.fq.gz print 40 21 > ManyStaphErrorReads250/run.sh

clean:
	rm -rf staph.fa* chr1.fa* staph.st* chr1.st*

.PHONY: all clean run-files
