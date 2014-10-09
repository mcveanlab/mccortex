
SAMTOOLS=samtools

all: chr22.fa chr22_17M_18M.fa

chr22.fa:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
	gzip -d chr22.fa.gz

chr22_17M_18M.fa: chr22.fa
	$(SAMTOOLS) faidx chr22.fa chr22:17000000-17999999 | tr -d 'n' | tr -d 'N' > $@
