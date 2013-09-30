
all: chr1.fa.gz staph.fa.gz

HGCHR1=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
STAPH=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=FN433596&format=fasta&style=raw

chr1.fa.gz:
	wget -q -0 - '$(HGCHR1)' | zcat | tr -d 'N' > chr1.fa
	gzip -c chr1.fa > chr1.fa.gz
	samtools faidx chr1.fa
	samtools faidx chr1.fa.gz

staph.fa.gz:
	wget -q -0 - '$(STAPH)' | tr -d 'N' | > staph.fa
	gzip -c staph.fa > staph.fa.gz
	samtools faidx staph.fa
	samtools faidx staph.fa.gz

clean:
	rm staph.fa* chr1.fa*

.PHONY: all clean
