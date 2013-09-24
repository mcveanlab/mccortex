
all:

chr1.fa.gz:
	wget -q -0 - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz | zcat | tr -d 'N' | gzip -c > chr1.fa.gz

staph.fa.gz:
	wget -q -0 - 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=FN433596&format=fasta&style=raw' | tr -d 'N' | gzip -c > staph.fa.gz

clean:


.PHONY: all clean
