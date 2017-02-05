
SHELL=/bin/bash -eou pipefail

SAMTOOLS=samtools
DNACAT=../../../libs/seq_file/bin/dnacat

all: chr22.fa chr22_17M_18M.fa chr22_28M_29M.fa

chr22.fa:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
	gzip -d chr22.fa.gz

chr22_17M_18M.fa: chr22.fa
	$(SAMTOOLS) faidx chr22.fa chr22:17000000-17999999 | \
	$(DNACAT) -M <(echo 'chr22_17M_18M') - > $@
	if $(DNACAT) -P $@ | grep -qi n; then false; fi

chr22_28M_29M.fa: chr22.fa
	$(SAMTOOLS) faidx chr22.fa chr22:28000000-28999999 | \
	$(DNACAT) -M <(echo 'chr22_28M_29M') - > $@
	if $(DNACAT) -P $@ | grep -qi n; then false; fi
