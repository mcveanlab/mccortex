
SHELL:=/bin/bash -euo pipefail

CTXDIR=../../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

all: ecoli.fa

ecoli.fa:
	wget -O - 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna' | \
	$(DNACAT) -M <(echo EColiK12) - > $@

.PHONY: all
