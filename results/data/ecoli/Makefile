
SHELL:=/bin/bash -euo pipefail

CTXDIR=../../../
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

all: ecoli.fa

NC_000913.fna:
	curl -o NC_000913.fna 'https://www.ncbi.nlm.nih.gov/nuccore/U00096.3?report=fasta&log$=seqview&format=text'

ecoli.fa: NC_000913.fna
	$(DNACAT) -M <(echo EColiK12) $< > $@

.PHONY: all
