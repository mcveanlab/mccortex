
all: kmer_agnostic.pdf

kmer_agnostic.pdf: kmer_agnostic.tex .force
	pdflatex kmer_agnostic
	pdflatex kmer_agnostic

clean:
	rm -rf kmer_agnostic.pdf kmer_agnostic.log kmer_agnostic.aux

.force:

.PHONY: all clean .force
