
all: kmer_agnostic.pdf

kmer_agnostic.pdf: kmer_agnostic.tex
	pdflatex kmer_agnostic
	pdflatex kmer_agnostic

clean:
	rm -rf kmer_agnostic.pdf kmer_agnostic.log kmer_agnostic.aux

.PHONY: all clean
