SHELL:=/bin/bash

all: confidence_plot.pdf confidence.pdf

confidence_plot.pdf: plot-confidence.R
	R --vanilla -f plot-confidence.R --args $@

confidence.pdf: confidence.tex confidence_plot.pdf .force
	pdflatex confidence
	pdflatex confidence

clean:
	rm -rf confidence_plot.pdf confidence.{pdf,aux,log}

.force:

.PHONY: all clean .force
