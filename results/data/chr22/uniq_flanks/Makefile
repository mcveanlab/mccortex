
# Makefile
# Isaac Turner
# 2014-09-05
# Using 1Mb of chr22 and error free reads with both perfect & stochstic coverage

SHELL:=/bin/bash

CTXDIR=../../../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
REF=../chr22_17M_18M.fa
OUT=chr22.1Mbp.uniq.fa

K0=GGTCGCACACAAATACTACGGGCATTGGATGCGACCATAAGTCTTGACAGGCTTTGTTCCCT
K1=ACAACTTAACCTGGAACTAGAACTAATTTATGAGCGAGCCAGAACAGGTAGTCTGAGGGAGT

K0s=$(shell echo $(K0) | awk '{print substr($$0,1,11)}')
K1s=$(shell echo $(K1) | awk '{print substr($$0,length($$0)-10,11)}')
K0r=$(shell echo $(K0s) | $(DNACAT) -P -r -)
K1r=$(shell echo $(K1s) | $(DNACAT) -P -r -)

all: $(OUT) check

$(OUT):
	(echo '>chr22_17M_18M.11bp.and.up.uniq.flanks'; \
	 echo $(K0); $(DNACAT) -P $(REF); echo $(K1);) | \
	  $(DNACAT) -F -w 80 - > $@

clean:
	rm -rf $(OUT)

check: $(OUT)
	@echo Command should only print two kmers:
	$(DNACAT) -P $(OUT) | grep -ioE '($(K0s)|$(K0r)|$(K1s)|$(K1r))' -
	@echo Check file statistics:
	$(DNACAT) -s $(OUT)

.PHONY: all clean check
