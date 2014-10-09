SHELL:=/bin/bash -euo pipefail
#
# Isaac Turner
# 2014-09-05
#

CTXDIR=../../../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
BWA=bwa
SAMTOOLS=samtools

REF=../chr22_17M_18M.fa

K0=GGTCGCACACAAATACTACGGGCATTGGATGCGACCATAAGTCTTGACAGGCTTTGTTCCCT
K1=ACAACTTAACCTGGAACTAGAACTAATTTATGAGCGAGCCAGAACAGGTAGTCTGAGGGAGT

K0s=$(shell echo $(K0) | awk '{print substr($$0,1,11)}')
K1s=$(shell echo $(K1) | awk '{print substr($$0,length($$0)-10,11)}')
K0r=$(shell echo $(K0s) | $(DNACAT) -P -r -)
K1r=$(shell echo $(K1s) | $(DNACAT) -P -r -)

all: chr22.1Mbp.uniq.fa chr22.1Mbp.uniq.fa.bwt chr22.1Mbp.uniq.fa.fai check

chr22.1Mbp.uniq.fa:
	(echo '>chr22_17M_18M.11bp.and.up.uniq.flanks'; \
	 echo $(K0); $(DNACAT) -P $(REF); echo $(K1);) | \
	  $(DNACAT) -F -w 80 - > $@

chr22.1Mbp.uniq.fa.bwt: chr22.1Mbp.uniq.fa
	$(BWA) index $<

chr22.1Mbp.uniq.fa.fai: chr22.1Mbp.uniq.fa
	$(SAMTOOLS) faidx $<

clean:
	rm -rf chr22.1Mbp.uniq.fa*

check: chr22.1Mbp.uniq.fa
	@echo; echo Command should only print two kmers:
	$(DNACAT) -P chr22.1Mbp.uniq.fa | grep -ioE '($(K0s)|$(K0r)|$(K1s)|$(K1r))' -
	@echo; echo Check file contains exactly 1Mbp+2*62:
	$(DNACAT) -s chr22.1Mbp.uniq.fa

.PHONY: all clean check
