SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
K=9

TGTS=ref.txt \
     bad.txt good.fa good.fq \
     rand.fq fix.fq \
     indels.bad.fq indels.good.fq \
     ref.k$(K).ctx

all: $(TGTS)

clean:
	rm -rf $(TGTS) good.fa.gz good.fq.gz fix.fq.gz indels.good.fq.gz

ref.txt:
	echo AGACAGGCATGTAGAGTTTTTTTTTTGGCTTGCACGAGGGAGAACCCATCAA > $@
	echo AGACAGGCATGTAGAGCTGGTGCGGAGGCTTGCACGAGGGAGAACCCATCAA >> $@

bad.txt:
	echo gGACAGGCATGcAGAGCTGGTGCGnnGGCTTGCACGAGGGAGAACCCATagNtcatacagata > $@
	echo gGACAGGCATGcAGAGnnnnnnnnnnGGCTTGCACGAGGGAGAACCCATagNtcaaacaTagaga >> $@
	echo nnnnnnnnnnnnnnAnnnnnnnnnnnnnnTaannnnnnaaaaaaaaaaaaN >> $@

indels.bad.fq:
	echo '@ins_T>TT' > $@
	echo AGACAGGCATGTAGAGCTGGTtGCGGAGGCTTGCACGAGGGAGAACCCATCAA >> $@
	echo + >> $@
	echo abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1 >> $@
	echo '@ins_T>TT_and_SNP_A>C' >> $@
	echo AGACAGGCATGTAGAGCTGGTtGCGGAGGCTTGCACGcGGGAGAACCCATCAA >> $@
	echo + >> $@
	echo abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1 >> $@
	echo '@ins_T>TTT_and_SNP_g>A' >> $@
	echo AGACAGGCATGTAGAGCTGGTttGCGGAGGCTTaCACGAGGGAGAACCCATCAA >> $@
	echo + >> $@
	echo abcdefghijklmnopqrstuvvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1 >> $@
	echo '@del_GG>G_and_extra' >> $@
	echo AGACAGGCATGTAGAGCTgTGCGGAGGCTTGCACGAGGGAGAACCCATgattgcaattgaa >> $@
	echo + >> $@
	echo abcdefghijklmnopqrsuwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ123456789ab >> $@
	echo '@extra_del_GG>G_and_GGG>G' >> $@
	echo accttgttACAGGCATGTAGAGCTgTGCGGAGGCTTGCACGAGAGAACCCATCAA >> $@
	echo + >> $@
	echo 654321abcdefghijklmnopqrsuwxyzABCDEFGHIJKLMPQRSTUVWXYZ1 >> $@

rand.fq:
	$(DNACAT) -n 40 -Q > $@

ref.k$(K).ctx: ref.txt
	$(CTX) build -m 10M -k $(K) -s ref -1 $< $@
	$(CTX) check -q $@

# Correct reads
# -P, --print-orig means add original sequence to read name 'orig=SEQ'
good.fa: bad.txt ref.k$(K).ctx
	$(CTX) correct -t 1 -m 10M -F FASTA --print-orig -1 bad.txt:good ref.k$(K).ctx
	gzip -d good.fa.gz
	cat ref.txt
	cat bad.txt
	cat good.fa

# input: plain output: fastq
good.fq: bad.txt ref.k$(K).ctx
	$(CTX) correct -t 1 -m 10M -F FASTQ --print-orig -1 bad.txt:good ref.k$(K).ctx
	gzip -d good.fq.gz
	cat ref.txt
	cat bad.txt
	cat good.fq

# input: fastq output: fastq
fix.fq: rand.fq ref.k$(K).ctx
	$(CTX) correct -t 1 -m 10M -F FASTQ --print-orig -1 rand.fq:fix ref.k$(K).ctx
	gzip -d fix.fq.gz
	cat ref.txt
	cat rand.fq
	cat fix.fq

indels.good.fq: indels.bad.fq
	$(CTX) correct -t 1 -m 10M -F FASTQ --print-orig -1 $<:indels.good ref.k$(K).ctx
	gzip -d indels.good.fq.gz
	@echo == in ==
	cat $<
	@echo == out ==
	cat $@

# Plots to help understand what is going on
plots: indel.AA.pdf snp.AT.pdf

indel.AA.pdf:
	printf 'CTGTTCCAAGAGTAACGTTA\nCTGTTCCAAGAAAGTAACGTTA\n' | \
	$(CTXDIR)/scripts/seq2pdf.sh $(K) - > $@

snp.AT.pdf:
	printf 'CTGTTCCAAGAGTAACGTTA\nCTGTTCCAAGTGTAACGTTA\n' | \
	$(CTXDIR)/scripts/seq2pdf.sh $(K) - > $@

.PHONY: all clean plots
