
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dna
K=9

TGTS=ref.txt \
     bad.txt good.fa good.fq \
     rand.fq fix.fq \
     ref.k$(K).ctx

all: $(TGTS)

clean:
	rm -rf $(TGTS)

ref.txt:
	echo AGACAGGCATGTAGAGTTTTTTTTTTGGCTTGCACGAGCGAGAACCCATCAA > $@
	echo AGACAGGCATGTAGAGCTGGTGCGGAGGCTTGCACGAGCGAGAACCCATCAA >> $@

bad.txt:
	echo gGACAGGCATGcAGAGCTGGTGCGnnGGCTTGCACGAGCGAGAACCCATagNtcatacagata > $@
	echo gGACAGGCATGcAGAGnnnnnnnnnnGGCTTGCACGAGCGAGAACCCATagNtcaaacaTagaga >> $@
	echo nnnnnnnnnnnnnnAnnnnnnnnnnnnnnTaannnnnnaaaaaaaaaaaaN >> $@

rand.fq:
	$(DNACAT) -n 40 -q > $@

ref.k$(K).ctx: ref.txt
	$(CTX) build -m 10M -k $(K) -s ref -1 $< $@
	$(CTX) check -q $@

# Correct reads
good.fa: bad.txt ref.k$(K).ctx
	$(CTX) correct -m 10M -F FASTA -1 bad.txt:good ref.k9.ctx
	gzip -d good.fa.gz
	cat ref.txt
	cat bad.txt
	cat good.fa

# input: plain output: fastq
good.fq: bad.txt ref.k$(K).ctx
	$(CTX) correct -m 10M -F FASTQ -1 bad.txt:good ref.k9.ctx
	gzip -d good.fq.gz
	cat ref.txt
	cat bad.txt
	cat good.fq

# input: fastq output: fastq
fix.fq: rand.fq ref.k$(K).ctx
	$(CTX) correct -m 10M -F FASTQ -1 rand.fq:fix ref.k9.ctx
	gzip -d fix.fq.gz
	cat ref.txt
	cat rand.fq
	cat fix.fq

.PHONY: all clean
