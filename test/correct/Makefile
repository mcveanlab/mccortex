
CTX=../../bin/ctx31
SEQRND=../../libs/seq_file/bin/seqrnd
K=9

TGTS=ref.fa good.fa bad.fa ref.k$(K).ctx

all: $(TGTS)

clean:
	rm -rf $(TGTS)

ref.fa:
	echo AGACAGGCATGTAGAGTTTTTTTTTTGGCTTGCACGAGCGAGAACCCATCAA > $@
	echo AGACAGGCATGTAGAGCTGGTGCGGAGGCTTGCACGAGCGAGAACCCATCAA >> $@

bad.fa:
	echo gGACAGGCATGcAGAGCTGGTGCGnnGGCTTGCACGAGCGAGAACCCATagNtcatacagata > $@
	echo gGACAGGCATGcAGAGnnnnnnnnnnGGCTTGCACGAGCGAGAACCCATagNtcaaacaTagaga >> $@
	echo nnnnnnnnnnnnnnAnnnnnnnnnnnnnnTaannnnnnaaaaaaaaaaaaN >> $@

ref.k$(K).ctx: ref.fa
	$(CTX) build -m 10M -k $(K) -s ref -1 ref.fa ref.k9.ctx

# Correct reads
good.fa: bad.fa ref.k$(K).ctx
	$(CTX) correct -m 10M -1 bad.fa:good ref.k9.ctx
	gzip -d good.fa.gz
	cat ref.fa
	cat bad.fa
	cat good.fa

.PHONY: all clean
