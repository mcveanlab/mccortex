SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNA=$(CTXDIR)/libs/seq_file/bin/dnacat
MCCORTEX=$(CTXDIR)/bin/mccortex31
K=9

READS=reads.fa reads.1.fa.gz reads.2.fa.gz reads.interleaved.fq.gz
TGTS=seq.fa seq.k$(K).ctx $(READS)
RESULTS=out/se.fq.gz out/se.1.fq.gz out/se.2.fq.gz \
        out/pe.fq.gz out/pe.1.fq.gz out/pe.2.fq.gz \
				out/ipe.fq.gz out/ipe.1.fq.gz out/ipe.2.fq.gz \
        out/pe.fa.gz out/pe.1.fa.gz out/pe.2.fa.gz \

all: $(TGTS) $(RESULTS)

seq.fa:
	echo ACGTTATTTAATCTGGTTACCGCCAGGTCAGGGCTATATGTGTAGACGAT > $@

reads.fa:
	printf \
'>r1 hit       |\n'\
'TTTTCGATTTAATCTCTATC\n'\
'>r2\n'\
'TGGCCGATGCAGCAGACGCC\n'\
'>r3\n'\
'GGGGGCCAGACCTCTAACAT\n'\
'>r4 hit  |\n'\
'TTATATGTGTACTGAGCCTC\n'\
'>r5\n'\
'TTGTGACCTCCTTTAGCCTG\n' > $@
reads.1.fa.gz:
	printf \
'>r1/1\n'\
'GCGAGTGGAACAGACGTTGA\n'\
'>r2/1 hit    |\n'\
'TAATCCCGCCAGGTACTCTC\n'\
'>r3/1\n'\
'CTTTGCGGAGCGATGAATGT\n'\
'>r4/1 hit     |\n'\
'GGGCACATCTGGTTACAGGG\n'\
'>r5/1\n'\
'GAGAAAACACAATCCCGAAC\n' | gzip -c > $@
reads.2.fa.gz:
	printf \
'>r1/2 hit          |\n'\
'GCGAGTGGAACGTAGACGAT\n'\
'>r2/2 hi|t\n'\
'ACGTTATTTGGACGACTCTA\n'\
'>r3/2\n'\
'CTTTGCGGAGCGATGAATGA\n'\
'>r4/2\n'\
'GGGCACATGGCCTCTCAGGC\n'\
'>r5/2\n'\
'GAGAAAACACAATCCCGAAG\n' | gzip -c > $@
reads.interleaved.fq.gz:
	$(DNA) --interleave reads.1.fa.gz reads.2.fa.gz | \
	  cat - <(printf '>hit\nTACCGCCAGGTCAGGGCT\n>moo\nACA') | gzip -c >$@

seq.k$(K).ctx: seq.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample Seq --seq $< $@

out/se.fq.gz: out/se.1.fq.gz
out/se.1.fq.gz: out/se.2.fq.gz
out/se.2.fq.gz: seq.k$(K).ctx $(READS)
	$(MCCORTEX) reads --seq reads.fa:out/se \
	             --seq reads.1.fa.gz:out/se.1 \
	             --seq reads.2.fa.gz:out/se.2 seq.k$(K).ctx

out/pe.fq.gz: out/pe.1.fq.gz
out/pe.1.fq.gz: out/pe.2.fq.gz
out/pe.2.fq.gz: seq.k$(K).ctx $(READS)
	$(MCCORTEX) reads --seq2 reads.1.fa.gz:reads.2.fa.gz:out/pe seq.k$(K).ctx

out/ipe.fq.gz: out/ipe.1.fq.gz
out/ipe.1.fq.gz: out/ipe.2.fq.gz
out/ipe.2.fq.gz: seq.k$(K).ctx reads.interleaved.fq.gz
	$(MCCORTEX) reads --seqi reads.interleaved.fq.gz:out/ipe seq.k$(K).ctx

out/pe.fa.gz: out/pe.1.fa.gz
out/pe.1.fa.gz: out/pe.2.fa.gz
out/pe.2.fa.gz: seq.k$(K).ctx $(READS)
	$(MCCORTEX) reads --format fa --seq2 reads.1.fa.gz:reads.2.fa.gz:out/pe seq.k$(K).ctx

clean:
	rm -rf $(TGTS) out

.PHONY: all clean
