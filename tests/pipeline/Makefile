SHELL=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
READSIM=$(CTXDIR)/libs/readsim/readsim
MKPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
VCFALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles
BWA=$(CTXDIR)/libs/bwa/bwa

SEQ_FILES=ref/ref.fa samples/ben.fa samples/jerry.fa \
          samples/ben.se.fa \
          samples/jerry.1.fa samples/jerry.2.fa

DIRS=samples ref

all: job.k29.k31.mk $(SEQ_FILES)
	$(MAKE) -f job.k29.k31.mk CTXDIR=$(CTXDIR) MEM=200M NKMERS=1M all
	$(BCFTOOLS) isec proj/vcfs/breakpoints.joint.links.k29.k31.vcf.gz proj/vcfs/bubbles.joint.links.k29.k31.vcf.gz -p isec
	@A=`$(VCFALLELES) isec/0000.vcf`; \
	 B=`$(VCFALLELES) isec/0001.vcf`; \
	 C=`$(VCFALLELES) isec/0002.vcf`; \
	(( $$A == 0 )) || ( echo "Bubbles missed $$A" && false; ); \
	(( $$B == 0 )) || ( echo "Breakpoint missed $$B" && false; ); \
	(( $$C == 6 )) || ( echo "Both callers missed some $$C != 6" && false; );
	@echo 'All correct!'

# 600 random bases generated with: dnacat -n 600
#.............................................................................V................................................V.......................................*...........................................................................................................VV.................................................................................................................................................*....................................................................................................................VVVV..........................................................
ref/ref.fa: | $(DIRS)
	( echo '>ref'; echo \
'GCGTTACAATATCGTATTGGGTTCGTGACCAACACTCCCATTTCTTGATATGACGCCATCAAACGACTACACGGAGACCCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCA-CCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTCCTTAAGATCAGGAACTGGGGTGTACGACGGCCTCGCTGACACGGTACCAGCCGTGCACCGATGCTGCTAGGCACCCGTCGCCTGCTCAAGAAATGGCTGGGTTCAATAGAAGCGTTTGTGAGTGCTTCGACTCGTTAGGATGTAATTAGGGCCAGTAGTCAACCAGCGCTAGTGAGAATATGATAGAGATTTCGCAAAGTCCTTGGTATACAGGATCTCAACCCACAGACTGCGGAGGCTGTGGTGCCATCATCGGACTCACTACGTCCTTGTCAGGCCTAACCTTTCAGGGCGGCAAGCTACGGTTACCTGACCGAAGTCTTATTCACAGTTCGGTAGCTCCAATCATTGCGAGGTTAGCTTAACGCCTGACATTACCTGGCAAACATGCTCCTTTCACGACCGTTTATCGGCGCGATTTGATATCCACTTG' | tr -d '-'; ) > $@

samples/ben.fa: | $(DIRS)
	( echo '>ben'; echo \
'GCGTTACAATATCGTATTGGGTTCGTGACCAACACTCCCATTTCTTGATATGACGCCATCAAACGACTACACGGAGA-CCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCAACCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTCCTTAAGATCAGGAACTGGGGTGTACGACGGCCTCGCTGACACGGTACCAGCCGTGCACCGATGCTGCTAGGCACCCGTCGCCTGCTCAAGAAATGGCTGGGTTCAAT--AAGCGTTTGTGAGTGCTTCGACTCGTTAGGATGTAATTAGGGCCAGTAGTCAACCAGCGCTAGTGAGAATATGATAGAGATTTCGCAAAGTCCTTGGTATACAGGATCTCAACCCACAGACTGCGGAGGCTGTGGTGCCATCATCGGACTCACTACGTCCTTGTCAGGCCTAACCTTTCAGGGCGGCAAGCTACGGTTACCTGACCGAAGTCTTATTCACAGTTCGGTAGCTCCAATCATTGCGAGGTTAGCTTAACGCCTGACATTACCTGGCAAACATGCTCCTTTCACGACCGTTTATCGGCGCGATTTGATATCCACTTG' | tr -d '-'; ) > $@

samples/jerry.fa:
	( echo '>jerry'; echo \
'GCGTTACAATATCGTATTGGGTTCGTGACCAACACTCCCATTTCTTGATATGACGCCATCAAACGACTACACGGAGACCCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCAACCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTaCTTAAGATCAGGAACTGGGGTGTACGACGGCCTCGCTGACACGGTACCAGCCGTGCACCGATGCTGCTAGGCACCCGTCGCCTGCTCAAGAAATGGCTGGGTTCAATAGAAGCGTTTGTGAGTGCTTCGACTCGTTAGGATGTAATTAGGGCCAGTAGTCAACCAGCGCTAGTGAGAATATGATAGAGATTTCGCAAAGTCCTTGGTATACAGGATCTCAACCCACAGACTGCGGAGGCTGTGGTGCCATCATCtGACTCACTACGTCCTTGTCAGGCCTAACCTTTCAGGGCGGCAAGCTACGGTTACCTGACCGAAGTCTTATTCACAGTTCGGTAGCTCCAATCATTGCGAGGTTAGCTTAACGCCTG----TACCTGGCAAACATGCTCCTTTCACGACCGTTTATCGGCGCGATTTGATATCCACTTG' | tr -d '-'; ) > $@

ref/ref.fa.bwt: ref/ref.fa
	$(BWA) index ref/ref.fa

bams: samples/ben.bam samples/jerry.bam

sams: samples/ben.sam samples/jerry.sam

samples/ben.sam: samples/ben.se.fa ref/ref.fa
	$(BWA) mem ref/ref.fa $< > $@

samples/jerry.sam: samples/jerry.1.fa samples/jerry.2.fa ref/ref.fa.bwt
	$(BWA) mem ref/ref.fa samples/jerry.1.fa samples/jerry.2.fa > $@

%.bam: %.sam
	samtools sort -O bam -T tmpsort $< > $@
	samtools index $@

%.se.fa: %.fa
	$(READSIM) -d 50 -l 100 -s -e 0.01 -r $< tmp
	gzip -fcd tmp.fa.gz > $@
	rm -rf tmp.fa.gz

# 400 bp fragment len, 100 bp reads
%.1.fa %.2.fa: %.fa
	$(READSIM) -d 50 -l 100 -i 75 -e 0.01 -r $< tmp
	gzip -fcd tmp.1.fa.gz > $*.1.fa
	gzip -fcd tmp.2.fa.gz > $*.2.fa
	rm -rf tmp.{1,2}.fa.gz

samples.txt:
	( echo "ben samples/ben.se.fa . ."; echo "jerry . samples/jerry.1.fa:samples/jerry.2.fa" ) > $@

job.k29.k31.mk: samples.txt
	$(MKPIPELINE) -r ref/ref.fa 29,31 proj $< > $@

$(DIRS):
	mkdir -p $@

clean:
	rm -rf samples.txt job.k29.k31.mk proj ref samples tmp.* isec

.PHONY: all clean sams bams
