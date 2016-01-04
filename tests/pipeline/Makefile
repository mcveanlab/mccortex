SHELL=/bin/bash -euo pipefail

CTXDIR=../..
MCCORTEX=$(CTXDIR)/bin/mccortex31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
READSIM=$(CTXDIR)/libs/readsim/readsim
MKPIPELINE=$(CTXDIR)/scripts/make-pipeline.pl
MKREPORT=$(CTXDIR)/scripts/make-report.pl
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools
BGZIP=$(CTXDIR)/libs/htslib/bgzip
VCFALLELES=$(CTXDIR)/libs/biogrok/vcf-count-alleles
VCFENTRIES=$(CTXDIR)/libs/biogrok/vcf-entries
BWA=$(CTXDIR)/libs/bwa/bwa

SEQ_FILES=ref/ref.fa samples/ben.fa samples/jerry.fa \
          samples/ben.se.fa \
          samples/jerry.1.fa samples/jerry.2.fa

CALLS=proj/vcfs/bubbles.joint.links.k29.k31.geno.vcf.gz \
      proj/vcfs/breakpoints.joint.links.k29.k31.geno.vcf.gz

DIRS=samples ref

DEPTH=50
ERRRATE=0.01
READLEN=100

all: $(CALLS) truth.vcf.gz report
	$(BCFTOOLS) isec truth.vcf.gz proj/vcfs/bubbles.joint.links.k29.k31.geno.vcf.gz -p isec_bub
	$(BCFTOOLS) isec truth.vcf.gz proj/vcfs/breakpoints.joint.links.k29.k31.geno.vcf.gz -p isec_brk
	@BUB_T=`$(VCFALLELES) isec_bub/0000.vcf`; \
	 BUB_C=`$(VCFALLELES) isec_bub/0001.vcf`; \
	 BRK_T=`$(VCFALLELES) isec_brk/0000.vcf`; \
	 BRK_C=`$(VCFALLELES) isec_brk/0001.vcf`; \
	(( $$BUB_T == 0 && $$BRK_T == 0 )) || ( echo "Bubbles missed $$BUB_T, Breakpoints missed $$BRK_T" && false; ); \
	(( $$BUB_C == 0 )) || ( echo "Bubbles called $$BUB_C extra"; ); \
	(( $$BRK_C == 0 )) || ( echo "Breakpoints called $$BRK_C extra"; );
	@echo 'All correct!'

$(CALLS): job.k29.k31.mk $(SEQ_FILES)
	$(MAKE) -f job.k29.k31.mk CTXDIR=$(CTXDIR) MEM=50M NKMERS=1M BRK_REF_KMERS=10 vcfs

report: report/Makefile

report/Makefile: $(CALLS)
	$(MKREPORT) proj report

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

truth.vcf.gz:
	(echo '##fileformat=VCFv4.1'; \
	 echo '##fileDate='`date "+%Y%m%d"`; \
	 echo '##reference=ref/ref.fa'; \
	 echo '##contig=<ID=ref,length=599>'; \
	 echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'; \
	 echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ben	jerryT'; \
	 echo 'ref	77	id0	AC	A	.	PASS	.	GT	1	0'; \
	 echo 'ref	125	id1	C	CA	.	PASS	.	GT	1	1'; \
	 echo 'ref	166	id2	C	A	.	PASS	.	GT	0	1'; \
	 echo 'ref	273	id3	TAG	T	.	PASS	.	GT	1	0'; \
	 echo 'ref	421	id4	G	T	.	PASS	.	GT	0	1'; \
	 echo 'ref	537	id5	GACAT	G	.	PASS	.	GT	0	1';) > truth.vcf
	$(BGZIP) truth.vcf
	$(BCFTOOLS) index truth.vcf.gz

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
	$(READSIM) -d $(DEPTH) -l $(READLEN) -s -e $(ERRRATE) -r $< tmp
	gzip -fcd tmp.fa.gz > $@
	rm -rf tmp.fa.gz

# 400 bp fragment len, 100 bp reads
%.1.fa %.2.fa: %.fa
	$(READSIM) -d $(DEPTH) -l $(READLEN) -i 75 -e $(ERRRATE) -r $< tmp
	gzip -fcd tmp.1.fa.gz > $*.1.fa
	gzip -fcd tmp.2.fa.gz > $*.2.fa
	rm -rf tmp.{1,2}.fa.gz

samples.txt:
	( echo "ben  samples/ben.se.fa  .  ."; echo "jerry .   samples/jerry.1.fa:samples/jerry.2.fa" ) > $@

job.k29.k31.mk: samples.txt
	$(MKPIPELINE) -r ref/ref.fa --ploidy 1 29,31 proj $< > $@

$(DIRS):
	mkdir -p $@

VIEWFILES=truth.vcf.gz \
          proj/vcfs/bubbles.joint.links.k29.k31.geno.vcf.gz \
          proj/vcfs/breakpoints.joint.links.k29.k31.geno.vcf.gz

view: $(VIEWFILES)
	$(VCFENTRIES) truth.vcf.gz
	$(VCFENTRIES) proj/vcfs/bubbles.joint.links.k29.k31.geno.vcf.gz
	$(VCFENTRIES) proj/vcfs/breakpoints.joint.links.k29.k31.geno.vcf.gz

clean:
	rm -rf samples.txt job.k29.k31.mk proj ref samples isec_bub isec_brk
	rm -rf tmp.* truth.* report

.PHONY: all clean sams bams view report
