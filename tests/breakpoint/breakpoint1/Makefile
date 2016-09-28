SHELL:=/bin/bash -euo pipefail

#
# chr1: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
# chr2: gCCGTAGGTAAGGGCGTTAGTGC
# chr3: GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG
# chr4: ATTCTACAGCAGGTCATGAGCAACCGGCACTCGAGCAGACGTACGGGAAA
#
# >mix
# gCCGTAGGTAAGGGCGTTAGaCGGGTTGGAGTTGGCCAAAGAAGTTCAAgg
# .1111111111111111111.3333333333333333333333333333..
#
# >repeat of 5bp
# CCCGTAGGTAAGGGCGTTAGTGCAGTGCAAGGCCACATTGGGACACGAGTTGAgg
# 111111111111111111111111xxxxx111111111111111111111111..
#
# >deletion of 1bp
# cAACACCCTTTGGTTTAAGCCG-GTTGGAGTTGGCCAAAGAAGTTCAAgg
# .333333333333333333333 3333333333333333333333333..
#
# >SNP
# ATTCTACAGCAGGTCATGAGCAACtGGCACTCGAGCAGACGTACGGGAAA
# 444444444444444444444444 4444444444444444444444444
#

K=11
CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
BRKCHCK=$(CTXDIR)/scripts/cortex_brkpnts_check_ref.pl
VCFSORT=$(CTXDIR)/libs/biogrok/vcf-sort
VCFRENAME=$(CTXDIR)/libs/biogrok/vcf-rename

BGZIP=$(CTXDIR)/libs/htslib/bgzip
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools

SEQS=sample.fa ref.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
TGTS=breakpoints.txt.gz breakpoints.norm.vcf.gz $(GRAPHS)
# join.k$(K).ctx

all: $(TGTS) cmp_breakpoint cmp_vcf

ref.fa:
	( echo '>chr1'; \
	  echo CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA; \
	  echo '>chr2'; \
	  echo gCCGTAGGTAAGGGCGTTAGTGC; \
	  echo '>chr3'; \
	  echo GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG; \
	  echo '>chr4'; \
	  echo ATTCTACAGCAGGTCATGAGCAACCGGCACTCGAGCAGACGTACGGGAAA; ) > $@

sample.fa:
	( echo '>contig0'; \
	  echo gCCGTAGGTAAGGGCGTTAGaCGGGTTGGAGTTGGCCAAAGAAGTTCAAgg; \
	  echo '>contig1'; \
	  echo CCCGTAGGTAAGGGCGTTAGTGCAGTGCAAGGCCACATTGGGACACGAGTTGAgg; \
	  echo '>contig2'; \
	  echo cAACACCCTTTGGTTTAAGCCGGTTGGAGTTGGCCAAAGAAGTTCAAgg; \
	  echo '>contig3'; \
	  echo ATTCTACAGCAGGTCATGAGCAACtGGCACTCGAGCAGACGTACGGGAAA; ) > $@

truth.vcf.gz:
	( printf "##fileformat=VCF4.1\n##fileDate="`date '+%Y%m%d'`"\n"; \
		printf "##reference=ref.fa\n"; \
		printf "##contig=<id=chr1,length=50>\n"; \
		printf "##contig=<id=chr2,length=23>\n"; \
		printf "##contig=<id=chr3,length=50>\n"; \
		printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"; \
		echo "chr1	18	.	T	TAGTGC	.	.	.	GT	1"; \
	  echo "chr3	21	.	CG	C	.	.	.	GT	1"; \
	  echo "chr4	25	.	C	T	.	.	.	GT	1"; ) > truth.vcf
	$(BGZIP) truth.vcf
	$(BCFTOOLS) index truth.vcf.gz

%.k$(K).ctx: %.fa
	$(MCCORTEX) build -m 10M -k $(K) --sample $* --seq $< $@ >& $@.log

breakpoints.txt.gz: sample.k$(K).ctx ref.fa
	$(MCCORTEX) breakpoints -t 1 -m 10M --minref 5 \
	                  --seq ref.fa --out $@ sample.k$(K).ctx >& $@.log

breakpoints.raw.vcf: breakpoints.txt.gz $(SEQS)
	$(MCCORTEX) calls2vcf -o $@ breakpoints.txt.gz ref.fa >& $@.log

breakpoints.sort.vcf: breakpoints.raw.vcf
	$(VCFSORT) $< > $@

breakpoints.norm.vcf.gz: breakpoints.sort.vcf
	$(BCFTOOLS) norm --site-win 5000 --fasta-ref ref.fa $< | \
	  $(BCFTOOLS) norm --rm-dup any --do-not-normalize > breakpoints.norm.vcf
	$(BGZIP) breakpoints.norm.vcf
	$(BCFTOOLS) index breakpoints.norm.vcf.gz

cmp_breakpoint: breakpoints.txt.gz ref.fa
	$(BRKCHCK) <(gzip -fcd breakpoints.txt.gz) ref.fa

# compare truth.vcf and breakpoints.norm.vcf.gz
# Check no entries private to either truth.vcf.gz or breakpoints.norm.vcf.gz
cmp_vcf: breakpoints.norm.vcf.gz truth.vcf.gz
	@$(BCFTOOLS) stats breakpoints.norm.vcf.gz truth.vcf.gz | \
	  grep '^SN\s*[01]\s' | grep -v 'number of samples' | \
		awk 'BEGIN{FS="\t"}{ if($$4 != 0){ print "Missing VCF entries!"; exit -1; } }'
	@echo 'VCF files match!'

join.k$(K).ctx: $(GRAPHS)
	$(MCCORTEX) join -o $@ $(GRAPHS)

join.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: join.k$(K).pdf

clean:
	rm -rf $(TGTS) $(SEQS)
	rm -rf ref.* breakpoints.* truth.* join.* *.log

.PHONY: all clean plots cmp_breakpoint cmp_vcf
