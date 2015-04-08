SHELL:=/bin/bash -euo pipefail

#
# chr1: CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA
# chr2: gCCGTAGGTAAGGGCGTTAGTGC
# chr3: GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG
#
# mx: gCCGTAGGTAAGGGCGTTAGaCGGGTTGGAGTTGGCCAAAGAAGTTCAAgg
#     .1111111111111111111.2222222222222222222222222222..
#
#     repeat of 5bp
#     gCCGTAGGTAAGGGCGTTAGTGCAGTGCAAGGCCACATTGGGACACGAGTTGAgg
#     111111111111111111111111xxxxx111111111111111111111111..
#     deletion of 1bp
#     cAACACCCTTTGGTTTAAGCCG-GTTGGAGTTGGCCAAAGAAGTTCAAgg
#     2222222222222222222222 2222222222222222222222222..

K=11
CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
BRKCHCK=$(CTXDIR)/scripts/cortex_brkpnts_check_ref.pl
VCFSORT=$(CTXDIR)/scripts/bash/vcf-sort
VCFRENAME=$(CTXDIR)/scripts/bash/vcf-rename

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
	  echo GAACACCCTTTGGTTTAAGCCGGGTTGGAGTTGGCCAAAGAAGTTCAACG; ) > $@

sample.fa:
	( echo '>sample0'; \
	  echo gCCGTAGGTAAGGGCGTTAGaCGGGTTGGAGTTGGCCAAAGAAGTTCAAgg; \
	  echo '>sample1'; \
	  echo gCCGTAGGTAAGGGCGTTAGTGCAGTGCAAGGCCACATTGGGACACGAGTTGAgg; \
	  echo '>sample2'; \
	  echo cAACACCCTTTGGTTTAAGCCGGTTGGAGTTGGCCAAAGAAGTTCAAgg; ) > $@

truth.vcf.gz:
	( printf "##fileformat=VCF4.1\n##fileDate="`date '+%Y%m%d'`"\n"; \
		printf "##reference=ref.fa\n"; \
		printf "##contig=<id=chr1,length=50>\n"; \
		printf "##contig=<id=chr2,length=23>\n"; \
		printf "##contig=<id=chr3,length=50>\n"; \
		printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"; \
		echo "chr1	18	truth0	T	TAGTGC	.	.	.	GT	1"; \
	  echo "chr3	21	truth1	CG	C	.	.	.	GT	1"; ) > truth.vcf
	$(BGZIP) truth.vcf
	$(BCFTOOLS) index truth.vcf.gz

%.k$(K).ctx: %.fa
	$(CTX) build -m 10M -k $(K) --sample $* --seq $< $@ >& $@.log

breakpoints.txt.gz: sample.k$(K).ctx ref.fa
	$(CTX) breakpoints -t 1 -m 10M --minref 5 --maxref 20 \
	                  --seq ref.fa --out $@ sample.k$(K).ctx >& $@.log

breakpoints.raw.vcf: breakpoints.txt.gz $(SEQS)
	$(CTX) calls2vcf -o $@ breakpoints.txt.gz ref.fa >& $@.log

breakpoints.sort.vcf: breakpoints.raw.vcf
	$(VCFSORT) $< > $@

breakpoints.norm.vcf.gz: breakpoints.sort.vcf
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref ref.fa --multiallelics +both $< | \
		$(VCFRENAME) > breakpoints.norm.vcf
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
	$(CTX) join -o $@ $(GRAPHS)

join.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

plots: join.k$(K).pdf

clean:
	rm -rf $(TGTS) $(SEQS)
	rm -rf ref.* breakpoints.* truth.* join.* *.log

.PHONY: all clean plots cmp_breakpoint cmp_vcf
