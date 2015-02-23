SHELL=/bin/bash -euo pipefail

#
# 1 het site in s0, 1 het between s0,s1
#
#                                                        v     v ref mismatch
# ref:GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACcTTGGGcCACGAGTTGATATCATCCGGCGTAGGAGGACT
#                                     x                  y
# s0: GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT
# s0: GTAGACTACTTGTCGCGACCCCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT
#
# s1: GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT
#

K=11
CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
CTX2DOT=$(CTXDIR)/scripts/cortex_to_graphviz.pl
CTXFLANKS=$(CTXDIR)/scripts/cortex_print_flanks.sh
VCFSORT=$(CTXDIR)/scripts/bash/vcf-sort
VCFRENAME=$(CTXDIR)/scripts/bash/vcf-rename

BWA=$(CTXDIR)/libs/bwa/bwa
BGZIP=$(CTXDIR)/libs/htslib/bgzip
BCFTOOLS=$(CTXDIR)/libs/bcftools/bcftools

SEQS=seq0.fa seq1.fa
GRAPHS=$(SEQS:.fa=.k$(K).ctx)
FLAT_GRAPHS=$(shell printf "0:%s\n" $(GRAPHS))
BUBBLES=bubbles.seq.txt.gz bubbles.join.txt.gz
PLOTS=seq.k$(K).pdf
TGTS=$(SEQS) $(GRAPHS) join.k$(K).ctx $(BUBBLES)

all: $(TGTS) cmp_bubbles cmp_vcf

test:
	echo $(SEQS)
	echo $(GRAPHS)

ref.fa:
	(echo '>chr1'; echo GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACcTTGGGcCACGAGTTGATATCATCCGGCGTAGGAGGACT;) > $@

seq0.fa:
	echo GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT > $@
	echo GTAGACTACTTGTCGCGACCCCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT >> $@

seq1.fa:
	echo GTAGACTACTTGTCGCGACCCCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGATATCATCCGGCGTAGGAGGACT > $@

truth.vcf.gz:
	( printf "##fileformat=VCF4.1\n##fileDate="`date '+%Y%m%d'`"\n"; \
		printf "##reference=ref.fa\n"; \
		printf "##contig=<id=chr1,length=90>\n"; \
		printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"; \
		echo "chr1	33	truth0	G	T	.	.	.	."; \
	  echo "chr1	52	truth1	C	A,T	.	.	.	."; ) > truth.vcf
	$(BGZIP) truth.vcf
	$(BCFTOOLS) index truth.vcf.gz

seq%.k$(K).ctx: seq%.fa
	$(CTX) build -m 10M -k $(K) --sample $< --seq $< $@ >& $@.log

join.k$(K).ctx: $(GRAPHS)
	$(CTX) join -o $@ $(FLAT_GRAPHS) >& $@.log

bubbles.seq.txt.gz: $(GRAPHS)
	$(CTX) bubbles -t 1 -m 10M -o $@ $(GRAPHS) >& $@.log

bubbles.join.txt.gz: join.k$(K).ctx
	$(CTX) bubbles -t 1 -m 10M -o $@ $(GRAPHS) >& $@.log

cmp_bubbles: bubbles.join.txt.gz bubbles.seq.txt.gz
	@L1=`gzip -cd bubbles.join.txt.gz | wc -l | grep -oE '[0-9]+'`; \
	L2=`gzip -cd bubbles.seq.txt.gz  | wc -l | grep -oE '[0-9]+'`; \
	echo bubbles.join.txt.gz [$$L1], bubbles.seq.txt.gz [$$L2] have same number of lines

seq.k$(K).pdf: join.k$(K).ctx
	$(CTX2DOT) --simplify $< | dot -Tpdf > $@

bubbles.raw.vcf: ref.fa bubbles.join.txt.gz
	$(BWA) index ref.fa
	$(CTXFLANKS) bubbles.join.txt.gz > bubbles.join.txt.flanks
	$(BWA) mem ref.fa bubbles.join.txt.flanks > bubbles.join.txt.sam
	$(CTX) calls2vcf -F bubbles.join.txt.sam -o $@ bubbles.join.txt.gz ref.fa >& $@.log

bubbles.sort.vcf: bubbles.raw.vcf
	$(VCFSORT) $< > $@

bubbles.norm.vcf.gz: bubbles.sort.vcf
	$(BCFTOOLS) norm --remove-duplicates --fasta-ref ref.fa --multiallelics +both $< | \
	$(VCFRENAME) > bubbles.norm.vcf
	$(BGZIP) -f bubbles.norm.vcf
	$(BCFTOOLS) index bubbles.norm.vcf.gz

# compare truth.vcf and bubbles.norm.vcf.gz
# Check no entries private to either truth.vcf.gz or bubbles.norm.vcf.gz
cmp_vcf: bubbles.norm.vcf.gz truth.vcf.gz
	$(BCFTOOLS) stats bubbles.norm.vcf.gz truth.vcf.gz | grep -E '^SN\s*[01]\s' | \
		awk 'BEGIN{FS="\t"}{ if($$4 != 0){ print "Missing VCF entries!"; exit -1; } }'
	@echo 'VCF files match!'

plots: $(PLOTS)

clean:
	rm -rf $(TGTS) $(PLOTS) *.log
	rm -rf ref.fa* bubbles.* truth.vcf*

.PHONY: all clean cmp_bubbles cmp_vcf plots
