SHELL:=/bin/bash -euo pipefail
# Makefile
# Isaac Turner
# 2014-10-01
# Using 1Mb of chr22 and error free reads to measure effect of kmer-size

CTXDIR=../..
CTX31=$(CTXDIR)/bin/ctx31
CTXK=$(CTXDIR)/bin/ctx
# READSIM=$(CTXDIR)/libs/readsim/readsim
ALLREADS=$(CTXDIR)/libs/seq_file/scripts/perfect_covg.sh
STRCHK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_substrings.pl
CONTIG_STATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl
# DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

BWA=bwa
SAMTOOLS=$(CTXDIR)/libs/samtools/samtools

REF=$(CTXDIR)/results/data/chr22/uniq_flanks/chr22.1Mbp.uniq.fa
READLEN=100
ERR_PROFILE=$(CTXDIR)/results/data/PhiX/PhiX.1.fq.gz

KMERS=15 21 31 41 51 63 75 99

# How many contigs to pull out to find median walk distance
NSEED_WALK=100

# Count number of kmers
NKMERS=$(shell echo $(KMERS) | tr ' ' '\n' | awk 'END{print NR}')
MIN_KMER=$(shell echo $(KMERS) | awk -F ' ' '{print $$1}')

cmds=$(foreach K,$(KMERS),$(shell echo $(CTXK)$$[ ($K+31)/32*32-1 ]))

PREFIXES=$(foreach K,$(KMERS),k/perf) $(foreach K,$(KMERS),k/stoch)

# Expand these for all k values
DIRS=reads logs $(foreach K,$(KMERS),k$(K)) plots
GRAPHS=$(foreach PRE,$(PREFIXES),k$(PRE).ctx)
PATHS=$(foreach PRE,$(PREFIXES),k$(PRE).se.ctp.gz)
CONTIGS=$(foreach PRE,$(PREFIXES),k$(PRE).contigs.links.fa) \
				$(foreach PRE,$(PREFIXES),k$(PRE).contigs.plain.fa)

UNITIGS=$(CONTIGS:.fa=.rmdup.fa)
ALLCONTIGS=$(CONTIGS) $(UNITIGS)

MEDWALKS=$(foreach PRE,$(PREFIXES),k$(PRE).medwalk.plain.txt) \
         $(foreach PRE,$(PREFIXES),k$(PRE).medwalk.links.txt)

PERF_LINKS_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.links.rmdup.csv)
PERF_PLAIN_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.plain.rmdup.csv)
STOCH_LINKS_CSVS=$(foreach K,$(KMERS),k$(K)/stoch.contigs.links.rmdup.csv)
STOCH_PLAIN_CSVS=$(foreach K,$(KMERS),k$(K)/stoch.contigs.plain.rmdup.csv)

# BAMS=$(UNITIGS:.fa=.bam)

JOINCSVS=perf.links.join.csv perf.plain.join.csv \
         stoch.links.join.csv stoch.plain.join.csv

# removed $(BAMS)
TGTS=reads/perfect.fa.gz reads/error.fa.gz \
     $(GRAPHS) $(PATHS) $(ALLCONTIGS) $(JOINCSVS) \
     links.lengths.txt $(MEDWALKS)


REQ=$(DIRS) $(cmds)

all: $(TGTS) make_plots

checks: contig_stats contig_check hierachy_check

reads/perf.fa.gz: $(REF) | $(REQ)
	$(ALLREADS) $(READLEN) $(REF) | gzip -c > reads/perfect.fa.gz

reads/stoch.fa.gz: $(REF) | $(REQ)
	$(READSIM) -l $(READLEN) -r $(REF) -s reads/error

# Cortex build k=$(K)
k%.ctx: reads/perfect.fa.gz | $(REQ)
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; \
	src=`echo $* | grep -oE '(perf|stoch)' `; echo $(CTXK)$$maxk; \
	$(CTXK)$$maxk build -m 50M -k $$k --sample everybase.chr22:17M-18M.errorfree --seq reads/$$src.fa.gz $@ >& $@.log

# Thread
k%.se.ctp.gz: k%.ctx reads/perfect.fa.gz | $(REQ)
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; \
	src=`echo $* | grep -oE '(perf|stoch)' `; echo $(CTXK)$$maxk; \
	$(CTXK)$$maxk thread -m 200M --seq reads/$$src.fa.gz --out $@ $< >& $@.log

# Contigs
k%.contigs.plain.fa: k%.ctx | $(REQ)
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; echo $(CTXK)$$maxk; \
	$(CTXK)$$maxk contigs -o $@ k$*.ctx

# % is 31/perf or 15/stoch
k%.contigs.links.fa: k%.ctx k%.se.ctp.gz | $(REQ)
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; echo $(CTXK)$$maxk; \
	$(CTXK)$$maxk contigs -o $@ -p k$*.se.ctp.gz $<;

k%.medwalk.plain.txt:
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; echo $(CTXK)$$maxk; \
	walk=$$($(CTXK)$$maxk contigs --reseed --ncontigs $(NSEED_WALK) k$*.ctx 2>&1 | \
			 	  grep -ioE 'Lengths:.*median: [0-9,]*' | grep -oE '[0-9,]+$$' | tr -d ','); \
	printf "med_walk,$$walk\n" > $@

k%.medwalk.links.txt:
	k=`echo $* | grep -oE '^\d+'`; maxk=$$[ ($$k+31)/32*32-1 ]; echo $(CTXK)$$maxk; \
	walk=$$($(CTXK)$$maxk contigs --reseed --ncontigs $(NSEED_WALK) -p k$*.se.ctp.gz k$*.ctx 2>&1 | \
	  	 	  grep -ioE 'Lengths:.*median: [0-9,]*' | grep -oE '[0-9,]+$$' | tr -d ','); \
	printf "med_walk,$$walk\n" > $@

# Remove duplicates
%.rmdup.fa: %.fa
	$(CTX31) rmsubstr -k 15 -m 100M -q $< > $@

%.rmdup.bam: %.rmdup.fa
	$(BWA) mem $(REF) $< | $(SAMTOOLS) sort -O bam -T tmp - > $@

perf.links.join.csv: $(PERF_LINKS_CSVS)
perf.plain.join.csv: $(PERF_PLAIN_CSVS)
stoch.links.join.csv: $(STOCH_LINKS_CSVS)
stoch.plain.join.csv: $(STOCH_PLAIN_CSVS)
%.join.csv:
	src=`echo $* | grep -oE '(perf|stoch)' `; \
	annot=`echo $* | grep -oE '(links|plain)' `; \
	colidx=$$(echo {1,$$[{1..$(NKMERS)}*2]} | tr ' ' ','); \
	hdr=$$(printf "metric,%s\n" $$(echo $(KMERS) | sed 's/ /,k/g')); \
 	(echo $$hdr; \
 	 printf "kmer,%s\n" $$(echo $(KMERS) | tr ' ' ','); \
	 paste -d, k*/$$src.contigs.$$annot.rmdup.csv | \
	 cut -d, -f $$colidx - | tail -n +2) > $@

links.lengths.txt:
	./list-contig-lengths.sh tmp_contig_lengths links $(KMERS) > $@
	rm -rf tmp_contig_lengths

make_plots: $(JOINCSVS) | $(DIRS)
	R --vanilla -f plot-results.R --args join.links.csv join.plain.csv

# Stats
%.csv: %.fa $(MEDWALKS)
	$(CONTIG_STATS) --print-csv $< > $@
	base=$$(echo $< | grep -oE 'k[0-9]+/(perf|stoch)'); \
	annot=$$(echo $< | grep -oE 'plain|links'); \
	cat $$base.medwalk.$$annot.txt >> $@

$(DIRS):
	mkdir -p $@

# Check all contigs for k=15,... exist in k=99 etc.
hierachy_check: $(UNITIGS)
	if [[ $(NKMERS) -gt 1 ]]; then \
		echo '== ascending kmer =='
		for i in `eval echo {1..$$[$(NKMERS)-1]}`; do \
			k1=$$(echo $(KMERS) | awk '{print $$'$$i'}'); \
			k2=$$(echo $(KMERS) | awk '{print $$'$$[$$i+1]'}'); \
			F1=k$$k1/perf.contigs.links.rmdup.fa; \
			F2=k$$k2/perf.contigs.links.rmdup.fa; \
			$(STRCHK) $(MIN_KMER) 0.1 $$F1 $$F2; \
		done; \
		echo '== descending kmer =='
		for i in `eval echo {1..$$[$(NKMERS)-1]}`; do \
			k1=$$(echo $(KMERS) | awk '{print $$'$$i'}'); \
			k2=$$(echo $(KMERS) | awk '{print $$'$$[$$i+1]'}'); \
			F1=k$$k1/perf.contigs.links.rmdup.fa; \
			F2=k$$k2/perf.contigs.links.rmdup.fa; \
			$(STRCHK) $(MIN_KMER) 0.1 $$F2 $$F1; \
		done; \
	fi

contig_stats: $(ALLCONTIGS)
	for f in $(ALLCONTIGS); do echo $$f; $(CONTIG_STATS) $$f; done

# Check error rate
contig_check: $(ALLCONTIGS)
	for f in $(ALLCONTIGS); do echo $$f; $(STRCHK) $(MIN_KMER) 0.1 $$f $(REF); done

clean:
	rm -rf $(DIRS) $(TGTS)

.INTERMEDIATE: $(PERF_LINKS_CSVS) $(PLAIN_CSVS)
	.INTERMEDIATE: $(PERF_LINKS_CSVS) $(PLAIN_CSVS)
.PHONY: all clean checks contig_stats contig_check make_plots hierachy_check
