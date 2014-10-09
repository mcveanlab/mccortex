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

KMERS=15 21 31 41 51 63 75 99

# How many contigs to pull out to find median walk distance
NSEED_WALK=50

# Count number of kmers
NKMERS=$(shell echo $(KMERS) | tr ' ' '\n' | awk 'END{print NR}')
MIN_KMER=$(shell echo $(KMERS) | awk -F ' ' '{print $$1}')

cmds=$(foreach K,$(KMERS),$(shell echo $(CTXK)$$[ ($K+31)/32*32-1 ]))

# Expand these for all k values
DIRS=reads logs $(foreach K,$(KMERS),k$(K)) plots
GRAPHS=$(foreach K,$(KMERS),k$(K)/perf.ctx)
PATHS=$(foreach K,$(KMERS),k$(K)/perf.se.ctp.gz)
CONTIGS=$(foreach K,$(KMERS),k$(K)/perf.contigs.links.fa) \
				$(foreach K,$(KMERS),k$(K)/perf.contigs.plain.fa)
UNITIGS=$(CONTIGS:.fa=.rmdup.fa)
ALLCONTIGS=$(CONTIGS) $(UNITIGS)

LINKS_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.links.rmdup.csv)
PLAIN_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.plain.rmdup.csv)

BAMS=$(UNITIGS:.fa=.bam)

JOINCSVS=join.links.csv join.plain.csv

TGTS=reads/perfect.fa.gz $(GRAPHS) $(PATHS) $(ALLCONTIGS) $(JOINCSVS) \
     links.lengths.txt $(BAMS)


REQ=$(DIRS) $(cmds)

all: $(TGTS) make_plots

checks: contig_stats contig_check hierachy_check

reads/perfect.fa.gz: $(REF) | $(REQ)
	$(ALLREADS) $(READLEN) $(REF) | gzip -c > reads/perfect.fa.gz

# Cortex build k=$(K)
k%/perf.ctx: reads/perfect.fa.gz | $(REQ)
	k=$$[ ($*+31)/32*32-1 ]; echo $(CTXK)$$k; \
	$(CTXK)$$k build -m 50M -k $* --sample everybase.chr22:17M-18M.errorfree --seq reads/perfect.fa.gz $@ >& logs/perf.build.k$*.log

# Thread
k%/perf.se.ctp.gz: k%/perf.ctx reads/perfect.fa.gz | $(REQ)
	k=$$[ ($*+31)/32*32-1 ]; echo $(CTXK)$$k; \
	$(CTXK)$$k thread -m 200M --seq reads/perfect.fa.gz --out $@ k$*/perf.ctx >& logs/perf.thread.k$*.log

# Contigs
k%/perf.contigs.plain.fa: k%/perf.ctx | $(REQ)
	k=$$[ ($*+31)/32*32-1 ]; echo $(CTXK)$$k; \
	$(CTXK)$$k contigs -o $@ k$*/perf.ctx; \
	walk=$$($(CTXK)$$k contigs --reseed --ncontigs $(NSEED_WALK) k$*/perf.ctx 2>&1 | \
			 	  grep -ioE 'Lengths:.*median: [0-9,]*' | grep -oE '[0-9,]+$$' | tr -d ','); \
	printf "med_walk,$$walk\n" > k$*/perf.medwalk.plain.txt

k%/perf.contigs.links.fa: k%/perf.ctx k%/perf.se.ctp.gz | $(REQ)
	k=$$[ ($*+31)/32*32-1 ]; echo $(CTXK)$$k; \
	$(CTXK)$$k contigs -o $@ -p k$*/perf.se.ctp.gz k$*/perf.ctx; \
	walk=$$($(CTXK)$$k contigs --reseed --ncontigs $(NSEED_WALK) -p k$*/perf.se.ctp.gz k$*/perf.ctx 2>&1 | \
	  	 	  grep -ioE 'Lengths:.*median: [0-9,]*' | grep -oE '[0-9,]+$$' | tr -d ','); \
	printf "med_walk,$$walk\n" > k$*/perf.medwalk.links.txt

# Remove duplicates
%.rmdup.fa: %.fa
	$(CTX31) rmsubstr -k 15 -m 100M -q $< > $@

%.rmdup.bam: %.rmdup.fa
	$(BWA) mem $(REF) $< | $(SAMTOOLS) sort -O bam -T tmp - > $@

join.links.csv: $(LINKS_CSVS)
join.plain.csv: $(PLAIN_CSVS)
join.%.csv:
	colidx=$$(echo {1,$$[{1..$(NKMERS)}*2]} | tr ' ' ','); \
	hdr=$$(printf "metric,%s\n" $$(echo $(KMERS) | sed 's/ /,k/g')); \
 	(echo $$hdr; \
 	 printf "kmer,%s\n" $$(echo $(KMERS) | tr ' ' ','); \
	 paste -d, k*/perf.contigs.$*.rmdup.csv | \
	 cut -d, -f $$colidx - | tail -n +2) > $@

links.lengths.txt:
	./list-contig-lengths.sh tmp_contig_lengths links $(KMERS) > $@
	rm -rf tmp_contig_lengths

make_plots: $(JOINCSVS) | $(DIRS)
	R --vanilla -f plot-results.R --args join.links.csv join.plain.csv

# Stats
%.csv: %.fa
	$(CONTIG_STATS) --print-csv $< > $@
	dir=$$(echo $< | grep -oE 'k[0-9]+'); \
	annot=$$(echo $< | grep -oE 'plain|links'); \
	cat $$dir/perf.medwalk.$$annot.txt >> $@

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

.INTERMEDIATE: $(LINKS_CSVS) $(PLAIN_CSVS)
.PHONY: all clean checks contig_stats contig_check make_plots hierachy_check
