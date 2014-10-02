# Makefile
# Isaac Turner
# 2014-10-01
# Using 1Mb of chr22 and error free reads to measure effect of kmer-size

CTXDIR=../..
CTX31=$(CTXDIR)/bin/ctx31
CTXK=$(CTXDIR)/bin/ctx
# READSIM=$(CTXDIR)/libs/readsim/readsim
ALLREADS=$(CTXDIR)/libs/readsim/perfect_covg.sh
STRCHK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_substrings.pl
CONTIG_STATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl

REF=$(CTXDIR)/results/data/chr22/uniq_flanks/chr22.1Mbp.uniq.fa
READLEN=100

KMERS=15 21 31 41 51 63 75 99

cmds=$(foreach K,$(KMERS),$(shell echo $(CTXK)$$[ ($K+31)/32*32-1 ]))

# Expand these for all k values
DIRS=reads logs $(foreach K,$(KMERS),k$(K)) plots
GRAPHS=$(foreach K,$(KMERS),k$(K)/perf.ctx)
PATHS=$(foreach K,$(KMERS),k$(K)/perf.se.ctp.gz)
CONTIGS=$(foreach K,$(KMERS),k$(K)/perf.contigs.links.fa) \
        $(foreach K,$(KMERS),k$(K)/perf.contigs.links.rmdup.fa) \
        $(foreach K,$(KMERS),k$(K)/perf.contigs.plain.fa) \
        $(foreach K,$(KMERS),k$(K)/perf.contigs.plain.rmdup.fa)

CSVS=$(CONTIGS:.fa=.csv)

LINKS_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.links.rmdup.csv)
PLAIN_CSVS=$(foreach K,$(KMERS),k$(K)/perf.contigs.plain.rmdup.csv)
JOINCSVS=join.links.csv join.plain.csv

TGTS=reads/perfect.fa.gz $(GRAPHS) $(PATHS) $(CONTIGS) $(CSVS) $(JOINCSVS) \
     links.lengths.txt
     

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
	$(CTXK)$$k contigs -o $@ k$*/perf.ctx

k%/perf.contigs.links.fa: k%/perf.ctx k%/perf.se.ctp.gz | $(REQ)
	k=$$[ ($*+31)/32*32-1 ]; echo $(CTXK)$$k; \
	$(CTXK)$$k contigs -o $@ -p k$*/perf.se.ctp.gz k$*/perf.ctx

# Remove duplicates
%.rmdup.fa: %.fa
	$(CTX31) rmsubstr -k 15 -m 100M -q $< > $@

join.links.csv: $(LINKS_CSVS)
join.plain.csv: $(PLAIN_CSVS)

join.%.csv:
	colidx=$$(eval echo $$(echo $(KMERS) | awk '{print "{1,$$[{1.."NF"}*2]}"}') | tr ' ' ','); \
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

$(DIRS):
	mkdir -p $@

# Check all contigs for k=15,... exist in k=99 etc.
hierachy_check: $(CONTIGS)
	for k in $(KMERS); do $(STRCHK) $$k 0.1 k$$k/perf.contigs.links.rmdup.fa k99/perf.contigs.links.rmdup.fa; done

contig_stats: $(CONTIGS)
	for f in $(CONTIGS); do echo $$f; $(CONTIG_STATS) $$f; done

# Check error rate
contig_check: $(CONTIGS)
	for f in $(CONTIGS); do echo $$f; $(STRCHK) 15 0.1 $$f $(REF); done

clean:
	rm -rf $(DIRS) $(TGTS)

.INTERMEDIATE: $(LINKS_CSVS) $(PLAIN_CSVS)
.PHONY: all clean checks contig_stats contig_check make_plots hierachy_check
