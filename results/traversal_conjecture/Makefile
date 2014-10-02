# Makefile
# Isaac Turner
# 2014-09-05
# Using 1Mb of chr22 and error free reads with both perfect & stochstic coverage

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx$(K)
READSIM=$(CTXDIR)/libs/readsim/readsim
ALLREADS=$(CTXDIR)/libs/readsim/perfect_covg.sh
STRCHK=$(CTXDIR)/libs/bioinf-perl/sim_mutations/sim_substrings.pl
CONTIG_STATS=$(CTXDIR)/libs/bioinf-perl/fastn_scripts/contig_stats.pl

REF=chr22.1Mbp.uniq.fa
READLEN=100 # if we start a 100bp read at each base, coverage is 100X
DEPTH=100
K=31
NRUNS=10000
ABDIST=500

DIRS=reads logs k$(K)

CONTIGS=k$(K)/perf.contigs.fa k$(K)/perf.contigs.rmdup.fa \
        k$(K)/stoch.contigs.fa k$(K)/stoch.contigs.rmdup.fa

TGTS=reads/stochstic.1.fa.gz reads/stochstic.2.fa.gz reads/perfect.fa.gz \
     k$(K)/stoch.k$(K).ctx k$(K)/stoch.k$(K).se.ctp.gz \
     k$(K)/perf.k$(K).ctx  k$(K)/perf.k$(K).se.ctp.gz \
     $(CONTIGS)

all: $(TGTS) test_abc contig_stats contig_check

# Simulate PE reads of each chrom each 50X, only use
reads/stochstic.1.fa.gz: reads/stochstic.2.fa.gz
reads/stochstic.2.fa.gz: $(REF) | $(DIRS)
	$(READSIM) -l $(READLEN) -r $(REF) -d $(DEPTH) reads/stochstic >& logs/stoch.reads.log

reads/perfect.fa.gz: $(REF) | $(DIRS)
	$(ALLREADS) $(READLEN) $(REF) | gzip -c > reads/perfect.fa.gz

# Cortex build k=$(K)
k$(K)/stoch.k$(K).ctx: reads/stochstic.1.fa.gz reads/stochstic.2.fa.gz | $(DIRS)
	$(CTX) build -m 50M -k $(K) --sample $(DEPTH)X.chr22:17M-18M.errorfree --seq reads/stochstic.1.fa.gz --seq reads/stochstic.2.fa.gz $@ >& logs/stoch.build.k$(K).log

k$(K)/perf.k$(K).ctx: reads/perfect.fa.gz | $(DIRS)
	$(CTX) build -m 50M -k $(K) --sample everybase.chr22:17M-18M.errorfree --seq reads/perfect.fa.gz $@ >& logs/perf.build.k$(K).log

# Thread
k$(K)/stoch.k$(K).se.ctp.gz: k$(K)/stoch.k$(K).ctx reads/stochstic.1.fa.gz reads/stochstic.2.fa.gz | $(DIRS)
	$(CTX) thread -m 50M --seq reads/stochstic.1.fa.gz --seq reads/stochstic.2.fa.gz --out $@ k$(K)/stoch.k$(K).ctx >& logs/stoch.thread.k$(K).log

k$(K)/perf.k$(K).se.ctp.gz: k$(K)/perf.k$(K).ctx reads/perfect.fa.gz | $(DIRS)
	$(CTX) thread -m 50M --seq reads/perfect.fa.gz --out $@ k$(K)/perf.k$(K).ctx >& logs/perf.thread.k$(K).log

# Contigs
k$(K)/%.contigs.fa: k$(K)/%.k$(K).ctx k$(K)/%.k$(K).se.ctp.gz | $(DIRS)
	$(CTX) contigs -o $@ -p k$(K)/$*.k$(K).se.ctp.gz k$(K)/$*.k$(K).ctx

# Remove duplicates
k$(K)/%.contigs.rmdup.fa: k$(K)/%.contigs.fa
	$(CTX) rmsubstr -m 50M -q $< > $@

test_abc: $(TGTS)
	$(CTX) exp_abc -N $(NRUNS) --max-AB-dist $(ABDIST) -p k$(K)/perf.k$(K).se.ctp.gz k$(K)/perf.k$(K).ctx
	$(CTX) exp_abc -N $(NRUNS) --max-AB-dist $(ABDIST) -p k$(K)/stoch.k$(K).se.ctp.gz k$(K)/stoch.k$(K).ctx

$(DIRS):
	mkdir -p $@

contig_stats: $(CONTIGS)
	$(CONTIG_STATS) k$(K)/perf.contigs.fa
	$(CONTIG_STATS) k$(K)/perf.contigs.rmdup.fa
	$(CONTIG_STATS) k$(K)/stoch.contigs.fa
	$(CONTIG_STATS) k$(K)/stoch.contigs.rmdup.fa

# Check error rate
contig_check: $(CONTIGS)
	$(STRCHK) $(K) 0.1 k$(K)/perf.contigs.fa $(REF)
	$(STRCHK) $(K) 0.1 k$(K)/perf.contigs.rmdup.fa $(REF)
	$(STRCHK) $(K) 0.1 k$(K)/stoch.contigs.fa $(REF)
	$(STRCHK) $(K) 0.1 k$(K)/stoch.contigs.rmdup.fa $(REF)

clean:
	rm -rf $(DIRS)

.PHONY: all clean contig_stats contig_check test_abc
