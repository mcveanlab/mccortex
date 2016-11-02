SHELL:=/bin/bash -euo pipefail

#
# Clean a multi-sample graph by loading all samples at once. Each sample has
# only one copy, remove <2 copies from population. 
# When only 2 samples are loaded we dump an empty graph.
# When all 3 samples are loaded we dump a full graph.
#
# works with all kmers
#

K=31
CTXDIR=../../..
MCCORTEX=$(CTXDIR)/bin/mccortex $(K)

SEQ_FILES=sample.1.fa sample.2.fa sample.3.fa
RAW_GRAPHS=pop2.k$(K).ctx pop3.k$(K).ctx
CLEAN_GRAPHS=clean2.k$(K).ctx clean2b.k$(K).ctx clean3.k$(K).ctx
TGTS= $(SEQ_FILES) $(RAW_GRAPHS) $(CLEAN_GRAPHS)

all: $(TGTS) check

sample.1.fa:
	echo A > $@

sample.2.fa:
	echo GCTTCTTATTTGGCATAATCCAACTTCCCTACGGAAGCCCAATAGGATTAAATTGAAGCT > $@

sample.3.fa:
	echo GCTTCTTATTTGGCATAATCCAACTTCCCTACGGAAGCCCAATAGGATTAAATTGAAGCT > $@

pop2.k$(K).ctx: $(SEQ_FILES)
	$(MCCORTEX) build -q -k $(K) --sample sample1 --seq sample.1.fa \
	                             --sample sample2 --seq sample.2.fa $@

pop3.k$(K).ctx: $(SEQ_FILES)
	$(MCCORTEX) build -q -k $(K) --sample sample1 --seq sample.1.fa \
	                             --sample sample2 --seq sample.2.fa \
	                             --sample sample2 --seq sample.3.fa $@

clean%.k$(K).ctx: pop%.k$(K).ctx
	$(MCCORTEX) clean -q -m 1M --unitigs=2 --tips=0 -o $@ $<

clean2b.k$(K).ctx: pop3.k$(K).ctx
	$(MCCORTEX) clean -q -m 1M --unitigs=2 --tips=0 -o $@ pop3.k$(K).ctx:0,1

check: $(RAW_GRAPHS) $(CLEAN_GRAPHS)
	diff -q <($(MCCORTEX) view -qk clean2.k$(K).ctx | cut -d' ' -f1 | sort) <(echo -n '')
	diff -q <($(MCCORTEX) view -qk clean2b.k$(K).ctx | cut -d' ' -f1 | sort) <(echo -n '')
	diff -q <($(MCCORTEX) view -qk clean3.k$(K).ctx | cut -d' ' -f1 | sort) \
	        <($(MCCORTEX) view -qk pop3.k$(K).ctx | cut -d' ' -f1 | sort)

clean:
	rm -rf $(TGTS)

.PHONY: all clean check
