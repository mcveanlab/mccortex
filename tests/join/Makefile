SHELL:=/bin/bash -euo pipefail

# Any kmer (K=) should work

K=7
CTXDIR=../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat

SAMPLES=$(shell echo in{,{0..2}}.k$(K).ctx)
MERGED=$(shell echo flatten013.k$(K).ctx merge.gaps.use{1..2}.k$(K).ctx)
GRAPHS=$(SAMPLES) $(MERGED) in.use2.k$(K).ctx
LOGS=$(addsuffix .log,$(GRAPHS))
TXTS=$(MERGED:.k$(K).ctx=.txt) in.txt in.use2.txt

all: $(GRAPHS) compare

seq%.fa:
	$(DNACAT) -F -n 100 > $@

in%.k$(K).ctx: seq%.fa
	$(MCCORTEX) build -m 1M -k $(K) --sample Sampe$* --seq $< $@ >& $@.log

# Output colours are {0,1,2,0+0,1+2,2}
in.k$(K).ctx: in0.k$(K).ctx in1.k$(K).ctx in2.k$(K).ctx
	$(MCCORTEX) join -o $@ 0:in0.k$(K).ctx 1:in1.k$(K).ctx 2:in2.k$(K).ctx 3:in0.k$(K).ctx 3:in0.k$(K).ctx 4:in1.k$(K).ctx 4:in2.k$(K).ctx 5:in2.k$(K).ctx >& $@.log

in.use2.k$(K).ctx: in0.k$(K).ctx in1.k$(K).ctx in2.k$(K).ctx
	$(MCCORTEX) join --ncols 2 -o $@ 0:in0.k$(K).ctx 1:in1.k$(K).ctx 2:in2.k$(K).ctx 3:in0.k$(K).ctx 3:in0.k$(K).ctx 4:in1.k$(K).ctx 4:in2.k$(K).ctx 5:in2.k$(K).ctx >& $@.log

flatten013.k$(K).ctx: in.k$(K).ctx
	$(MCCORTEX) join -o flatten013.k$(K).ctx 0:in.k$(K).ctx:1 0:in.k$(K).ctx:0 0:in.k$(K).ctx:3-3 >& $@.log

# Output colours are: {1,0,empty,empty,3} => {1,0,empty,empty,0+0}
merge.gaps.use1.k$(K).ctx: in.k$(K).ctx
	$(MCCORTEX) join --ncols 1 -o merge.gaps.use1.k$(K).ctx 1:in.k$(K).ctx:0 0:in.k$(K).ctx:1 4:in.k$(K).ctx:3 >& $@.log
merge.gaps.use2.k$(K).ctx: in.k$(K).ctx
	$(MCCORTEX) join --ncols 2 -o merge.gaps.use2.k$(K).ctx 1:in.k$(K).ctx:0 0:in.k$(K).ctx:1 4:in.k$(K).ctx:3 >& $@.log

%.txt: %.k$(K).ctx
	$(MCCORTEX) view -q --kmers $< | sort > $@

compare: $(TXTS)
	diff -q in.txt in.use2.txt
	diff -q merge.gaps.use*.txt

clean:
	rm -rf $(GRAPHS) $(TXTS) seq*.fa $(LOGS)

.PHONY: all clean compare
