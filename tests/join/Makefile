SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
CTX=$(CTXDIR)/bin/mccortex31

SAMPLES=$(shell echo in{,{0..2}}.ctx)
MERGED=$(shell echo flatten013.ctx merge.gaps.use{1..2}.ctx)
GRAPHS=$(SAMPLES) $(MERGED) in.use2.ctx
TXTS=$(MERGED:.ctx=.txt) in.txt in.use2.txt

all: $(GRAPHS) compare

seq%.fa:
	$(DNACAT) -F -n 100 > $@

in%.ctx: seq%.fa
	$(CTX) build -m 1M -k 7 --sample Sampe$* --seq $< $@

# Output colours are {0,1,2,0+0,1+2,2}
in.ctx: in0.ctx in1.ctx in2.ctx
	$(CTX) join -o $@ 0:in0.ctx 1:in1.ctx 2:in2.ctx 3:in0.ctx 3:in0.ctx 4:in1.ctx 4:in2.ctx 5:in2.ctx

in.use2.ctx: in0.ctx in1.ctx in2.ctx
	$(CTX) join --ncols 2 -o $@ 0:in0.ctx 1:in1.ctx 2:in2.ctx 3:in0.ctx 3:in0.ctx 4:in1.ctx 4:in2.ctx 5:in2.ctx

flatten013.ctx: in.ctx
	$(CTX) join -o flatten013.ctx 0:in.ctx:1 0:in.ctx:0 0:in.ctx:3-3

# Output colours are: {1,0,empty,empty,3} => {1,0,empty,empty,0+0}
merge.gaps.use1.ctx: in.ctx
	$(CTX) join --ncols 1 -o merge.gaps.use1.ctx 1:in.ctx:0 0:in.ctx:1 4:in.ctx:3
merge.gaps.use2.ctx: in.ctx
	$(CTX) join --ncols 2 -o merge.gaps.use2.ctx 1:in.ctx:0 0:in.ctx:1 4:in.ctx:3

%.txt: %.ctx
	$(CTX) view --kmers $< | sort > $@

compare: $(TXTS)
	diff -q in.txt in.use2.txt
	diff -q merge.gaps.use*.txt

clean:
	rm -rf $(GRAPHS) $(TXTS) seq*.fa

.PHONY: all clean compare
