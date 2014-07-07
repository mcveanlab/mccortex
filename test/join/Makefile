SHELL := /bin/bash

CTX=../../bin/ctx31
SEQRND=../../libs/seq_file/bin/seqrnd

TGTS=$(shell echo in.ctx in{0..2}.ctx flatten013.use{1..2}.ctx merge013.use{1..2}.ctx merge.gaps.use{1..2}.ctx)
TXT=in.txt flatten013.txt merge013.txt merge.gaps.txt

all: $(TGTS) $(TXT)
	diff -q flatten013.use*.ctx
	diff -q merge.gaps.use*.ctx
	diff -q merge013.use*.ctx
	head *.txt

seq%.fa:
	$(SEQRND) 100 > $@

in%.ctx: seq%.fa
	$(CTX) build -k 7 --sample Sampe$* --seq $< $@

# Output colours are {0,1,2,0+0,1+2,2}
in.ctx: in0.ctx in1.ctx in2.ctx
	$(CTX) join --overlap -o in.ctx 0:in0.ctx 1:in1.ctx 2:in2.ctx 3:in0.ctx 3:in0.ctx 4:in1.ctx 4:in2.ctx 5:in2.ctx

# flatten013.ctx and merge013.ctx should be identical (int.ctx:0,1,3 flattened)
flatten013.use1.ctx: in.ctx
	$(CTX) join --flatten -o flatten013.use1.ctx in.ctx:1 10:in.ctx:0 in.ctx:3-3
flatten013.use2.ctx: in.ctx
	$(CTX) join --flatten --ncols 2 -o flatten013.use2.ctx in.ctx:1 10:in.ctx:0 in.ctx:3-3

merge013.use1.ctx: in.ctx
	$(CTX) join --overlap -o merge013.use1.ctx in.ctx:0 in.ctx:1 in.ctx:3
merge013.use2.ctx: in.ctx
	$(CTX) join --overlap --ncols 2 -o merge013.use2.ctx in.ctx:0 in.ctx:1 in.ctx:3

# Output colours are: {1,0,empty,empty,3} => {1,0,empty,empty,0+0}
merge.gaps.use1.ctx: in.ctx
	$(CTX) join --overlap -o merge.gaps.use1.ctx 1:in.ctx:0 0:in.ctx:1 4:in.ctx:3
merge.gaps.use2.ctx: in.ctx
	$(CTX) join --overlap --ncols 2 -o merge.gaps.use2.ctx 1:in.ctx:0 0:in.ctx:1 4:in.ctx:3

in.txt: in.ctx
	$(CTX) view --kmers in.ctx | sort > in.txt
flatten013.txt: flatten013.use1.ctx
	$(CTX) view --kmers flatten013.use1.ctx | sort > flatten013.txt
merge013.txt: merge013.use1.ctx
	$(CTX) view --kmers merge013.use1.ctx | sort > merge013.txt
merge.gaps.txt: merge.gaps.use1.ctx
	$(CTX) view --kmers merge.gaps.use1.ctx | sort > merge.gaps.txt

clean:
	rm -rf $(TGTS) $(TXT) seq*.fa

.PHONY: all clean
