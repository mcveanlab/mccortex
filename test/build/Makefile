
CTX=../../bin/ctx31
SEQRND=../../libs/seq_file/bin/seqrnd
K=9

all: seq.k$(K).ctx

clean:
	rm -rf seq.fa seq.k*.ctx test.ctx A:2.fa

seq.fa:
	$(SEQRND) 50 > seq.fa

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample Wallace \
	                     --sample Gromit --seq seq.fa \
	                     --sample Trousers --seq seq.fa --seq2 seq.fa:seq.fa $@
	$(CTX) view $@

test.ctx:
	$(SEQRND) 20 > A:2.fa
	$(CTX) build -k $(K) --sample Colon -1 A:2.fa $@
	rm A:2.fa

.PHONY: all clean
