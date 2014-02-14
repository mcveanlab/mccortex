
SEQRND=../../libs/seq_file/bin/seqrnd
FACAT=../../libs/seq_file/bin/facat
CTX=../../bin/ctx31
BIOINF=../../libs/bioinf-perl

all: test

seq.fa:
	$(SEQRND) 1001 | $(FACAT) > seq.fa

seq.k9.ctx: seq.fa
	$(CTX) build -k 9 --sample MrsSample --seq seq.fa seq.k9.ctx

test: seq.fa seq.k9.ctx
	 $(CTX) contigs --ncontigs 10 --print seq.k9.ctx | $(BIOINF)/sim_mutations/sim_substrings.pl 9 0.1 - seq.fa

clean:
	rm -rf seq.fa seq.k9.ctx

.PHONY: all clean test
