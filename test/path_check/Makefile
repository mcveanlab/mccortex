
K=7
GENOME=200
KEEP=$(shell echo genome.fa genome.k$(K).{ctx,ctp,dot,pdf} reads.1.fa.gz reads.2.fa.gz reads.{se,pe}.k$(K).ctp)

all: $(KEEP)

clean:
	rm -rf $(KEEP)
	rm -rf gap_sizes.*.csv

# Sample random genome
genome.fa:
	seqrnd $(GENOME) > genome.fa

genome.k$(K).ctx:
	ctx31 build -m 10M -k $(K) --sample MssrGenome --seq genome.fa genome.k$(K).ctx

genome.k$(K).dot:
	ctx31 supernodes --graphviz genome.k$(K).ctx > genome.k$(K).dot

genome.k$(K).pdf:
	dot -Tpdf genome.k$(K).dot > genome.k$(K).pdf

reads.1.fa.gz reads.2.fa.gz:
	../../libs/readsim/readsim -r genome.fa -l 10 -i 20 -v 0.1 -d 2 reads

genome.k$(K).ctp:
	ctx31 thread -m 10M -a 1 -t 1 --col 0 --seq genome.fa genome.k$(K).ctp genome.k$(K).ctx

reads.se.k$(K).ctp:
	ctx31 thread -m 10M -a 1 -t 1 --col 0 --seq reads.1.fa.gz --seq reads.2.fa.gz reads.se.k$(K).ctp genome.k$(K).ctx

reads.pe.k$(K).ctp:
	ctx31 thread -m 10M -a 1 -t 1 --col 0 --seq2 reads.1.fa.gz reads.2.fa.gz reads.pe.k$(K).ctp genome.k$(K).ctx

.PHONY: all clean
