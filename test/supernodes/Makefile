
K=7
GENOME=200
KEEP=$(shell echo genome.fa genome.k$(K).{ctx,dot,pdf} genome.k$(K).old.pdf)

all: $(KEEP)

clean:
	rm -rf $(KEEP)

# Sample random genome
genome.fa:
	seqrnd $(GENOME) > genome.fa

genome.k$(K).ctx:
	ctx31 build -k $(K) --sample MssrGenome --seq genome.fa genome.k$(K).ctx

genome.k$(K).dot:
	ctx31 supernodes --graphviz --points genome.k$(K).ctx > genome.k$(K).dot

genome.k$(K).pdf:
	dot -Tpdf genome.k$(K).dot > genome.k$(K).pdf

genome.k$(K).old.pdf:
	../../scripts/cortex_to_graphviz.pl --simplify genome.k$(K).ctx | dot -Tpdf > genome.k$(K).old.pdf

.PHONY: all clean
