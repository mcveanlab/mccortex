SHELL:=/bin/bash

CTXDIR=../..
CTX=$(CTXDIR)/bin/ctx31
DNACAT=$(CTXDIR)/libs/seq_file/bin/dnacat
K=11

# 1) Stream to output file
# 2) Edit in place
# 3) File to output file

# 1) generate sequence, construct graph
# 2) generate graph without edges

# 1) pop edges
# 2) all edges

GRAPHS=seq.k$(K).ctx noedges.k$(K).ctx \
       mix.00.k$(K).ctx mix.01.k$(K).ctx mix.10.k$(K).ctx mix.11.k$(K).ctx \

FIXED=$(shell echo fix.1{a..f}.ctx fix.{2,3,4}{a,b}.ctx)

# inferedges --pop test using CAAGG kmer
CAAGG_TXTS=$(shell echo CAAGG.infer.k5.col{0..4}.txt CAAGG.k5.col0.txt CAAGG.{left,right}.edges.k5.txt empty.k5.txt)

TXTS=$(GRAPHS:.ctx=.txt) $(FIXED:.ctx=.txt) $(CAAGG_TXTS)

all: $(TXTS) $(FIXED) $(GRAPHS) test

seq.fa:
	$(DNACAT) -F -n 100 > $@

seq.k$(K).ctx: seq.fa
	$(CTX) build -k $(K) --sample PinkPanther --seq $< $@
	$(CTX) inferedges --all $@

noedges.k$(K).ctx: seq.k$(K).ctx
	set -euo pipefail; \
	$(CTX) build -k $(K) -m 1M --sample Puma \
	             --seq <($(CTX) view -q -k $< | awk '{print $$1}') $@

mix.11.k$(K).ctx: seq.k$(K).ctx noedges.k$(K).ctx
	$(CTX) join -m 1M -o $@ seq.k$(K).ctx seq.k$(K).ctx

mix.10.k$(K).ctx: seq.k$(K).ctx noedges.k$(K).ctx
	$(CTX) join -m 1M -o $@ seq.k$(K).ctx noedges.k$(K).ctx

mix.01.k$(K).ctx: seq.k$(K).ctx noedges.k$(K).ctx
	$(CTX) join -m 1M -o $@ noedges.k$(K).ctx seq.k$(K).ctx

mix.00.k$(K).ctx: seq.k$(K).ctx noedges.k$(K).ctx
	$(CTX) join -m 1M -o $@ noedges.k$(K).ctx noedges.k$(K).ctx

# 1) Stream
fix.1a.ctx: mix.10.k$(K).ctx
	set -euo pipefail; cat $< | $(CTX) inferedges -m 1M --all - > $@
fix.1b.ctx: mix.10.k$(K).ctx
	set -euo pipefail; cat $< | $(CTX) inferedges -m 1M --pop - > $@

# Edit file in place (mmap)
fix.1c.ctx: mix.10.k$(K).ctx
	cp mix.10.k$(K).ctx $@
	$(CTX) inferedges --all $@
fix.1d.ctx: mix.10.k$(K).ctx
	cp mix.10.k$(K).ctx $@
	$(CTX) inferedges --pop $@

# Separate input + output files
fix.1e.ctx: mix.10.k$(K).ctx
	$(CTX) inferedges --all -o $@ $<
fix.1f.ctx: mix.10.k$(K).ctx
	$(CTX) inferedges --pop -o $@ $<

# 2) fix on 01 is same as 10
fix.2a.ctx: mix.01.k$(K).ctx
	$(CTX) inferedges --all -o $@ $<
fix.2b.ctx: mix.01.k$(K).ctx
	$(CTX) inferedges --pop -o $@ $<

# 3) fix on 11 is just 11
fix.3a.ctx: mix.11.k$(K).ctx
	$(CTX) inferedges --all -o $@ $<
fix.3b.ctx: mix.11.k$(K).ctx
	$(CTX) inferedges --pop -o $@ $<

# fix on 00 is just 00
fix.4a.ctx: mix.00.k$(K).ctx
	$(CTX) inferedges --all -o $@ $<
fix.4b.ctx: mix.00.k$(K).ctx
	$(CTX) inferedges --pop -o $@ $<

#
#  CAAGG TEST: Tests inferring population kmers
#
# ACAAG             AAGGA
# CCAAG*-> CAAGG -> AAGGC
# GCAAG             AAGGG
# TCAAG            *AAGGT  * = edges
#
CAAGG.prev.txt:
	printf "ACAAG\nCCAAG\nGCAAG\nTCAAG\n" > $@

CAAGG.next.txt:
	printf "AAGGA\nAAGGC\nAAGGG\nAAGGT\n" > $@

CAAGG.left.kmers.k5.ctx: CAAGG.prev.txt
	$(CTX) build -m 1M -k 5 --sample LeftKmers --seq <(echo CAAGG; cat CAAGG.prev.txt) $@

CAAGG.right.kmers.k5.ctx: CAAGG.next.txt
	$(CTX) build -m 1M -k 5 --sample RightKmers --seq <(echo CAAGG; cat CAAGG.next.txt) $@

CAAGG.left.edges.k5.ctx:
	$(CTX) build -m 1M -k 5 --sample LeftEdges --seq <(printf "ACAAG\nCCAAGG\nGCAAG\nTCAAG") $@

CAAGG.right.edges.k5.ctx:
	$(CTX) build -m 1M -k 5 --sample RightEdges --seq <(printf "AAGGA\nAAGGC\nAAGGG\nCAAGGT\n") $@

empty.k5.ctx:
	$(CTX) build -m 1M -k 5 --sample MsEmpty --seq <(echo) $@

# Colour 0 has all kmers and 2 edges
# Colour 1 has left kmers but no edges (gets 1 from infer)
# Colour 2 has right kmers but no edges (gets 1 from infer)
# Colour 3 has all kmers but no edges (gets them from infer)
# Colour 4 has no kmers
CAAGG.k5.ctx: CAAGG.left.kmers.k5.ctx CAAGG.right.kmers.k5.ctx CAAGG.left.edges.k5.ctx CAAGG.right.edges.k5.ctx empty.k5.ctx
	$(CTX) join --out $@ \
	  0:CAAGG.left.edges.k5.ctx 0:CAAGG.right.edges.k5.ctx \
    1:CAAGG.left.kmers.k5.ctx \
    2:CAAGG.right.kmers.k5.ctx \
    3:CAAGG.left.kmers.k5.ctx 3:CAAGG.right.kmers.k5.ctx \
    4:empty.k5.ctx \

CAAGG.infer.k5.ctx: CAAGG.k5.ctx
	$(CTX) inferedges --pop -o $@ $<

CAAGG.infer.k5.col%.txt: CAAGG.infer.k5.ctx
	set -euo pipefail; $(CTX) view --kmers $<:$* | sort > $@

CAAGG.k5.col%.txt: CAAGG.k5.ctx
	set -euo pipefail; $(CTX) view --kmers $<:$* | sort > $@

%.txt: %.ctx
	set -euo pipefail; $(CTX) view --kmers $< | sort > $@

test: $(TXTS) $(CAAGG_TXTS)
	# diff -q fix.1a.txt mix.11.k$(K).txt
	diff -q fix.1b.txt mix.11.k$(K).txt
	diff -q fix.1c.txt fix.1a.txt
	diff -q fix.1d.txt mix.11.k$(K).txt
	diff -q fix.1e.txt fix.1a.txt
	diff -q fix.1f.txt mix.11.k$(K).txt

	diff -q fix.2a.txt fix.1a.txt
	diff -q fix.2b.txt mix.11.k$(K).txt

	diff -q fix.3a.txt fix.1a.txt
	diff -q fix.3b.txt mix.11.k$(K).txt

	diff -q fix.4a.txt fix.1a.txt
	diff -q fix.4b.txt mix.00.k$(K).txt

	diff -q CAAGG.infer.k5.col0.txt CAAGG.k5.col0.txt
	diff -q CAAGG.infer.k5.col1.txt CAAGG.left.edges.k5.txt
	diff -q CAAGG.infer.k5.col2.txt CAAGG.right.edges.k5.txt
	diff -q CAAGG.infer.k5.col3.txt CAAGG.k5.col0.txt
	diff -q CAAGG.infer.k5.col4.txt empty.k5.txt

clean:
	rm -rf seq.fa *.txt *.ctx

.PHONY: all clean
