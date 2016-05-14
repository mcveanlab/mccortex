SHELL=/bin/bash -euo pipefail

# Test pop bubbles with a single sample input file

CTXDIR=../../..
CTX=$(CTXDIR)/bin/mccortex31
K=21

SEQS=seq.fa truth.fa
GRAPHS=seq.ctx popped.ctx truth.ctx

all: popped.ctx truth.ctx check

seq.fa:
	( echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGAGATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; \
		echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGAGATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; \
	  echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGAcATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; ) > $@

truth.fa:
	( echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGAGATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; \
		echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGAGATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; \
		echo CCTAGGGTGCAGTCAATTGCCAACGGTCGGGA; echo ATAACTTCTCCAAACCAGGTTCATGACAGCCAACCAA; ) > $@

%.ctx: %.fa
	$(CTX) build -k $(K) --sample $* --seq $< $@

popped.ctx: seq.ctx
	$(CTX) popbubbles --out $@ $<

check: truth.ctx popped.ctx
	diff -q <($(CTX) view -q -k popped.ctx | sort) <($(CTX) view -q -k truth.ctx | sort) && \
	echo "Kmers match."

clean:
	rm -rf $(SEQS) $(GRAPHS)

.PHONY: all clean check
