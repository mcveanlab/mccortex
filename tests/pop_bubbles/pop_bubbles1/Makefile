SHELL=/bin/bash -euo pipefail

# Test pop bubbles with a single sample input file

K=21
CTXDIR=../../..
MCCORTEX=$(shell echo $(CTXDIR)/bin/mccortex$$[(($(K)+31)/32)*32 - 1])

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
	$(MCCORTEX) build -q -k $(K) --sample $* --seq $< $@

popped.ctx: seq.ctx
	$(MCCORTEX) popbubbles -q --out $@ $<

check: truth.ctx popped.ctx
	diff -q <($(MCCORTEX) view -qk popped.ctx | sort) <($(MCCORTEX) view -qk truth.ctx | sort) && \
	echo "Kmers match."

clean:
	rm -rf $(SEQS) $(GRAPHS)

.PHONY: all clean check
