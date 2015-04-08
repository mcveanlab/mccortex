SHELL:=/bin/bash -euo pipefail

CTXDIR=../..
CTX=$(CTXDIR)/bin/mccortex31

LAST=5
INPUT_FILES=$(shell echo input.{0..$(LAST)}.fa)
OUTPUT_FILES=$(shell echo output.{0..$(LAST)}.fa)
RESULT_FILES=$(shell echo results.{0..$(LAST)}.fa)
TESTS=$(shell echo test.{0..$(LAST)})

all: $(INPUT_FILES) $(OUTPUT_FILES) $(RESULT_FILES) test

input.0.fa:
	printf '>a\nAAA\n'   > $@
	printf '>b\nAAA\n'  >> $@

results.0.fa:
	printf '>a\nAAA\n'   > $@

input.1.fa:
	printf '>a\naaat\n'  > $@
	printf '>b\nAAT\n'  >> $@

results.1.fa:
	printf '>a\naaat\n'  > $@

input.2.fa:
	printf '>a\nAAT\n'   > $@
	printf '>b\naaat\n' >> $@

results.2.fa:
	printf '>b\naaat\n' > $@

input.3.fa:
	printf '>a\nTAATAATGACCCTGTGTGACGGAA\n'        > $@
	printf '>b\nTAATAATGACCCTGTGTGACGGAAGGAAAG\n' >> $@
	printf '>c\nTAATAATGACCCTGTGTGACGGAAGGAAAG\n' >> $@
	printf '>d\nATATCGAATAGTGA\n'                 >> $@
	printf '>e\nCCCTGTGTGACGGAAGGAAAG\n'          >> $@
	printf '>f\nTAATAATGACCCTGTGTGAC\n'           >> $@
	printf '>g\nTAATAATGACCCTGTGTGACGGAAGGAAAG\n' >> $@
	printf '>h\nCGAATAGTGAGGAGTT\n'               >> $@
	printf '>i\nATATCGAATAGTGAGGAGTT\n'           >> $@

results.3.fa:
	printf '>b\nTAATAATGACCCTGTGTGACGGAAGGAAAG\n' > $@
	printf '>i\nATATCGAATAGTGAGGAGTT\n'          >> $@

input.4.fa:
	printf '>a\nNAACAT\n' > $@
	printf '>b\nAAC\n'   >> $@
	printf '>c\nnAAC\n'  >> $@
	printf '>d\nbAAC\n'  >> $@
	printf '>e\nNACAT\n' >> $@

results.4.fa:
	printf '>a\nNAACAT\n' > $@
	printf '>d\nbAAC\n'  >> $@
	printf '>e\nNACAT\n' >> $@

input.5.fa:
	printf '>a\nACACGATCTT\n'            > $@
#     read a is  <--------|
	printf '>b\nTACAAGATCGTGTGTTGCCT\n' >> $@
	printf '>c\nTAACAGCAACCCTGACC\n'    >> $@
	printf '>d\nACACGATCTTGGATATACAG\n' >> $@

results.5.fa:
	printf '>b\nTACAAGATCGTGTGTTGCCT\n'  > $@
	printf '>c\nTAACAGCAACCCTGACC\n'    >> $@
	printf '>d\nACACGATCTTGGATATACAG\n' >> $@

output.%.fa: input.%.fa
	$(CTX) rmsubstr -q -n 1024 -k 3 $< > $@

test.%: output.%.fa results.%.fa
	diff -q output.$*.fa results.$*.fa

test: $(TESTS)

clean:
	rm -rf $(INPUT_FILES) $(OUTPUT_FILES) $(RESULT_FILES)

.PHONY: all clean test
