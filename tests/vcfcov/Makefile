SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

# call0: test with empty VCF
# call1: overlapping SNPs at ref:1,50,199, chr1:50. Handle missing contig= header.
# call2: SNP every base, across two chromosomes
# call3: blocks of overlapping variants (on one chrom)
# call4: variants exactly k-1 bases apart (on one chrom)
# call5: test for large indels

all:
	cd calls0 && $(MAKE) clean && $(MAKE)
	cd calls1 && $(MAKE) clean && $(MAKE)
	cd calls2 && $(MAKE) clean && $(MAKE)
	cd calls3 && $(MAKE) clean && $(MAKE)
	cd calls4 && $(MAKE) clean && $(MAKE)
	cd calls5 && $(MAKE) clean && $(MAKE)
	@echo "All looks good."

clean:
	cd calls0 && $(MAKE) clean
	cd calls1 && $(MAKE) clean
	cd calls2 && $(MAKE) clean
	cd calls3 && $(MAKE) clean
	cd calls4 && $(MAKE) clean
	cd calls5 && $(MAKE) clean

.PHONY: all clean view
