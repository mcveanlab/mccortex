SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

# We show that a sequence can be losslessly regenerated from its linked dBG
# lossless0: random sequence in a single read is regenerated
# lossless1: minimum required read length to reproduce exact contig

all:
	cd lossless0 && $(MAKE)
	cd lossless1 && $(MAKE)
	@echo "All looks good."

clean:
	cd lossless0 && $(MAKE) clean
	cd lossless1 && $(MAKE) clean

.PHONY: all clean
