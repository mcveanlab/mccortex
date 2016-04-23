SHELL=/bin/bash -euo pipefail

# Test bubble calling and decomposing to VCF
# Each test case is in a separate sub-directory

all:
	cd bubbles1 && $(MAKE)
	cd bubbles2 && $(MAKE)
	cd bubbles3 && $(MAKE)
	@echo "All looks good."

clean:
	cd bubbles1 && $(MAKE) clean
	cd bubbles2 && $(MAKE) clean
	cd bubbles3 && $(MAKE) clean

.PHONY: all clean
