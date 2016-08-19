SHELL=/bin/bash -euo pipefail

# Test bubble calling and decomposing to VCF
# Each test case is in a separate sub-directory

all:
	cd bubbles1 && $(MAKE)
	cd bubbles2 && $(MAKE)
	cd bubbles3 && $(MAKE)
	cd bubbles4 && $(MAKE)
	cd bubbles5 && $(MAKE)
	@echo "bubbles: All looks good."

clean:
	cd bubbles1 && $(MAKE) clean
	cd bubbles2 && $(MAKE) clean
	cd bubbles3 && $(MAKE) clean
	cd bubbles4 && $(MAKE) clean
	cd bubbles5 && $(MAKE) clean

.PHONY: all clean
