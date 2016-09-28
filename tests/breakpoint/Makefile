SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

# breakpoint0: empty breakpoint calls (K=31)
# breakpoint1: test we get the correct VCF (K=11)
# breakpoint2: check we don't call any ref bubbles

all:
	cd breakpoint0 && $(MAKE)
	cd breakpoint1 && $(MAKE)
	cd breakpoint2 && $(MAKE)
	@echo "All looks good."

clean:
	cd breakpoint0 && $(MAKE) clean
	cd breakpoint1 && $(MAKE) clean
	cd breakpoint2 && $(MAKE) clean

.PHONY: all clean
