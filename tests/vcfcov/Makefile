SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

all:
	cd calls0 && $(MAKE)
	cd calls1 && $(MAKE)
	@echo "All looks good."

clean:
	cd calls0 && $(MAKE) clean
	cd calls1 && $(MAKE) clean

.PHONY: all clean view
