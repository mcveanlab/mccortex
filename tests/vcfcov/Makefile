SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

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
