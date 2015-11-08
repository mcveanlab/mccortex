SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

all:
	cd bubbles1 && $(MAKE) clean && $(MAKE)
	cd bubbles2 && $(MAKE) clean && $(MAKE)
	@echo "All looks good."

clean:
	cd bubbles1 && $(MAKE) clean
	cd bubbles2 && $(MAKE) clean

.PHONY: all clean view
