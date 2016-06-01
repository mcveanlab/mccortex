SHELL=/bin/bash -euo pipefail

# Test bubble popping
# Each test case is in a separate sub-directory

all:
	cd pop_bubbles1 && $(MAKE)
	cd pop_bubbles2 && $(MAKE)
	@echo "pop_bubbles: All looks good."

clean:
	cd pop_bubbles1 && $(MAKE) clean
	cd pop_bubbles2 && $(MAKE) clean

.PHONY: all clean
