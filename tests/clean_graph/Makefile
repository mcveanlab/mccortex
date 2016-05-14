SHELL=/bin/bash -euo pipefail

# Test graph cleaning
# Each test case is in a separate sub-directory

all:
	cd clean1 && $(MAKE)
	cd clean2 && $(MAKE)
	cd clean3 && $(MAKE)
	@echo "clean_graph: All looks good."

clean:
	cd clean1 && $(MAKE) clean
	cd clean2 && $(MAKE) clean
	cd clean3 && $(MAKE) clean

.PHONY: all clean
