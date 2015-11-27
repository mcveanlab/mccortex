SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

# build0: random sequence, sort graph, reassemble sequence
# build1: test --intersection and --graph arguments 

all:
	cd build0 && $(MAKE)
	cd build1 && $(MAKE)
	@echo "All looks good."

clean:
	cd build0 && $(MAKE) clean
	cd build1 && $(MAKE) clean

.PHONY: all clean
