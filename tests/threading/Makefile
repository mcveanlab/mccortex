SHELL=/bin/bash -euo pipefail

# Test read threading through the graph
# Each test case is in a separate sub-directory
# threading1:
# threading2: paired-end threading
# threading3: paired-end threading with short reads
# threading4:

all:
	cd threading1 && $(MAKE)
	cd threading2 && $(MAKE)
	cd threading3 && $(MAKE)
	cd threading4 && $(MAKE)
	@echo "threading: All looks good."

clean:
	cd threading1 && $(MAKE) clean
	cd threading2 && $(MAKE) clean
	cd threading3 && $(MAKE) clean
	cd threading4 && $(MAKE) clean

.PHONY: all clean
