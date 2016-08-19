SHELL=/bin/bash -euo pipefail

# Each test case is in a separate sub-directory

# pjoin0:
# pjoin1:

all:
	cd pjoin0 && $(MAKE)
	cd pjoin1 && $(MAKE)
	@echo "All looks good."

clean:
	cd pjoin0 && $(MAKE) clean
	cd pjoin1 && $(MAKE) clean

.PHONY: all clean
