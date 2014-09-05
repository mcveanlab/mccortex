NWORDS=1
SHELL:=/bin/bash

all: revcmp

clean:
	rm -rf revcmp

revcmp: revcmp.c
	$(CC) -O4 -Wall -Wextra -DNUM_BKMER_WORDS=$(NWORDS) -o $@ $<

profile:
	for i in {1..5}; do for m in {0..4}; do time ./revcmp -m $$m -n 1000000000; done; done

.PHONY: all clean profile
