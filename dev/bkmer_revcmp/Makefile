NWORDS=1

all: revcmp

clean:
	rm -rf revcmp

revcmp: revcmp.c
	$(CC) -O4 -Wall -Wextra -DNUM_BKMER_WORDS=$(NWORDS) -o $@ $<

.PHONY: all clean
