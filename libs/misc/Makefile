# These files are all compiled by the main Makefile,
# this is just here as documentation

PLATFORM := $(shell uname)
COMPILER := $(shell ($(CC) -v 2>&1) | tr A-Z a-z )

# clang Link Time Optimisation (lto) seems to have issues atm
ifneq (,$(findstring clang,$(COMPILER)))
	OPT = -O3
else
	OPT = -O4 -flto
endif

CFLAGS=-Wall -Wextra $(OPT)

all: city.o mem_size.o

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $*.c -o $*.o

clean:
	rm -rf *.o

.PNOHY: all clean
