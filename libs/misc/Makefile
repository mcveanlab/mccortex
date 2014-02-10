CFLAGS=-Wall -Wextra -O3

# These files are all compiled by the main Makefile,
# this is just here as documentation

all: city.o mem_size.o

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $*.c -o $*.o

clean:
	rm -rf *.o

.PNOHY: all clean
