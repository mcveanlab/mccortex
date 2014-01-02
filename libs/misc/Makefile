CFLAGS=-Wall -Wextra -O3

all: city.o mem_size.o

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $*.c -o $*.o

clean:
	rm -rf *.o

.PNOHY: all clean
