all: city.o lookup3.o

CFLAGS=-Wall -Wextra -O3

city.o: city.c city.h
	$(CC) $(CFLAGS) -c city.c -o city.o

lookup3.o: lookup3.c lookup3.h
	$(CC) $(CFLAGS) -c lookup3.c -o lookup3.o

clean:
	rm -rf *.o

.PNOHY: all clean
