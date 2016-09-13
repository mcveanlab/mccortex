#ifndef GLOBAL_H_
#define GLOBAL_H_

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE
// #define _GNU_SOURCE

// Request PRIu64 etc. from inttypes.h
#define __STDC_FORMAT_MACROS

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <stddef.h> // defines ptrdiff_t
#include <string.h>
#include <strings.h> // strcasecmp
#include <inttypes.h> // fixed size integers e.g. uint64_t
#include <stdbool.h> // defines bool, true, false
#include <limits.h>
#include <pthread.h>
#include <zlib.h>
#include <errno.h>

#include "bit_array/bit_macros.h" // from bit_array
#include "string_buffer/string_buffer.h"

#include "version.h" // defines CTX_VERSION
#include "cortex_types.h" // common basic data types
#include "kmer_size.h"

// Makefile sets:
//  NDEBUG=1         Turns off debugging
//  CTXCHECKS=1      Turns on heavy checks
//  MIN_KMER_SIZE    Min kmer-size compiled e.g. 3 for maxk=31, 33 for maxk=63
//  MAX_KMER_SIZE    Max kmer-size compiled e.g. 31 for maxk=31, 63 for maxk=63
//  USE_CITY_HASH=1  Use Google's CityHash instead of Bob Jenkin's lookup3
//  USE_XXHASH=1     Use xxHash instead of Bob Jenkin's lookup3

#define ONE_MEGABYTE (1<<20)
#define MAX_IO_THREADS 10
#define DEFAULT_IO_BUFSIZE (4*ONE_MEGABYTE)
#define MCCORTEX_URL "https://github.com/mcveanlab/mccortex"

#include "ctx_assert.h"
#include "ctx_alloc.h" // Wrappers for malloc, calloc etc.
#include "ctx_output.h" // Printing status messages

#include "htslib/version.h"
#define LIBS_VERSION "zlib="ZLIB_VERSION" htslib="HTS_VERSION

// Must include hash.h to use this!
#define VERSION_STATUS_STR "mccortex="CTX_VERSION" "LIBS_VERSION" "ASSERTSTR" hash="HASH_NAME_STR" "CHECKSTR

//
// Setup / clear up of library functions
//
void cortex_init();
void cortex_destroy();

//
// Common MACROs
//

#define QUOTE_MACRO(str) #str
#define QUOTE_VALUE(str) QUOTE_MACRO(str)

#define SWAP(x,y) do { __typeof(x) _tmp = (x); (x) = (y); (y) = _tmp; } while(0)

// Swap macro swaps a byte at a time. x and y must not be overlapping.
#define SWAPCPY(x,y) do {                                             \
  char *_a = &(x), *_b = &(y), *_end = _a + sizeof(x), _tmp;          \
  for(; _a < _end; _a++, _b++) { _tmp = *_a; *_a = *_b; *_b = _tmp; } \
} while(0)

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))

#define ABSDIFF(a,b) ((a) > (b) ? (a)-(b) : (b)-(a))

#define cmp(a,b) (((a) > (b)) - ((b) > (a)))

// Number of reads to hold in the msg pool
#define MSGPOOLSIZE 2048
#define USE_MSG_POOL MSGP_LOCK_MUTEX

// MSGP_LOCK_SPIN
// MSGP_LOCK_YIELD

#endif /* GLOBAL_H_ */
