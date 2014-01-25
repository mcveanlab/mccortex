#ifndef GLOBAL_H_
#define GLOBAL_H_

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <stddef.h> // defines ptrdiff_t
#include <string.h>
#include <strings.h> // strcasecmp
#include <inttypes.h>
#include <limits.h>
#include <assert.h>
#include <zlib.h>
#include <pthread.h>

#include "bit_macros.h" // from bit_array

// #define CTXVERSIONSTR "0.0"
#include "version.h"
#include "cortex_types.h"

// Number of reads to hold in the msg pool
#define MSGPOOLRSIZE 20

// set to NULL to turn off message printing
extern FILE *ctx_msg_out;
extern pthread_mutex_t biglock;

#define QUOTE_MACRO(str) #str
#define QUOTE_VALUE(str) QUOTE_MACRO(str)

#define SWAP(x,y,tmp) ((tmp) = (x), (x) = (y), (y) = (tmp))

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))

#define ABSDIFF(a,b) ((a) > (b) ? (a)-(b) : (b)-(a))
#define MEDIAN(arr,len) \
        (!(len)?0:((len)&1?(arr)[(len)/2]:((arr)[(len)/2-1]+(arr)[(len)/2])/2.0))

// dynamic memory allocation with checks
#define malloc2(mem) ctx_malloc(mem,__FILE__,__LINE__)
#define calloc2(nel,elsize) ctx_calloc(nel,elsize,__FILE__,__LINE__)
#define realloc2(ptr,mem) ctx_realloc(ptr,mem,__FILE__,__LINE__)

void* ctx_malloc(size_t mem, const char *file, int line);
void* ctx_calloc(size_t nel, size_t elsize, const char *file, int line);
void* ctx_realloc(void *ptr, size_t mem, const char *file, int line);

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define warn(fmt, ...) call_warn(__FILE__, __LINE__, fmt, ##__VA_ARGS__)

void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

void call_warn(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)));

void fmessage(FILE *fh, const char *fmt, ...)
__attribute__((format(printf, 2, 3)));

void ftimestamp(FILE *fh);

void fstatus(FILE *fh, const char *fmt, ...)
__attribute__((format(printf, 2, 3)));

#define message(fmt,...) fmessage(ctx_msg_out,fmt, ##__VA_ARGS__)
#define timestamp()      ftimestamp(ctx_msg_out)
#define status(fmt,...)  fstatus(ctx_msg_out,fmt, ##__VA_ARGS__)

void print_usage(const char *msg, const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 2, 3)));

void cortex_init();
void cortex_destroy();

#endif /* GLOBAL_H_ */
