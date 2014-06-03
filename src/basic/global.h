#ifndef GLOBAL_H_
#define GLOBAL_H_

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE
// #define _GNU_SOURCE

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

// Makefile sets:
//  NDEBUG=1         Turns off debugging
//  CTXCHECKS=1      Turns on heavy checks
//  MIN_KMER_SIZE    Min kmer-size compiled e.g. 3 for maxk=31, 33 for maxk=63
//  MAX_KMER_SIZE    Max kmer-size compiled e.g. 31 for maxk=31, 63 for maxk=63
//  USE_CITY_HASH=1  Use Google's CityHash instead of Bob Jenkin's lookup3

#define MAX_IO_THREADS 10

#if defined(CTXCHECKS) && CTXCHECKS == 0
#  undef CTXCHECKS
#endif

#ifdef NDEBUG
#  define ASSERTSTR "ASSERTS=OFF"
#else
#  define ASSERTSTR "ASSERTS=ON"
#endif
#ifdef CTXCHECKS
#  define CHECKSTR "CHECKS=ON"
#else
#  define CHECKSTR "CHECKS=OFF"
#endif

#include "htslib/version.h"
#define LIBS_VERSION "zlib="ZLIB_VERSION" htslib="HTS_VERSION
#define VERSION_STATUS_STR "ctx="CTX_VERSION" "LIBS_VERSION" " \
                           ASSERTSTR" "CHECKSTR" " \
                           "k="QUOTE_VALUE(MIN_KMER_SIZE)".."QUOTE_VALUE(MAX_KMER_SIZE)

//
// Dynamic memory allocation with checks
//

volatile size_t ctx_num_allocs, ctx_num_frees;

// Macros for memory management
// `ptr` can be NULL
#define ctx_malloc(mem) call_malloc(mem,__FILE__,__func__,__LINE__)
#define ctx_calloc(nel,elsize) call_calloc(nel,elsize,__FILE__,__func__,__LINE__)
#define ctx_realloc(ptr,mem) call_realloc(ptr,mem,__FILE__,__func__,__LINE__)
#define ctx_recalloc(ptr,old,new) call_recalloc(ptr,old,new,__FILE__,__func__,__LINE__)
#define ctx_free(ptr) call_free(ptr)

// Memory allocation functions
// `ptr` can be NULL
void* call_malloc(size_t mem, const char *file, const char *func, int line);
void* call_calloc(size_t nel, size_t elsize, const char *file, const char *func, int line);
void* call_realloc(void *ptr, size_t mem, const char *file, const char *func, int line);

// Resize memory, zero new memory
void* call_recalloc(void *ptr, size_t oldsize, size_t newsize,
                    const char *file, const char *func, int line);

void call_free(void *ptr);

//
// Internal Integrity Checks: ctx_check(), ctx_assert(), ctx_assume()
//
// ctx_assert() behaves like assert()
// ctx_assume() behaves like assert() when NDEBUG=1,
//              otherwise interpret as a truth statement for the compiler
// ctx_check()  behaves like assert() when CTXCHECKS=1 but for heavy checks
//              without CTXCHECKS, does nothing
//
//            | NDEBUG=1 |      Default
//----------------------------------------------
// ctx_assert | nothing  |  fast check + abort()
// ctx_assume | optimise |  fast check + abort()
//
//            | Default  |     CTXCHECKS=1
//----------------------------------------------
// ctx_check  | nothing  |  slow check + abort()

void call_assert(const char *file, const char *func, int line,
                 const char *assert, const char *fmt, ...)
__attribute__((format(printf, 5, 6)))
__attribute__((noreturn));

void call_assert_no_abort(const char *file, const char *func, int line,
                          const char *assert, const char *fmt, ...)
__attribute__((format(printf, 5, 6)));

// ctx_assert()
#ifdef NDEBUG
  #define ctx_assert2(x,msg,...) do {} while(0)
#else
  #define ctx_assert2(x,msg,...) \
((x) ? (void)0 : call_assert(__FILE__,__func__,__LINE__,QUOTE_VALUE(x),msg,##__VA_ARGS__))
#endif

// ctx_assume(): tells the compiler a condition that always holds
#ifdef NDEBUG
  #define ctx_assume2(x,msg,...) do { if(x) (void)0; else __builtin_unreachable(); } while(0)
#else
  #define ctx_assume2(x,msg,...) ctx_assert2(x,msg,##__VA_ARGS__)
#endif

// ctx_check():
#ifdef CTXCHECKS
  #define ctx_check2(x,msg,...) ctx_assert2(x,msg,##__VA_ARGS__)
#else
  #define ctx_check2(x,msg,...) do {} while(0)
#endif

// Check is turned on with CTXCHECKS=1 -> heavy lifting involved
// assert -> no action if NDEBUG=1
// assume -> declares !x impossible (helps with optimisations)
#define ctx_check(x)  ctx_check2(x,NULL)
#define ctx_assert(x) ctx_assert2(x,NULL)
#define ctx_assume(x) ctx_assume2(x,NULL)

// Return false if a condition fails, rather than aborting
// Note: doesn't depend on CTXCHECKS
#define check_ret2(x,msg,...) do { if(!(x)) {                               \
  call_assert_no_abort(__FILE__,__func__,__LINE__,QUOTE_VALUE(x),msg,##__VA_ARGS__);\
  return false;                                                                \
}} while(0)

#define check_ret(x) check_ret2(x,NULL)

//
// Exit with an error, give a warning: die() / warn()
//

// Output destination and lock for output
// set to NULL to turn off message printing
extern FILE *ctx_msg_out;
extern pthread_mutex_t biglock;

#define die(fmt, ...) call_die(__FILE__, __func__, __LINE__, fmt, ##__VA_ARGS__)
#define warn(fmt, ...) call_warn(__FILE__, __func__, __LINE__, fmt, ##__VA_ARGS__)

void call_die(const char *file, const char *func, int line, const char *fmt, ...)
__attribute__((format(printf, 4, 5)))
__attribute__((noreturn));

void call_warn(const char *file, const char *func, int line, const char *fmt, ...)
__attribute__((format(printf, 4, 5)));

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

// Swap macro swaps a byte at a time
#define SWAPCPY(x,y) do {                                             \
  char *_a = &(x), *_b = &(y), *_end = _a + sizeof(x), _tmp;          \
  for(; _a < _end; _a++, _b++) { _tmp = *_a; *_a = *_b; *_b = _tmp; } \
} while(0)

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))

#define ABSDIFF(a,b) ((a) > (b) ? (a)-(b) : (b)-(a))
#define MEDIAN(arr,len) \
        (!(len)?0:((len)&1?(arr)[(len)/2]:((arr)[(len)/2-1]+(arr)[(len)/2])/2.0))

// Number of reads to hold in the msg pool
#define MSGPOOLSIZE 2048
#define USE_MSG_POOL MSGP_LOCK_MUTEX

// MSGP_LOCK_SPIN
// MSGP_LOCK_YIELD

#endif /* GLOBAL_H_ */
