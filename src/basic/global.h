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

typedef signed char boolean;

#ifndef true
#define true 1
#define false 0
#endif

#define CTXVERSIONSTR "0.0"

#define QUOTE_MACRO(str) #str
#define QUOTE_VALUE(str) QUOTE_MACRO(str)

#define SWAP(x,y,tmp) ((tmp) = (x), (x) = (y), (y) = (tmp))

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))

#define ABSDIFF(a,b) ((a) > (b) ? (a)-(b) : (b)-(a))
#define MEDIAN(arr,len) (len&0x1?arr[len/2]:(arr[(len-1)/2]+arr[len/2])/2.0)

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) (0x1UL << (64 - __builtin_clzl(x)))
#endif

#define SAFE_ADD(dst,add,max) do {                 \
  size_t _full = ((max)-(dst))<(add);              \
  (dst) = (_full)*(max) + (!_full)*((dst)+(add));  \
} while(0)

#define rot32(x,r) (((x)<<(r)) | ((x)>>(32-(r))))
#define rot64(x,r) (((x)<<(r)) | ((x)>>(64-(r))))

// Unaligned memory access
// #define NO_UNALIGNED_ACCESS 1
// src, tgt are both ptrs
#ifdef NO_UNALIGNED_ACCESS
  #define unaligned_16(tgt,src) memcpy(tgt, src, 2)
  #define unaligned_32(tgt,src) memcpy(tgt, src, 4)
  #define unaligned_64(tgt,src) memcpy(tgt, src, 8)
#else
  #define unaligned_16(tgt,src) (*tgt = *(uint16_t*)(src))
  #define unaligned_32(tgt,src) (*tgt = *(uint32_t*)(src))
  #define unaligned_64(tgt,src) (*tgt = *(uint64_t*)(src))
#endif

// Use bit sets with word size 64
#define round_bits_to_bytes(bits)   (((bits)+7)/8)
#define round_bits_to_words64(bits) (((bits)+63)/64)

#define bitset2_has(arr,idx,offset)    (((arr)[idx] >> (offset)) & 0x1UL)
#define bitset2_set(arr,idx,offset)     ((arr)[idx] |= (0x1UL << (offset)))
#define bitset2_del(arr,idx,offset)     ((arr)[idx] &=~(0x1UL << (offset)))
#define bitset2_cpy(arr,idx,offset,bit) ((arr)[idx] |= ((uint64_t)(bit) << (offset)))

#define bitset_has(arr,pos) \
        bitset2_has(arr, (pos)/(sizeof(*(arr))*8), (pos)%(sizeof(*(arr))*8))
#define bitset_set(arr,pos) \
        bitset2_set(arr, (pos)/(sizeof(*(arr))*8), (pos)%(sizeof(*(arr))*8))
#define bitset_del(arr,pos) \
        bitset2_del(arr, (pos)/(sizeof(*(arr))*8), (pos)%(sizeof(*(arr))*8))
#define bitset_cpy(arr,pos,bit) \
        bitset2_cpy(arr, (pos)/(sizeof(*(arr))*8), (pos)%(sizeof(*(arr))*8), (bit))

#define bitset_clear_word(arr,pos) ((arr)[(pos) / (sizeof(*(arr))*8)] = 0)

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

// A function for standard output
void message(const char *fmt, ...)
__attribute__((format(printf, 1, 2)));

void print_usage(const char *msg, const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 2, 3)));

#endif /* GLOBAL_H_ */
