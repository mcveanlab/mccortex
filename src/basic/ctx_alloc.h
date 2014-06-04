#ifndef CTX_ALLOC_H_
#define CTX_ALLOC_H_

#include <stdlib.h>

//
// Dynamic memory allocation with checks
//

// Macros for memory management
// `ptr` can be NULL
#define ctx_malloc(mem) alloc_malloc(mem,__FILE__,__func__,__LINE__)
#define ctx_calloc(nel,elsize) alloc_calloc(nel,elsize,__FILE__,__func__,__LINE__)
#define ctx_realloc(ptr,mem) alloc_realloc(ptr,mem,__FILE__,__func__,__LINE__)
#define ctx_reallocarray(ptr,nel,elsize) alloc_reallocarray(ptr,nel,elsize,__FILE__,__func__,__LINE__)
#define ctx_recallocarray(ptr,oldnel,newnel,elsize) alloc_recallocarray(ptr,oldnel,newnel,elsize,__FILE__,__func__,__LINE__)
#define ctx_free(ptr) alloc_free(ptr)

// Memory allocation functions
// `ptr` can be NULL
void* alloc_malloc(size_t mem, const char *file, const char *func, int line);
void* alloc_calloc(size_t nel, size_t elsize, const char *file, const char *func, int line);
void* alloc_realloc(void *ptr, size_t mem, const char *file, const char *func, int line);

void* alloc_reallocarray(void *ptr, size_t nel, size_t elsize,
                        const char *file, const char *func, int line);

// Resize memory, zero new memory
void* alloc_recallocarray(void *ptr, size_t oldnel, size_t newnel, size_t elsize,
                          const char *file, const char *func, int line);

void alloc_free(void *ptr);

size_t alloc_get_num_allocs();
size_t alloc_get_num_frees();

#endif /* CTX_ALLOC_H_ */
