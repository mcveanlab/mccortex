#ifndef CTX_ALLOC_H_
#define CTX_ALLOC_H_

#include <stdlib.h>

//
// Wrappers for dynamic memory allocation functions
//
// We do this to provide counters of all alloc+free calls so we can report
// if there are any memory leaks on exit. It also allows us to check if malloc()
// etc. have returned NULL and exit with an informative message with line number
// of offending call.
//

// Macros for memory management
// `ptr` can be NULL
#define ctx_malloc(mem) alloc_mem(NULL,1,mem,false,__FILE__,__func__,__LINE__)
#define ctx_calloc(nel,elsize) alloc_mem(NULL,nel,elsize,true,__FILE__,__func__,__LINE__)
#define ctx_realloc(ptr,mem) alloc_mem(ptr,1,mem,false,__FILE__,__func__,__LINE__)
#define ctx_reallocarray(ptr,nel,elsize) alloc_mem(ptr,nel,elsize,false,__FILE__,__func__,__LINE__)
#define ctx_recallocarray(ptr,oldnel,newnel,elsize) alloc_recallocarray(ptr,oldnel,newnel,elsize,__FILE__,__func__,__LINE__)
#define ctx_free(ptr) alloc_free(ptr)

// Allocate / reallocate memory. `ptr` can be NULL
// Prints error message and calls exit() if out of memory / cannot alloc
void* alloc_mem(void *ptr, size_t nel, size_t elsize, bool zero,
                const char *file, const char *func, int line);

// Allocate / resize memory, ensure all new memory is zero'ed
void* alloc_recallocarray(void *ptr, size_t oldnel, size_t newnel, size_t elsize,
                          const char *file, const char *func, int line);

// Free allocated memory, `ptr` is allowed to be NULL
void alloc_free(void *ptr);

// Get number of allocations / frees
size_t alloc_get_num_allocs();
size_t alloc_get_num_frees();

#endif /* CTX_ALLOC_H_ */
