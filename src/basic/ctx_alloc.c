#include "global.h"
#include "ctx_alloc.h"
#include "util.h"

static volatile size_t ctx_num_allocs = 0, ctx_num_frees = 0;

static inline void _oom(void *ptr, size_t nel, size_t elsize,
                        const char *file, const char *func, int line)
__attribute__((noreturn));

static inline void _oom(void *ptr, size_t nel, size_t elsize,
                        const char *file, const char *func, int line)
{
  char memstr[50];
  if(SIZE_MAX / elsize < nel) strcpy(memstr, "overflow");
  else bytes_to_str(nel*elsize, 1, memstr);
  dief(file, func, line, "Out of memory (%p, %zu x %zu = %s)",
       ptr, nel, elsize, memstr);
}

// If `zero` is true and ptr is NULL call calloc, otherwise realloc
void* alloc_mem(void *ptr, size_t nel, size_t elsize, bool zero,
                const char *file, const char *func, int line)
{
  void *ptr2;
  if(nel && elsize && SIZE_MAX / elsize < nel)
    _oom(ptr, nel, elsize, file, func, line);

  if(ptr || !zero)
    ptr2 = realloc(ptr, nel*elsize);
  else
    ptr2 = calloc(nel, elsize);

  if(ptr2 == NULL) _oom(ptr, nel, elsize, file, func, line);
  if(ptr == NULL) __sync_add_and_fetch(&ctx_num_allocs, 1); // ++ctx_num_allocs

  return ptr2;
}

// Resize memory, zero new memory
void* alloc_recallocarray(void *ptr, size_t oldnel, size_t newnel, size_t elsize,
                          const char *file, const char *func, int line)
{
  void *ptr2 = alloc_mem(ptr, newnel, elsize, true, file, func, line);
  if(ptr != NULL && newnel > oldnel)
    memset((char*)ptr2+oldnel*elsize, 0, (newnel-oldnel)*elsize);
  return ptr2;
}

// `ptr` can be NULL
void alloc_free(void *ptr)
{
  free(ptr);
  if(ptr != NULL)
    __sync_add_and_fetch(&ctx_num_frees, 1); // ++ctx_num_frees
}

size_t alloc_get_num_allocs()
{
  return (size_t)ctx_num_allocs;
}

size_t alloc_get_num_frees()
{
  return (size_t)ctx_num_frees;
}
