#include "global.h"
#include "ctx_alloc.h"
#include "util.h"

static volatile size_t ctx_num_allocs = 0, ctx_num_frees = 0;

void* alloc_malloc(size_t mem, const char *file, const char *func, int line)
{
  void *ptr = malloc(mem);
  if(ptr == NULL) {
    char memstr[100];
    bytes_to_str(mem, 1, memstr);
    dief(file, func, line, "Out of memory (malloc %s)", memstr);
  }
  __sync_add_and_fetch(&ctx_num_allocs, 1); // ++ctx_num_allocs
  return ptr;
}

void* alloc_calloc(size_t nel, size_t elsize,
                  const char *file, const char *func, int line)
{
  void *ptr = calloc(nel, elsize);
  if(ptr == NULL) {
    char nelstr[100], elsizestr[100], memstr[100];
    ulong_to_str(nel, nelstr);
    bytes_to_str(elsize, 1, elsizestr);
    bytes_to_str(nel * elsize, 1, memstr);
    dief(file, func, line, "Out of memory (calloc %s x %s = %s)",
         nelstr, elsizestr, memstr);
  }
  __sync_add_and_fetch(&ctx_num_allocs, 1); // ++ctx_num_allocs
  return ptr;
}

void* alloc_reallocarray(void *ptr, size_t nel, size_t elsize,
                        const char *file, const char *func, int line)
{
  void *ptr2;
  if(SIZE_MAX / elsize < nel || (ptr2 = realloc(ptr, nel*elsize)) == NULL)
  {
    char nelstr[50], elsizestr[50], memstr[50];
    bytes_to_str(nel, 1, nelstr);
    bytes_to_str(elsize, 1, elsizestr);
    if(SIZE_MAX / elsize < nel) strcpy(memstr, "overflow");
    else bytes_to_str(nel*elsize, 1, memstr);
    dief(file, func, line, "Out of memory (reallocarray %sx%s=%s)",
         nelstr, elsizestr, memstr);
  }
  if(ptr == NULL) __sync_add_and_fetch(&ctx_num_allocs, 1); // ++ctx_num_allocs
  return ptr2;
}

// Resize memory
void* alloc_realloc(void *ptr, size_t mem,
                   const char *file, const char *func, int line)
{
  return alloc_reallocarray(ptr, 1, mem, file, func, line);
}

// Resize memory, zero new memory
void* alloc_recallocarray(void *ptr, size_t oldnel, size_t newnel, size_t elsize,
                          const char *file, const char *func, int line)
{
  if(ptr == NULL) return alloc_calloc(newnel, elsize, file, func, line);
  else {
    ptr = alloc_reallocarray(ptr, newnel, elsize, file, func, line);
    if(newnel > oldnel)
      memset((char*)ptr+oldnel*elsize, 0, (newnel-oldnel)*elsize);
    return ptr;
  }
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
