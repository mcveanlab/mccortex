#include "global.h"
#include "util.h"

// uncomment next line to silence most (stdout) output from cortex
//#define GLOBAL_MSG_SILENT

void* ctx_malloc(size_t mem, const char *file, int line)
{
  void *ptr = malloc(mem);
  if(ptr == NULL) {
    char memstr[100];
    bytes_to_str(mem, 1, memstr);
    call_die(file, line, "Out of memory (malloc %s)", memstr);
  }
  return ptr;
}

void* ctx_calloc(size_t nel, size_t elsize, const char *file, int line)
{
  void *ptr = calloc(nel, elsize);
  if(ptr == NULL) {
    char nelstr[100], elsizestr[100], memstr[100];
    ulong_to_str(nel, nelstr);
    bytes_to_str(elsize, 1, elsizestr);
    bytes_to_str(nel * elsize, 1, memstr);
    call_die(file, line, "Out of memory (calloc %s x %s = %s)",
             nelstr, elsizestr, memstr);
  }
  return ptr;
}

void* ctx_realloc(void *ptr, size_t mem, const char *file, int line)
{
  void *ptr2 = realloc(ptr, mem);
  if(ptr2 == NULL) {
    char memstr[100];
    bytes_to_str(mem, 1, memstr);
    call_die(file, line, "Out of memory (realloc %s)", memstr);
  }
  return ptr2;
}

void call_die(const char *file, int line, const char *fmt, ...)
{
  fflush(stdout);

  // Print error
  fprintf(stderr, "[%s:%i] Error: ", file, line);

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt + strlen(fmt) - 1) != '\n')
  {
    fprintf(stderr, "\n");
  }

  exit(EXIT_FAILURE);
}

void call_warn(const char *file, int line, const char *fmt, ...)
{
  fflush(stdout);

  // Print warning
  fprintf(stderr, "[%s:%i] Warning: ", file, line);

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt + strlen(fmt) - 1) != '\n')
  {
    fprintf(stderr, "\n");
  }

  fflush(stderr);
}

// A function for standard output
void message(const char *fmt, ...)
{
#ifdef GLOBAL_MSG_SILENT
  (void)fmt;
#else
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stdout, fmt, argptr);
  va_end(argptr);
  fflush(stdout);
#endif
}

void print_usage(const char *msg, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
  }

  fputs(msg, stderr);
  exit(EXIT_FAILURE);
}
