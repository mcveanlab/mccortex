#include "global.h"
#include "ctx_assert.h"

//
// Checks and asserts
//

static void ctx_assertf2(const char *file, const char *func, int line,
                         const char *asserttxt, const char *fmt, va_list argptr)
{
  pthread_mutex_lock(&ctx_biglock);
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Assert Failed %s(): %s", file, line, func, asserttxt);

  if(fmt != NULL) {
    fputs(": ", stderr);
    vfprintf(stderr, fmt, argptr);
  }

  // Print a timestamp so we know when the crash occurred
  fprintf(stderr, "\n");
  timestampf(stderr);
  fputs(" Assert Error\n", stderr);
  fflush(stderr);
  pthread_mutex_unlock(&ctx_biglock);
}

void ctx_assertf_no_abort(const char *file, const char *func, int line,
                          const char *asserttxt, const char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  ctx_assertf2(file, func, line, asserttxt, fmt, argptr);
  va_end(argptr);
}

void ctx_assertf(const char *file, const char *func, int line,
                 const char *asserttxt, const char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  ctx_assertf2(file, func, line, asserttxt, fmt, argptr);
  va_end(argptr);
  abort();
}
