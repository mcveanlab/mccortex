#include "global.h"
#include "all_tests.h"

// Common functions here
FILE *ctx_tst_out = NULL;

// A function for standard output
void test_status(const char *fmt, ...)
{
  va_list argptr;

  if(ctx_tst_out != NULL) {
    pthread_mutex_lock(&biglock);
    timestamp(ctx_tst_out);
    if(fmt[0] != ' ' && fmt[0] != '[') fputc(' ', ctx_tst_out);
    va_start(argptr, fmt);
    vfprintf(ctx_tst_out, fmt, argptr);
    va_end(argptr);
    if(fmt[strlen(fmt)-1] != '\n') fputc('\n', ctx_tst_out);
    fflush(ctx_tst_out);
    pthread_mutex_unlock(&biglock);
  }
}
