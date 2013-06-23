#include "global.h"

#include <stddef.h> // defines ptrdiff_t

char print_debug = 0;

// uncomment next line to silence most (stdout) output from cortex
//#define GLOBAL_MSG_SILENT

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
