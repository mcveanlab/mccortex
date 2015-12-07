#include "global.h"
#include "ctx_output.h"
#include "util.h"

#include <libgen.h> // basename

FILE *ctx_msg_out = NULL;
pthread_mutex_t ctx_biglock;
char ctx_cmdcode[4] = "000";

// Get output lock before calling this function
static void err_msg(const char *type, const char *path, const char *func,
                    int line, const char *fmt, va_list argptr)
{
  (void)func;
  // Get filename from path
  char *filename = strdup(path);
  fflush(stdout);
  fflush(stderr);
  timestampf(stderr);
  fprintf(stderr, "[%s:%i] %s: ", basename(filename), line, type);
  vfprintf(stderr, fmt, argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  free(filename);
}

void dief(const char *file, const char *func, int line, const char *fmt, ...)
{
  va_list argptr;
  pthread_mutex_lock(&ctx_biglock); // we never release this lock
  va_start(argptr, fmt);
  err_msg("Fatal Error", file, func, line, fmt, argptr);
  va_end(argptr);
  exit(EXIT_FAILURE);
  // abort();
}

void warnf(const char *file, const char *func, int line, const char *fmt, ...)
{
  va_list argptr;
  pthread_mutex_lock(&ctx_biglock);
  va_start(argptr, fmt);
  err_msg("Warn", file, func, line, fmt, argptr);
  va_end(argptr);
  pthread_mutex_unlock(&ctx_biglock);
}

// A function for standard output
void messagef(FILE *fh, const char *fmt, ...)
{
  if(fh != NULL) {
    va_list argptr;
    va_start(argptr, fmt);
    vfprintf(fh, fmt, argptr);
    va_end(argptr);
    fflush(fh);
  }
}

void timestampf(FILE *fh)
{
  time_t t;
  char tstr[200];

  if(fh != NULL) {
    time(&t);
    strftime(tstr, 100, "[%d %b %Y %H:%M:%S", localtime(&t));
    fprintf(fh, "%s-%s]", tstr, ctx_cmdcode);
  }
}

// A function for standard output
void statusf(FILE *fh, const char *fmt, ...)
{
  va_list argptr;

  if(fh != NULL) {
    pthread_mutex_lock(&ctx_biglock);
    timestampf(fh);
    if(fmt[0] != ' ' && fmt[0] != '[') fputc(' ', fh);
    va_start(argptr, fmt);
    vfprintf(fh, fmt, argptr);
    va_end(argptr);
    if(fmt[strlen(fmt)-1] != '\n') fputc('\n', fh);
    fflush(fh);
    pthread_mutex_unlock(&ctx_biglock);
  }
}

void print_usage(const char *msg, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    pthread_mutex_lock(&ctx_biglock);
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fputc('\n', stderr);
  }

  fputs(msg, stderr);
  exit(EXIT_FAILURE);
}

void ctx_output_init()
{
  if(!ctx_msg_out) ctx_msg_out = stderr;

  static const char consonants[]
    = "bcdfghjklmnpqrstvwxyzBCDFGHJKLMNPQRSTVWXYZ";
  static const char vowels[] = "aeiouAEIOU";

  if(pthread_mutex_init(&ctx_biglock, NULL) != 0) {
    printf("%s:%i: mutex init failed\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // set unique code for this output
  ctx_cmdcode[0] = consonants[rand() % strlen(consonants)];
  ctx_cmdcode[1] = vowels[rand() % strlen(vowels)];
  ctx_cmdcode[2] = consonants[rand() % strlen(consonants)];
  ctx_cmdcode[3] = '\0';
}

void ctx_output_destroy()
{
  pthread_mutex_destroy(&ctx_biglock);
}

void ctx_update(const char *job_name, size_t niter)
{
  if(niter % CTX_UPDATE_REPORT_RATE == 0)
  {
    char num_str[100];
    long_to_str(niter, num_str);
    status("[%s] Read %s entries (reads / read pairs)", job_name, num_str);
  }
}

// If `nold`...`nnew` crosses `nreport` value, print update status message
// e.g. 45 -> 113, update=100 => print message
// e.g. 45 -> 313, update=100 => print message
// e.g. 45 -> 99,  update=100 => don't print message
void ctx_update2(const char *job_name, size_t nold, size_t nnew, size_t nreport)
{
  char num_str[100];
  ctx_assert2(nnew >= nold, "%zu %zu", nold, nnew);
  if((nnew - nold) >= nreport - (nold % nreport)) {
    long_to_str(nnew, num_str);
    status("[%s] Read %s entries (reads / read pairs)", job_name, num_str);
  }
}
