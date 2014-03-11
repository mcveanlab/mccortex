#include "global.h"
#include <time.h>
#include <sys/time.h> // for seeding random
#include <unistd.h> // getpid
#include <pthread.h>

#include "util.h"
#include "jenkins.h" // hash functions

FILE *ctx_msg_out = NULL;
pthread_mutex_t biglock;
char cmdcode[4];

void* ctx_malloc(size_t mem, const char *file, const char *func, int line)
{
  void *ptr = malloc(mem);
  if(ptr == NULL) {
    char memstr[100];
    bytes_to_str(mem, 1, memstr);
    call_die(file, func, line, "Out of memory (malloc %s)", memstr);
  }
  return ptr;
}

void* ctx_calloc(size_t nel, size_t elsize, const char *file, const char *func, int line)
{
  void *ptr = calloc(nel, elsize);
  if(ptr == NULL) {
    char nelstr[100], elsizestr[100], memstr[100];
    ulong_to_str(nel, nelstr);
    bytes_to_str(elsize, 1, elsizestr);
    bytes_to_str(nel * elsize, 1, memstr);
    call_die(file, func, line, "Out of memory (calloc %s x %s = %s)",
             nelstr, elsizestr, memstr);
  }
  return ptr;
}

void* ctx_realloc(void *ptr, size_t mem, const char *file, const char *func, int line)
{
  void *ptr2 = realloc(ptr, mem);
  if(ptr2 == NULL) {
    char memstr[100];
    bytes_to_str(mem, 1, memstr);
    call_die(file, func, line, "Out of memory (realloc %s)", memstr);
  }
  return ptr2;
}

// Resize memory, zero new memory
void* ctx_recalloc(void *ptr, size_t oldsize, size_t newsize,
                   const char *file, const char *func, int line)
{
  ptr = ctx_realloc(ptr, newsize, file, func, line);
  memset((char*)ptr+oldsize, 0, newsize-oldsize);
  return ptr;
}

//
// Checks and asserts
//

void call_assert2(const char *file, const char *func, int line,
                  const char *asserttxt, const char *fmt, va_list argptr)
{
  pthread_mutex_lock(&biglock);
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Assert Failed %s(): %s", file, line, func, asserttxt);

  if(fmt != NULL) {
    fputs(": ", stderr);
    vfprintf(stderr, fmt, argptr);
  }

  // Print a timestamp so we know when the crash occurred
  fprintf(stderr, "\n");
  ftimestamp(stderr);
  fputs(" Assert Error\n", stderr);
  fflush(stderr);
  pthread_mutex_unlock(&biglock);
}

void call_assert_no_abort(const char *file, const char *func, int line,
                          const char *asserttxt, const char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  call_assert2(file, func, line, asserttxt, fmt, argptr);
  va_end(argptr);
}

void call_assert(const char *file, const char *func, int line,
                 const char *asserttxt, const char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  call_assert2(file, func, line, asserttxt, fmt, argptr);
  va_end(argptr);
  abort();
}

void call_die(const char *file, const char *func, int line, const char *fmt, ...)
{
  pthread_mutex_lock(&biglock);
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Error %s(): ", file, line, func);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  // Print a timestamp so we know when the crash occurred
  ftimestamp(stderr);
  fputs(" Fatal Error\n", stderr);
  exit(EXIT_FAILURE);
}

void call_warn(const char *file, const char *func, int line, const char *fmt, ...)
{
  pthread_mutex_lock(&biglock);
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Warning %s(): ", file, line, func);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  fflush(stderr);
  pthread_mutex_unlock(&biglock);
}

// A function for standard output
void fmessage(FILE *fh, const char *fmt, ...)
{
  if(fh != NULL) {
    va_list argptr;
    va_start(argptr, fmt);
    vfprintf(fh, fmt, argptr);
    va_end(argptr);
    fflush(fh);
  }
}

void ftimestamp(FILE *fh)
{
  time_t t;
  char tstr[200];

  if(fh != NULL) {
    time(&t);
    strftime(tstr, 100, "[%d %b %Y %H:%M:%S", localtime(&t));
    fprintf(fh, "%s-%s]", tstr, cmdcode);
  }
}

// A function for standard output
void fstatus(FILE *fh, const char *fmt, ...)
{
  va_list argptr;

  if(fh != NULL) {
    pthread_mutex_lock(&biglock);
    ftimestamp(fh);
    if(fmt[0] != ' ' && fmt[0] != '[') fputc(' ', fh);
    va_start(argptr, fmt);
    vfprintf(fh, fmt, argptr);
    va_end(argptr);
    if(fmt[strlen(fmt)-1] != '\n') fputc('\n', fh);
    fflush(fh);
    pthread_mutex_unlock(&biglock);
  }
}

void print_usage(const char *msg, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    pthread_mutex_lock(&biglock);
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

void seed_random()
{
  struct timeval now;
  gettimeofday(&now, NULL);

  uint32_t h;
  h = strhash_fast_mix(0, (uint32_t)now.tv_sec);
  h = strhash_fast_mix(h, (uint32_t)now.tv_usec);
  h = strhash_fast_mix(h, (uint32_t)getpid());

  srand(h);
  srand48(~h);
}

void cortex_init()
{
  // #define ALPHALEN 62
  // static const char alphabet[ALPHALEN]
  //   = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
  static const char consonants[]
    = "bcdfghjklmnpqrstvwxyzBCDFGHJKLMNPQRSTVWXYZ";
  static const char vowels[] = "aeiouAEIOU";

  seed_random();
  if(pthread_mutex_init(&biglock, NULL) != 0) {
    printf("%s:%i: mutex init failed\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // set unique code for this output
  cmdcode[0] = consonants[rand() % strlen(consonants)];
  cmdcode[1] = vowels[rand() % strlen(vowels)];
  cmdcode[2] = consonants[rand() % strlen(consonants)];
  cmdcode[3] = '\0';

  #undef ALPHALEN
}

void cortex_destroy()
{
  pthread_mutex_destroy(&biglock);
}
