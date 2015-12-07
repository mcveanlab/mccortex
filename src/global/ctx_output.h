#ifndef CTX_OUTPUT_H_
#define CTX_OUTPUT_H_

//
// Exit with an error, give a warning: die() / warn()
//

// Output destination and lock for output
// set to NULL to turn off message printing
extern FILE *ctx_msg_out;
extern pthread_mutex_t ctx_biglock;
extern char ctx_cmdcode[4];

void dief(const char *file, const char *func, int line, const char *fmt, ...)
__attribute__((format(printf, 4, 5)))
__attribute__((noreturn));

void warnf(const char *file, const char *func, int line, const char *fmt, ...)
__attribute__((format(printf, 4, 5)));

void messagef(FILE *fh, const char *fmt, ...)
__attribute__((format(printf, 2, 3)));

void timestampf(FILE *fh);

void statusf(FILE *fh, const char *fmt, ...)
__attribute__((format(printf, 2, 3)));

#define die(fmt, ...) dief(__FILE__, __func__, __LINE__, fmt, ##__VA_ARGS__)
#define warn(fmt, ...) warnf(__FILE__, __func__, __LINE__, fmt, ##__VA_ARGS__)

#define message(fmt,...) messagef(ctx_msg_out,fmt, ##__VA_ARGS__)
#define timestamp()      timestampf(ctx_msg_out)
#define status(fmt,...)  statusf(ctx_msg_out,fmt, ##__VA_ARGS__)

void print_usage(const char *msg, const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 2, 3)));

void ctx_output_init();
void ctx_output_destroy();

// Print progress every 5M reads
#define CTX_UPDATE_REPORT_RATE 5000000

void ctx_update(const char *job_name, size_t niter);

// If `nold`...`nnew` crosses `nreport` value, print update status message
// e.g. 45 -> 113, update=100 => print message
// e.g. 45 -> 313, update=100 => print message
// e.g. 45 -> 99,  update=100 => don't print message
void ctx_update2(const char *job_name, size_t nold, size_t nnew, size_t nreport);

#endif /* CTX_OUTPUT_H_ */
