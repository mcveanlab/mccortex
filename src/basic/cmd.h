#ifndef CMD_H_
#define CMD_H_

#include <getopt.h> // struct option

#include "seq_file/seq_file.h"

// Constants

#define CMD "mccortex"QUOTE_VALUE(MAX_KMER_SIZE)

#define DEFAULT_NTHREADS 2
#define DEFAULT_MEM 1UL<<29 /*512MB*/
#define DEFAULT_NKMERS 1UL<<22 /*4Million*/

// cmd line storing

void cmd_init(int argc, char **argv);
void cmd_destroy();

void cmd_set_usage(const char *usage);
const char*  cmd_get_usage();
const char*  cmd_get_cmdline();
const char*  cmd_get_cwd();
const char** cmd_get_argv();
int          cmd_get_argc();

// Print status updates:
// [cmd] ...
// [cwd] ...
// [version] ...
void cmd_print_status_header();

//
// General argument parsing
//

// if !x print error message
#define cmd_check(x,cmd) do {                               \
  if(!(x)) cmd_print_usage("%s given twice", cmd);          \
} while(0)

void cmd_get_longopt_str(const struct option *longs, char shortopt,
                         char *cmd, size_t buflen);
void cmd_long_opts_to_short(const struct option *longs,
                            char *opts, size_t buflen);

double cmd_udouble(const char *cmd, const char *arg);
double cmd_udouble_nonzero(const char *cmd, const char *arg);
uint8_t cmd_uint8(const char *cmd, const char *arg);
int32_t cmd_int32(const char *cmd, const char *arg);
uint32_t cmd_uint32(const char *cmd, const char *arg);
uint32_t cmd_uint32_nonzero(const char *cmd, const char *arg);
size_t cmd_size(const char *cmd, const char *arg);
size_t cmd_size_nonzero(const char *cmd, const char *arg);
size_t cmd_parse_arg_mem(const char *cmd, const char *arg);
size_t cmd_kmer_size(const char *cmdstr, const char *arg);
size_t cmd_bases(const char *cmdstr, const char *arg);// bases_to_integer

seq_format cmd_parse_format(const char *cmd, const char *arg);

// Beware: we don't escape spaces etc.
// Remember to free return value
char* cmd_concat_args(int argc, char **argv);

// Remember to free pointers and strings
//    char *strmem;
//    const char **args = cmd_clone_args(argc, argv, &strmem);
//    free(strmem);
//    free(args);
const char** cmd_clone_args(int argc, char **argv, char **strmem);

void cmd_print_usage(const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 1, 2)));

#endif
