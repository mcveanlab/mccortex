#ifndef CMD_H_
#define CMD_H_

#include <getopt.h> // struct option

// Constants

#define CTXCMD "ctx"QUOTE_VALUE(MAX_KMER_SIZE)
#define CMD "ctx"QUOTE_VALUE(MAX_KMER_SIZE)

#define DEFAULT_NTHREADS 2
#define DEFAULT_MEM 1UL<<29 /*512MB*/
#define DEFAULT_NKMERS 1UL<<22 /*4Million*/

// cmd line storing

void cmd_init(int argc, char **argv);
void cmd_destroy();

void cmd_set_usage(const char *usage);
const char* cmd_get_usage();
const char* cmd_get_cmdline();
const char* cmd_get_cwd();

//
// General argument parsing
//

void cmd_get_longopt_str(const struct option *longs, char shortopt,
                         char *cmd, size_t buflen);
void cmd_long_opts_to_short(const struct option *longs,
                            char *opts, size_t buflen);

double cmd_parse_arg_udouble_nonzero(const char *cmd, const char *arg);
uint8_t cmd_parse_arg_uint8(const char *cmd, const char *arg);
int32_t cmd_parse_arg_int32(const char *cmd, const char *arg);
uint32_t cmd_parse_arg_uint32(const char *cmd, const char *arg);
uint32_t cmd_parse_arg_uint32_nonzero(const char *cmd, const char *arg);
size_t cmd_parse_arg_mem(const char *cmd, const char *arg);

// Remember to free return value
char* cmd_concat_args(int argc, char **argv);

void cmd_print_usage(const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 1, 2)));

//
// Old
//

typedef struct
{
  char *cmdline;
  int cmdidx; // command specified
  bool print_help;
  // kmer, mem, ncols
  bool num_kmers_set, mem_to_use_set, num_threads_set, use_ncols_set;
  size_t num_kmers, mem_to_use, use_ncols;
  // Threads
  bool max_io_threads_set, max_work_threads_set;
  size_t max_io_threads, max_work_threads;
  // Input/output files
  bool output_file_set;
  char *output_file;
  // ctp files
  size_t num_ctp_files;
  char **ctp_files;
  // arguments not including command:
  int argc;
  char **argv;
} CmdArgs;

// Defaults
#define CMD_ARGS_INIT_MACRO { \
  .cmdline = NULL, .cmdidx = -1, .print_help = false, \
  .num_kmers_set = false, .num_kmers = DEFAULT_NKMERS, \
  .mem_to_use_set = false, .mem_to_use = DEFAULT_MEM, \
  .max_io_threads_set = false, .max_io_threads = 4, \
  .max_work_threads_set = false, .max_work_threads = DEFAULT_NTHREADS, \
  .use_ncols_set = false, .use_ncols = 1, \
  .output_file_set = false, .output_file = NULL, \
  .num_ctp_files = 0, .ctp_files = NULL, \
  .argc = 0, .argv = NULL}

void cmd_alloc(CmdArgs *args, int argc, char **argv);
void cmd_free(CmdArgs *args);

// accptopts is a string of valid args,
// e.g. "tk" accepts kmer-size and number of threads
// NULL means anything valid, "" means no args valid
void cmd_accept_options(const CmdArgs *args, const char *accptopts,
                        const char *usage);
void cmd_require_options(const CmdArgs *args, const char *requireopts,
                         const char *usage);

#endif
