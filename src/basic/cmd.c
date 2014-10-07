#include "global.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "seqout.h"

#include <ctype.h>

// Global object that clones input arguments
struct Cmd {
  // These arrays are allocated
  char *merged_args;
  const char **init_argv;
  int init_argc;
  char *cmd_line_str, *cwd;
  // This is just a pointer -- not allocated
  const char *usage;
  bool printed_header; // If we have called cmd_print_status_header()
} cmd;

void cmd_init(int argc, char **argv)
{
  memset(&cmd, 0, sizeof(cmd));

  cmd.init_argc = argc;
  cmd.init_argv = cmd_clone_args(argc, argv, &cmd.merged_args);
  cmd.cmd_line_str = cmd_concat_args(argc, argv);

  cmd.cwd = ctx_malloc(PATH_MAX+1);
  if(futil_get_current_dir(cmd.cwd) == NULL) {
    warn("Couldn't get current working dir path");
    strcpy(cmd.cwd, "./");
  }
}

void cmd_destroy()
{
  ctx_free(cmd.merged_args);
  ctx_free(cmd.init_argv);
  ctx_free(cmd.cmd_line_str);
  ctx_free(cmd.cwd);
  memset(&cmd, 0, sizeof(cmd));
}

void cmd_set_usage(const char *usage) { cmd.usage = usage; }
const char*  cmd_get_usage() { return cmd.usage; }
const char*  cmd_get_cmdline() { return cmd.cmd_line_str; }
const char*  cmd_get_cwd() { return cmd.cwd; }
const char** cmd_get_argv() { return cmd.init_argv; }
int          cmd_get_argc() { return cmd.init_argc; }

// Print status updates:
// [cmd] ...
// [cwd] ...
// [version] ...
void cmd_print_status_header()
{
  // Print 1) command used 2) current working directory 3) version info
  cmd.printed_header |= (ctx_msg_out != NULL);
  status("[cmd] %s", cmd_get_cmdline());
  status("[cwd] %s", cmd_get_cwd());
  status("[version] "VERSION_STATUS_STR" k=%i..%i",
         get_min_kmer_size(), get_max_kmer_size());
}

void cmd_get_longopt_str(const struct option *longs, char shortopt,
                         char *cmdstr, size_t buflen)
{
  const char def_str[] = "-X, --Unknown";
  ctx_assert(buflen >= strlen(def_str)+1);

  string_safe_ncpy(cmdstr, def_str, buflen);
  cmdstr[1] = shortopt;

  size_t i;
  for(i = 0; longs[i].name != NULL; i++) {
    if(longs[i].val == shortopt) {
      cmdstr[0] = '-';
      cmdstr[1] = shortopt;
      string_safe_ncpy(cmdstr+6, longs[i].name, buflen-6);
      break;
    }
  }
}

void cmd_long_opts_to_short(const struct option *longs,
                            char *opts, size_t buflen)
{
  ctx_assert(buflen >= 2);
  size_t i, j;
  opts[0] = '+'; // stop as soon as we hit an arg that is not one of the options
  for(i = 0, j = 1; longs[i].name != NULL; i++) {
    if(isprint(longs[i].val)) {
      ctx_assert(j+4 <= buflen); // check we have space
      opts[j++] = longs[i].val;
      if(longs[i].has_arg) opts[j++] = ':';
      if(longs[i].has_arg == optional_argument) opts[j++] = ':';
    }
  }
  opts[j] = '\0';
}

double cmd_udouble(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  double tmp;
  if(!parse_entire_double(arg, &tmp))
    cmd_print_usage("%s requires a double >= 0: %s", cmdstr, arg);
  return tmp;
}

double cmd_udouble_nonzero(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  double tmp;
  if(!parse_entire_double(arg, &tmp) || tmp <= 0)
    cmd_print_usage("%s requires a double > 0: %s", cmdstr, arg);
  return tmp;
}

uint8_t cmd_uint8(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int 0 <= x < 255: %s", cmdstr, arg);
  return tmp;
}

int32_t cmd_int32(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  int tmp;
  if(!parse_entire_int(arg, &tmp))
    cmd_print_usage("%s requires an int: %s", cmdstr, arg);
  return tmp;
}

uint32_t cmd_uint32(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int x >= 0: %s", cmdstr, arg);
  return tmp;
}

uint32_t cmd_uint32_nonzero(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  uint32_t n = cmd_uint32(cmdstr, arg);
  if(n == 0) cmd_print_usage("%s <N> must be > 0: %s", cmdstr, arg);
  return n;
}

size_t cmd_size(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  size_t tmp;
  if(!parse_entire_size(arg, &tmp))
    cmd_print_usage("%s requires an int x >= 0: %s", cmdstr, arg);
  return tmp;
}

size_t cmd_size_nonzero(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  size_t n = cmd_size(cmdstr, arg);
  if(n == 0) cmd_print_usage("%s <N> must be > 0: %s", cmdstr, arg);
  return n;
}

size_t cmd_kmer_size(const char *cmdstr, const char *arg)
{
  size_t min_kmer_size = get_min_kmer_size();
  size_t max_kmer_size = get_max_kmer_size();
  size_t kmer_size = cmd_size_nonzero(cmdstr, arg);
  if(kmer_size < min_kmer_size || kmer_size > max_kmer_size)
    die("Please recompile with correct kmer size (%zu)", kmer_size);
  if(!(kmer_size&1)) {
    die("Invalid kmer-size (%zu): requires odd number %zu <= k <= %zu",
        kmer_size, min_kmer_size, max_kmer_size);
  }
  return kmer_size;
}

size_t cmd_parse_arg_mem(const char *cmdstr, const char *arg)
{
  ctx_assert(cmdstr != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmdstr);
  size_t mem;
  if(!mem_to_integer(arg, &mem))
    cmd_print_usage("%s %s valid options: 1024 2MB 1G", cmdstr, arg);
  return mem;
}

seq_format cmd_parse_format(const char *cmdstr, const char *arg)
{
  int tmp = seqout_str2fmt(optarg);
  if(tmp != SEQ_FMT_FASTQ && tmp != SEQ_FMT_FASTA && tmp != SEQ_FMT_PLAIN)
    cmd_print_usage("Invalid %s {FASTA,FASTQ,PLAIN} option: %s", cmdstr, arg);
  return (seq_format)tmp;
}

// Beware: we don't escape spaces etc.
// Remember to free return value
char* cmd_concat_args(int argc, char **argv)
{
  int i; size_t len = 0; char *cmdline, *str;

  for(i = 0; i < argc; i++)
    len += strlen(argv[i]) + 1;

  cmdline = ctx_malloc(len);

  for(str = cmdline, i = 0; i < argc; i++) {
    len = strlen(argv[i]);
    memcpy(str, argv[i], len);
    str[len] = ' ';
    str += len+1; // length + space
  }
  // remove last space, terminate string
  *(--str) = '\0';

  return cmdline;
}

// Remember to free pointers and strings
//    char *strmem;
//    const char **args = cmd_clone_args(argc, argv, &strmem);
//    free(strmem);
//    free(args);
const char** cmd_clone_args(int argc, char **argv, char **strmem)
{
  int i; size_t len = 0; char *str;

  for(i = 0; i < argc; i++)
    len += strlen(argv[i]) + 1;

  *strmem = ctx_malloc(len);
  const char **cmd_ptrs = ctx_malloc(argc * sizeof(char*));

  for(str = *strmem, i = 0; i < argc; i++) {
    cmd_ptrs[i] = str;
    len = strlen(argv[i]);
    memcpy(str, argv[i], len);
    str[len] = '\0';
    str += len+1; // length + null byte
  }

  return cmd_ptrs;
}


void cmd_print_usage(const char *errfmt,  ...)
{
  ctx_msg_out = stderr;

  if(!cmd.printed_header)
    cmd_print_status_header();

  pthread_mutex_lock(&ctx_biglock); // lock is never released

  if(errfmt != NULL) {
    fprintf(stderr, "\nError: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fputc('\n', stderr);
    fputc('\n', stderr);
  }

  fputs(cmd.usage, stderr);
  exit(EXIT_FAILURE);
}
