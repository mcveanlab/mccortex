#include "global.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"

#include <ctype.h>

const char *cmd_usage = NULL;
char *cmd_line_given = NULL, *cmd_cwd = NULL;

void cmd_init(int argc, char **argv)
{
  cmd_line_given = cmd_concat_args(argc, argv);
  cmd_cwd = ctx_malloc(PATH_MAX+1);
  if(futil_get_current_dir(cmd_cwd) == NULL) {
    warn("Couldn't get current working dir path");
    strcpy(cmd_cwd, "./");
  }
}

void cmd_destroy()
{
  ctx_free(cmd_line_given);
  ctx_free(cmd_cwd);
  cmd_usage = cmd_line_given = cmd_cwd = NULL;
}

void cmd_set_usage(const char *usage) { cmd_usage = usage; }
const char* cmd_get_usage() { return cmd_usage; }
const char* cmd_get_cmdline() { return cmd_line_given; }
const char* cmd_get_cwd() { return cmd_cwd; }

void cmd_get_longopt_str(const struct option *longs, char shortopt,
                         char *cmd, size_t buflen)
{
  const char def_str[] = "-X, --Unknown";
  ctx_assert(buflen >= strlen(def_str)+1);

  string_safe_ncpy(cmd, def_str, buflen);
  cmd[1] = shortopt;

  size_t i;
  for(i = 0; longs[i].name != NULL; i++) {
    if(longs[i].val == shortopt) {
      cmd[0] = '-';
      cmd[1] = shortopt;
      string_safe_ncpy(cmd+6, longs[i].name, buflen-6);
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

double cmd_udouble_nonzero(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  double tmp;
  if(!parse_entire_double(arg, &tmp) || tmp <= 0)
    cmd_print_usage("%s requires a double > 0: %s", cmd, arg);
  return tmp;
}

uint8_t cmd_uint8(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int 0 <= x < 255: %s", cmd, arg);
  return tmp;
}

int32_t cmd_int32(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  int tmp;
  if(!parse_entire_int(arg, &tmp))
    cmd_print_usage("%s requires an int: %s", cmd, arg);
  return tmp;
}

uint32_t cmd_uint32(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int x >= 0: %s", cmd, arg);
  return tmp;
}

uint32_t cmd_uint32_nonzero(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  uint32_t n = cmd_uint32(cmd, arg);
  if(n == 0) cmd_print_usage("%s <N> must be > 0: %s", cmd, arg);
  return n;
}

size_t cmd_parse_arg_mem(const char *cmd, const char *arg)
{
  ctx_assert(arg != NULL);
  size_t mem;
  if(!mem_to_integer(arg, &mem))
    cmd_print_usage("%s %s valid options: 1024 2MB 1G", cmd, arg);
  return mem;
}

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

void cmd_print_usage(const char *errfmt,  ...)
{
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

  fputs(cmd_usage, stderr);
  exit(EXIT_FAILURE);
}


//
// Old
//

void cmd_accept_options(const CmdArgs *args, const char *accptopts,
                        const char *usage)
{
  ctx_assert(accptopts != NULL);
  ctx_assert(usage != NULL);
  if(args->print_help) print_usage(usage, NULL);
  if(args->num_kmers_set && strchr(accptopts,'n') == NULL)
    print_usage(usage, "-n <hash-entries> argument not valid for this command");
  if(args->mem_to_use_set && strchr(accptopts,'m') == NULL)
    print_usage(usage, "-m <memory> argument not valid for this command");
  if(args->max_io_threads_set && strchr(accptopts,'a') == NULL)
    print_usage(usage, "-a <iothreads> argument not valid for this command");
  if(args->max_work_threads_set && strchr(accptopts,'t') == NULL)
    print_usage(usage, "-t <threads> argument not valid for this command");
  if(args->use_ncols_set && strchr(accptopts,'c') == NULL)
    print_usage(usage, "-s <ncols> argument not valid for this command");
  if(args->output_file_set && strchr(accptopts,'o') == NULL)
    print_usage(usage, "-o <file> argument not valid for this command");
  if(args->num_ctp_files > 0 && strchr(accptopts,'p') == NULL)
    print_usage(usage, "-p <in.ctp> argument not valid for this command");
  // Check for programmer error - all options should be valid
  const char *p;
  for(p = accptopts; *p != '\0'; p++) {
    if(strchr("nmkatcfop", *p) == NULL)
      die("Invalid option: %c", *p);
  }
}

void cmd_require_options(const CmdArgs *args, const char *requireopts,
                         const char *usage)
{
  if(requireopts == NULL) return;
  if(args->argc == 0 && *requireopts != '\0') print_usage(usage, NULL);
  for(; *requireopts != '\0'; requireopts++)
  {
    if(*requireopts == 'n') {
      if(!args->num_kmers_set)
        die("-n <hash-entries> argument required for this command");
    }
    else if(*requireopts == 'm') {
      if(!args->mem_to_use_set)
        die("-m <memory> argument required for this command");
    }
    else if(*requireopts == 'a') {
      if(!args->max_io_threads_set)
        die("-a <iothreads> argument required for this command");
    }
    else if(*requireopts == 't') {
      if(!args->max_work_threads_set)
        die("-t <threads> argument required for this command");
    }
    else if(*requireopts == 's') {
      if(!args->use_ncols_set)
        die("-s <ncols> argument required for this command");
    }
    else if(*requireopts == 'o') {
      if(!args->output_file_set)
        die("-o <file> argument required for this command");
    }
    else if(*requireopts == 'p') {
      if(args->num_ctp_files == 0)
        die("-p <in.ctp> argument required for this command");
    }
    else warn("Ignored required option: %c", *requireopts);
  }
}

void cmd_alloc(CmdArgs *args, int argc, char **argv)
{
  CmdArgs tmp = CMD_ARGS_INIT_MACRO;
  memcpy(args, &tmp, sizeof(CmdArgs));

  args->ctp_files = ctx_malloc((size_t)(argc/2) * sizeof(char*));
  args->num_ctp_files = 0;

  // Get command index
  bool is_ctx_cmd = (strstr(argv[0],"ctx") != NULL);

  args->argc = 0;
  args->argv = ctx_malloc((size_t)argc * sizeof(char*));

  // Get cmdline string
  size_t len = (size_t)argc - 1; // spaces
  int i;
  for(i = 0; i < argc; i++) len += strlen(argv[i]);
  args->cmdline = ctx_malloc((len+1) * sizeof(char));
  sprintf(args->cmdline, "%s", argv[0]);
  for(i = 1, len = strlen(argv[0]); i < argc; len += strlen(argv[i]) + 1, i++)
    sprintf(args->cmdline+len, " %s", argv[i]);

  // argv[0] = ctx, argv[1] = cmd, argv[2..] = args
  for(i = is_ctx_cmd ? 2 : 1; i < argc; i++)
  {
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--h") || !strcmp(argv[i], "--help"))
    {
      args->print_help = true;
    }
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nkmers") == 0)
    {
      if(i + 1 == argc) die("%s <hash-kmers> requires an argument", argv[i]);
      if(args->num_kmers_set) die("-n <hash-kmers> given more than once");
      if(!mem_to_integer(argv[i+1], &args->num_kmers) || args->num_kmers == 0)
        die("Invalid hash size: %s", argv[i+1]);
      args->num_kmers_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--memory") == 0)
    {
      if(i + 1 == argc) die("%s <memory> requires an argument", argv[i]);
      if(args->mem_to_use_set) die("-m <memory> given more than once");
      if(!mem_to_integer(argv[i+1], &args->mem_to_use) || args->mem_to_use == 0)
        die("Invalid memory argument: %s", argv[i+1]);
      args->mem_to_use_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--asyncio") == 0)
    {
      if(i + 1 == argc) die("%s <iothreads> requires an argument", argv[i]);
      if(args->max_io_threads_set) die("-a <iothreads> given more than once");
      if(!parse_entire_size(argv[i+1], &args->max_io_threads) || !args->max_io_threads)
        die("Invalid number of io threads: %s", argv[i+1]);
      args->max_io_threads_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--threads") == 0)
    {
      if(i + 1 == argc) die("%s <threads> requires an argument", argv[i]);
      if(args->max_work_threads_set) die("-t <threads> given more than once");
      if(!parse_entire_size(argv[i+1], &args->max_work_threads) || !args->max_work_threads)
        die("Invalid number of worker threads: %s", argv[i+1]);
      args->max_work_threads_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--ncols") == 0)
    {
      if(i + 1 == argc) die("%s <ncols> requires an argument", argv[i]);
      if(args->use_ncols_set) die("-s <ncols> given more than once");
      if(!parse_entire_size(argv[i+1], &args->use_ncols) || args->use_ncols == 0)
        die("Invalid number of colours to use: %s", argv[i+1]);
      args->use_ncols_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--out") == 0)
    {
      if(i + 1 == argc) die("%s <file> requires an argument", argv[i]);
      if(args->output_file_set) die("%s <file> given more than once", argv[i]);
      args->output_file = argv[i+1];
      args->output_file_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--path") == 0 ||
            strcmp(argv[i], "--paths") == 0)
    {
      if(i + 1 == argc) die("%s <in.ctp> requires an argument", argv[i]);
      args->ctp_files[args->num_ctp_files++] = argv[i+1];
      i++;
    }
    else {
      args->argv[args->argc++] = argv[i];
    }
  }
}

void cmd_free(CmdArgs *args)
{
  ctx_free(args->ctp_files);
  ctx_free(args->argv);
  ctx_free(args->cmdline);
}

