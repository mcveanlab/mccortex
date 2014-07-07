#include "global.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"

#include <ctype.h>

const char *cmd_usage = NULL;
char *cmd_line_given = NULL, *cmd_cwd = NULL;
bool printed_header = false; // If we have called cmd_print_status_header()

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

// Print status updates:
// [cmd] ...
// [cwd] ...
// [version] ...
void cmd_print_status_header()
{
  // Print 1) command used 2) current working directory 3) version info
  printed_header |= (ctx_msg_out != NULL);
  status("[cmd] %s", cmd_get_cmdline());
  status("[cwd] %s", cmd_get_cwd());
  status("[version] "VERSION_STATUS_STR"");
}

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

double cmd_udouble(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  double tmp;
  if(!parse_entire_double(arg, &tmp))
    cmd_print_usage("%s requires a double >= 0: %s", cmd, arg);
  return tmp;
}

double cmd_udouble_nonzero(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  double tmp;
  if(!parse_entire_double(arg, &tmp) || tmp <= 0)
    cmd_print_usage("%s requires a double > 0: %s", cmd, arg);
  return tmp;
}

uint8_t cmd_uint8(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int 0 <= x < 255: %s", cmd, arg);
  return tmp;
}

int32_t cmd_int32(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  int tmp;
  if(!parse_entire_int(arg, &tmp))
    cmd_print_usage("%s requires an int: %s", cmd, arg);
  return tmp;
}

uint32_t cmd_uint32(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  unsigned int tmp;
  if(!parse_entire_uint(arg, &tmp))
    cmd_print_usage("%s requires an int x >= 0: %s", cmd, arg);
  return tmp;
}

uint32_t cmd_uint32_nonzero(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
  uint32_t n = cmd_uint32(cmd, arg);
  if(n == 0) cmd_print_usage("%s <N> must be > 0: %s", cmd, arg);
  return n;
}

size_t cmd_parse_arg_mem(const char *cmd, const char *arg)
{
  ctx_assert(cmd != NULL);
  ctx_assert2(arg != NULL, "cmd: %s", cmd);
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
  ctx_msg_out = stderr;

  if(!printed_header)
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

  fputs(cmd_usage, stderr);
  exit(EXIT_FAILURE);
}
