#include "global.h"
#include <time.h>

#include "commands.h"
#include "util.h"
#include "file_util.h"

// To add a new command to ctx31 <cmd>:
// 0. create a file src/commands/ctx_X.c
// 1. add function and usage declaration to src/commands/commands.h
// 2. add entry to cmdobjs below
// 3. that's it!  Write the function I guess...

typedef struct
{
  const char *cmd, *blurb, *usage, *optargs, *reqargs;
  int minargs, maxargs; // counts AFTER standard args taken
  int hide; // set hide to >0 to remove from listings
  int (*func)(int argc, char **argv);
} CtxCmd;

CtxCmd cmdobjs[] = {
{
  .cmd = "build", .func = ctx_build, .hide = false,
  .blurb = "construct cortex graph from FASTA/FASTQ/BAM",
  .usage = build_usage
},
{
  .cmd = "view", .func = ctx_view, .hide = false,
  .blurb = "text view of a cortex graph file (.ctx)",
  .usage = view_usage
},
{
  .cmd = "check", .func = ctx_health_check, .hide = false,
  .blurb = "load and check graph (.ctx) and path (.ctp) files",
  .usage = health_usage
},
{
  .cmd = "clean", .func = ctx_clean, .hide = false,
  .blurb = "clean errors from a graph",
  .usage = clean_usage
},
{
  .cmd = "join", .func = ctx_join, .hide = false,
  .blurb = "combine graphs, filter graph intersections",
  .usage = join_usage
},
{
  .cmd = "supernodes", .func = ctx_supernodes, .hide = false,
  .blurb = "pull out supernodes",
  .usage = supernodes_usage
},
{
  .cmd = "subgraph", .func = ctx_subgraph, .hide = false,
  .blurb = "filter a subgraph using seed kmers",
  .usage = subgraph_usage
},
{
  .cmd = "reads", .func = ctx_reads, .hide = false,
  .blurb = "filter reads against a graph",
  .usage = reads_usage
},
{
  .cmd = "contigs", .func = ctx_contigs, .hide = false,
  .blurb = "pull out contigs for a sample",
  .usage = contigs_usage
},
{
  .cmd = "inferedges", .func = ctx_infer_edges, .hide = false,
  .blurb = "infer graph edges between kmers before calling `thread`",
  .usage = inferedges_usage
},
{
  .cmd = "thread", .func = ctx_thread, .hide = false,
  .blurb = "thread reads through cleaned graph",
  .usage = thread_usage,
},
{
  .cmd = "correct", .func = ctx_correct, .hide = false,
  .blurb = "error correct reads",
  .usage = correct_usage
},
{
  .cmd = "pjoin", .func = ctx_pjoin, .hide = false,
  .blurb = "merge path files (.ctp)",
  .usage = pjoin_usage
},
{
  .cmd = "bubbles", .func = ctx_bubbles, .hide = false,
  .blurb = "find bubbles in graph which are potential variants",
  .usage = bubbles_usage
},
{
  .cmd = "breakpoints", .func = ctx_breakpoints, .hide = false,
  .blurb = "use a trusted assembled genome to call large events",
  .usage = breakpoints_usage
},
{
  .cmd = "unique", .func = ctx_unique, .hide = true,
  .blurb = "remove duplicated bubbles, produce VCF",
  .usage = unique_usage
},
{
  .cmd = "place", .func = ctx_place, .hide = true,
  .blurb = "place variants against a reference",
  .usage = place_usage
},
{
  .cmd = "coverage", .func = ctx_coverage, .hide = false,
  .blurb = "print contig coverage",
  .usage = coverage_usage
},
{
  .cmd = "rmsubstr", .func = ctx_rmsubstr, .hide = false,
  .blurb = "reduce set of strings to remove substrings",
  .usage = rmsubstr_usage
},
{
  .cmd = "calls2vcf", .func = ctx_calls2vcf, .hide = true,
  .blurb = "reduce set of strings to remove substrings",
  .usage = calls2vcf_usage
}
};


//
// Command listing
//

static const char header[] =
"\n"
"usage: "CMD" <command> [options] <args>\n"
"version: "VERSION_STATUS_STR"\n"
"\n";

static const char options[] =
"  Type a command with no arguments to see help.\n"
"\n"
"Common Options:\n"
"  -m, --memory <M>      Memory e.g. 1GB [default: 1GB]\n"
"  -n, --nkmers <H>      Hash entries [default: 4M, ~4 million]\n"
"  -t, --threads <T>     Limit on proccessing threads [default: 2]\n"
"  -o, --out <file>      Output file\n"
"  -p, --paths <in.ctp>  Assembly file to load (can specify multiple times)\n"
"\n";

static int ctxcmd_cmp(const void *aa, const void *bb)
{
  const CtxCmd *a = (const CtxCmd*)aa, *b = (const CtxCmd*)bb;
  return strcmp(a->cmd, b->cmd);
}

static void print_help(FILE *out, const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 2, 3)));

static void print_help(FILE *out, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(out, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(out, errfmt, argptr);
    va_end(argptr);
  }

  size_t i, j, s, maxlen = 0, n = sizeof(cmdobjs) / sizeof(CtxCmd);

  qsort(cmdobjs, n, sizeof(CtxCmd), ctxcmd_cmp);

  for(i = 0; i < n; i++)
    if(!cmdobjs[i].hide) maxlen = MAX2(maxlen, strlen(cmdobjs[i].cmd));

  fprintf(out, "%s", header);

  fprintf(out, "Commands:   ");
  for(i = 0, j = 0; i < n; i++) {
    if(!cmdobjs[i].hide) {
      fprintf(out, "%s%s", j > 0 ? "            " : "", cmdobjs[i].cmd);
      s = maxlen - strlen(cmdobjs[i].cmd) + 2;
      while(--s != SIZE_MAX) fputc(' ', out);
      fprintf(out, "%s\n", cmdobjs[i].blurb);
      j++;
    }
  }

  fprintf(out, "\n%s", options);
  exit(EXIT_FAILURE);
}

static const CtxCmd* ctx_get_command(const char* cmd)
{
  size_t i, n = sizeof(cmdobjs) / sizeof(CtxCmd);
  for(i = 0; i < n; i++)
    if(strcasecmp(cmdobjs[i].cmd,cmd) == 0)
      return &cmdobjs[i];
  return NULL;
}

// Print status updates:
// [cmd] ...
// [cwd] ...
// [version] ...
static void print_status_header()
{
  // Print 1) command used 2) current working directory 3) version info
  status("[cmd] %s", cmd_get_cmdline());
  status("[cwd] %s", cmd_get_cwd());
  status("[version] "VERSION_STATUS_STR"");
}

int main(int argc, char **argv)
{
  time_t start, end;
  time(&start);

  ctx_msg_out = stderr;
  cortex_init();
  cmd_init(argc, argv);

  if(argc == 1) print_help(stderr, NULL);
  const CtxCmd *cmd = ctx_get_command(argv[1]);
  if(cmd == NULL) print_help(stderr, "Unrecognised command: %s", argv[1]);

  // Once we have set cmd_usage, we can call cmd_print_usage() from anywhere
  cmd_set_usage(cmd->usage);

  // If no arguments after command, print help
  if(argc == 2) cmd_print_usage(NULL);

  // Print status header
  print_status_header();

  SWAP(argv[1],argv[0]);
  int ret = cmd->func(argc-1, argv+1);

  time(&end);
  cmd_destroy();

  // Warn if more allocations than deallocations
  size_t still_alloced = alloc_get_num_allocs() - alloc_get_num_frees();
  if(still_alloced) warn("%zu allocates not free'd.", still_alloced);

  status(ret == 0 ? "Done." : "Fail.");

  // Print time taken
  double diff = difftime(end,start);
  if(diff < 60) status("[time] %.2lf seconds\n", diff);
  else {
    char timestr[100];
    seconds_to_str((size_t)diff, timestr);
    status("[time] %.2lf seconds (%s)\n", diff, timestr);
  }

  cortex_destroy();

  return ret;
}
