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
  int (*func)(CmdArgs *cmd_args);
} CtxCmd;

CtxCmd cmdobjs[] = {
{
  .cmd = "build", .func = ctx_build, .hide = 0,
  .minargs = 3, .maxargs = INT_MAX, .optargs = "atmnk", .reqargs = "k",
  .blurb = "construct cortex graph from FASTA/FASTQ/BAM",
  .usage = build_usage
},
{
  .cmd = "view", .func = ctx_view, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "", .reqargs = "",
  .blurb = "view and check a cortex graph file (.ctx)",
  .usage = view_usage
},
{
  .cmd = "check", .func = ctx_health_check, .hide = 0,
  .minargs = 1, .maxargs = 2, .optargs = "pmn", .reqargs = "",
  .blurb = "load and check graph (.ctx) and path (.ctp) files",
  .usage = health_usage
},
{
  .cmd = "clean", .func = ctx_clean, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "mnco", .reqargs = "o",
  .blurb = "clean errors from a graph",
  .usage = clean_usage
},
{
  .cmd = "join", .func = ctx_join, .hide = 0,
  .minargs = 2, .maxargs = INT_MAX, .optargs = "mnc", .reqargs = "",
  .blurb = "combine graphs, filter graph intersections",
  .usage = join_usage
},
{
  .cmd = "supernodes", .func = ctx_supernodes, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "mnpo", .reqargs = "",
  .blurb = "pull out supernodes",
  .usage = supernodes_usage
},
{
  .cmd = "subgraph", .func = ctx_subgraph, .hide = 0,
  .minargs = 4, .maxargs = INT_MAX, .optargs = "mnco", .reqargs = "",
  .blurb = "filter a subgraph using seed kmers",
  .usage = subgraph_usage
},
{
  .cmd = "reads", .func = ctx_reads, .hide = 0,
  .minargs = 4, .maxargs = INT_MAX, .optargs = "mn", .reqargs = "",
  .blurb = "filter reads against a graph",
  .usage = reads_usage
},
{
  .cmd = "extend", .func = ctx_extend, .hide = 1,
  .minargs = 4, .maxargs = 4, .optargs = "mn", .reqargs = "",
  .blurb = "extend contigs using a graph",
  .usage = extend_usage
},
{
  .cmd = "contigs", .func = ctx_contigs, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "mnpo", .reqargs = "",
  .blurb = "pull out contigs for a sample",
  .usage = contigs_usage
},
{
  .cmd = "inferedges", .func = ctx_infer_edges, .hide = 0,
  .minargs = 1, .maxargs = 3, .optargs = "mno", .reqargs = "",
  .blurb = "infer graph edges before calling `thread`",
  .usage = inferedges_usage
},
{
  .cmd = "thread", .func = ctx_thread, .hide = 0,
  .minargs = 2, .maxargs = INT_MAX, .optargs = "atmnpo", .reqargs = "o",
  .blurb = "thread reads through cleaned graph",
  .usage = thread_usage,
},
{
  .cmd = "correct", .func = ctx_correct, .hide = 0,
  .minargs = 2, .maxargs = INT_MAX, .optargs = "atmnp", .reqargs = "",
  .blurb = "error correct reads",
  .usage = correct_usage
},
{
  .cmd = "pview", .func = ctx_pview, .hide = 0,
  .minargs = 1, .maxargs = 3, .optargs = "mn", .reqargs = "",
  .blurb = "view read threading information",
  .usage = pview_usage
},
{
  .cmd = "pjoin", .func = ctx_pjoin, .hide = 0,
  .minargs = 3, .maxargs = INT_MAX, .optargs = "mnf", .reqargs = "",
  .blurb = "merge path files (.ctp)",
  .usage = pjoin_usage
},
{
  .cmd = "call", .func = ctx_call, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "tmnpo", .reqargs = "o",
  .blurb = "call variants with bubble caller",
  .usage = call_usage
},
{
  .cmd = "unique", .func = ctx_unique, .hide = 0,
  .minargs = 2, .maxargs = 2, .optargs = "", .reqargs = "",
  .blurb = "remove duplicated bubbles, produce VCF",
  .usage = unique_usage
},
{
  .cmd = "place", .func = ctx_place, .hide = 0,
  .minargs = 3, .maxargs = INT_MAX, .optargs = "o", .reqargs = "",
  .blurb = "place variants and genotype",
  .usage = place_usage
},
{
  .cmd = "breakpoints", .func = ctx_breakpoints, .hide = 0,
  .minargs = 1, .maxargs = INT_MAX, .optargs = "opmn", .reqargs = "o",
  .blurb = "Use trusted assembled genome to call large events",
  .usage = breakpoints_usage
},
{
  .cmd = "coverage", .func = ctx_coverage, .hide = 0,
  .minargs = 3, .maxargs = INT_MAX, .optargs = "omn", .reqargs = "",
  .blurb = "Get contig coverage",
  .usage = coverage_usage
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
"  -m --memory <M>      Memory e.g. 1GB [default: 1GB]\n"
"  -n --nkmers <H>      Hash entries [default: 4M, ~4 million]\n"
"  -c --ncols <C>       Number of graph colours to load at once [default: 1]\n"
"  -a --asyncio <A>     Limit on file reading threads [default: 4]\n"
"  -t --threads <T>     Limit on proccessing threads [default: 2]\n"
"  -k --kmer <K>        Kmer size [default: read from graph files]\n"
"  -f --file <file>     Input file\n"
"  -o --out <file>      Output file\n"
"  -p --paths <in.ctp>  Assembly file to load (can specify multiple times)\n"
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

static void print_header(const char *cmdline)
{
  char abspath[PATH_MAX+1];
  status("[cmd] %s", cmdline);
  if(futil_get_current_dir(abspath) != NULL) status("[cwd] %s", abspath);
  status("[version] "VERSION_STATUS_STR"");
}

int main(int argc, char **argv)
{
  CmdArgs args;
  time_t start, end;
  time(&start);

  cortex_init();
  ctx_msg_out = stderr;

  if(argc == 1) print_help(stderr, NULL);

  cmd_alloc(&args, argc, argv);

  const CtxCmd *cmd = ctx_get_command(argv[1]);
  if(cmd == NULL) print_help(stderr, "Unrecognised command: %s", argv[1]);

  // Once we have set cmd_usage, we can call cmd_print_usage() from anywhere
  cmd_usage = cmd->usage;

  // If no arguments after command, print help
  if(argc == 2) cmd_print_usage(NULL);

  print_header(args.cmdline);

  // Check number of args, required args, optional args
  cmd_accept_options(&args, cmd->optargs, cmd->usage);
  cmd_require_options(&args, cmd->reqargs, cmd->usage);

  if(args.argc < cmd->minargs) cmd_print_usage("Too few arguments");
  if(args.argc > cmd->maxargs) cmd_print_usage("Too many arguments");

  // Run command
  int ret = cmd->func(&args);

  cmd_free(&args);
  time(&end);

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
