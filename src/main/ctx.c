#include "global.h"
#include <time.h>
#include "cmd.h"
#include "util.h"

#include <execinfo.h>
#include <signal.h>

static const char usage[] =
"\n"
"usage: "CMD" <command> [options] <args>\n"
"version: "CTXVERSIONSTR"; zlib: "ZLIB_VERSION"\n"
"\n"
"Command:  build       FASTA/FASTQ/BAM -> binary graph file\n"
"          view        view and check a cortex graph file (.ctx)\n"
"          healthcheck load and check a cortex graph file (.ctx)\n"
"          clean       clean errors from a graph\n"
"          join        combine graphs, filter graph intersections\n"
"          supernodes  pull out supernodes\n"
"          subgraph    filter a subgraph\n"
"          reads       filter reads against a graph\n"
"          extend      extend contigs using a graph\n"
"          contigs     pull out contigs for a sample\n"
"          inferedges  infer edges before calling `thread`\n"
"          thread      thread reads through cleaned population\n"
"          pview       view read threading information\n"
"          pmerge      merge path files (.ctp)\n"
"          call        call variants\n"
"          unique      remove duplicated bubbles, produce VCF\n"
"          place       place variants and genotype\n"
// "          diverge     path divergence caller               (unfinished)\n"
// "          covg        add covg to a VCF file               (unfinished)\n"
"\n"
"  Type a command with no arguments to see help.\n"
"\n"
"Common Options:\n"
"  -m <M>       Memory e.g. 1GB [default: 1GB]\n"
"  -h <H>       Hash entries [default: 4M (~4 million)]\n"
"  -g <G>       Species genome size [default: 3.1Gbp]\n"
"  -t <T>       Number of threads [default: 2]\n"
"  -k <K>       Kmer size [default: read from binaries]\n"
"  -f <file>    Input file\n"
"  -p <in.ctp>  Assembly file\n"
"\n";

int main(int argc, char **argv)
{
  signal(SIGSEGV, errhandler);   // install our handler

  CmdArgs args;
  time_t start, end;

  time(&start);
  if(argc == 1) print_usage(usage, NULL);
  cmd_alloc(&args, argc, argv);
  if(args.cmdidx == -1) print_usage(usage, "Unrecognised command: %s", argv[1]);

  boolean usestderr = (strcmp(argv[1], "view") == 0 ||
                       strcmp(argv[1], "pview") == 0 ||
                       strcmp(argv[1], "place") == 0 ||
                       strcmp(argv[1], "contigs") == 0);

  ctx_msg_out = usestderr ? stderr : stdout;

  status("[cmd] %s\n", args.cmdline);
  int ret = ctx_funcs[args.cmdidx](&args);
  cmd_free(&args);
  time(&end);

  status(ret == 0 ? "Done." : "Fail.");

  if(strcmp(argv[1],"view") != 0)
  {
    // Print time taken
    double diff = difftime (end,start);
    if(diff < 60) status("[time] %.2lf seconds\n", diff);
    else {
      char timestr[100];
      seconds_to_str(diff, timestr);
      status("[time] %.2lf seconds (%s)\n", diff, timestr);
    }
  }

  return ret;
}
