#include "global.h"
#include "util.h"
#include "cmd.h"

static const char usage[] =
"\n"
"Commands:\n" // priority show
"  build     [FASTA/FASTQ/BAM -> binary]\n"
"  clean     [clean population / sample against merged]  (unfinished)\n" // 4
"  join      [combine binaries]\n"
"  intersect [dump intersection of a.ctx with b.ctx]\n"
"  subgraph  [filter a subgraph]\n"
"  reads     [filter reads]\n"
"  extend    [extend contigs using a population graph]\n"
"  contigs   [pull out contigs for a sample]             (unfinished)\n" // 2
"  thread    [thread reads through cleaned population]\n"
"  pview     [view read threading information]\n"
"  pmerge    [merge path files (.ctp)]                   (unfinished)\n" // 5
"  call      [call variants]\n"
"  diverge   [path divergence caller]                    (unfinished)\n" // 3
"  unique    [remove duplicated bubbles -> VCF]\n"
"  covg      [add covg to a VCF file]                    (unfinished)\n" // 1
"  place     [place variants and genotype]\n"
"\n"
"  Type a command with no arguments to see help.\n"
"\n"
"Common Options:\n"
"  -m <M>    Memory e.g. 1GB [default: 1GB]\n"
"  -h <H>    Hash entries [default: 4M (~4 million)]\n"
"  -g <G>    Species genome size [default: 3.1Gbp]\n"
"  -t <T>    Number of threads [default: 2]\n"
"  -k <K>    Kmer size [default: read from binaries]\n"
"\n";

const int num_cmds = 16;

const char *cmds[num_cmds]
  = {"build", "clean", "join", "intersect", "subgraph", "reads", "extend",
     "contigs", "thread", "pview", "pmerge", "call", "diverge", "unique",
     "covg", "place"};

int (*funcs[num_cmds])(CmdArgs *cmd_args)
  = {ctx_build, ctx_clean, ctx_join, ctx_intersect, ctx_subgraph, ctx_reads,
     ctx_extend, ctx_contigs, ctx_thread, ctp_view, ctp_merge, ctx_call,
     ctx_diverge, ctx_unique, ctx_covg, ctx_place};

// Commands not implemented yet
int ctp_merge(CmdArgs *args) { die("Not implemented yet [cmd: %s]", args->cmdline); }

static void cmd_init(CmdArgs *args, int argc, char **argv)
{
  args->genome_size = 3100000000;
  args->num_kmers = 4UL << 20; // 4M
  args->mem_to_use = 1UL << 30; // 1G
  args->num_threads = 2;

  args->genome_size_set = false;
  args->mem_to_use_set = false;
  args->mem_to_use_set = false;
  args->num_threads_set = false;

  args->argc = 0;
  args->argv = malloc(argc * sizeof(char**));

  size_t len = argc - 1; // spaces
  int i;
  for(i = 0; i < argc; i++) len += strlen(argv[i]);
  args->cmdline = malloc((len+1) * sizeof(char));
  sprintf(args->cmdline, "%s", argv[0]);
  for(i = 1, len = strlen(argv[0]); i < argc; len += strlen(argv[i]) + 1, i++)
    sprintf(args->cmdline+len, " %s", argv[i]);

  // argv[0] = ctx, argv[1] = cmd, argv[2..] = args
  for(i = 2; i < argc; i++)
  {
    if(strcmp(argv[i], "-m") == 0)
    {
      if(i + 1 == argc) die("-m requires an argument");
      if(!mem_to_integer(argv[i+1], &args->mem_to_use) || args->mem_to_use == 0)
        die("Invalid memory argument: %s", argv[i+1]);
      args->mem_to_use_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0)
    {
      if(i + 1 == argc) die("-g requires an argument");
      if(!bases_to_integer(argv[i+1], &args->genome_size) || args->genome_size == 0)
        die("Invalid genome size argument: %s", argv[i+1]);
      args->genome_size_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0)
    {
      if(i + 1 == argc) die("-t requires an argument");
      if(!parse_entire_uint(argv[i+1], &args->num_threads) || args->num_threads == 0)
        die("Invalid number of threads: %s", argv[i+1]);
      args->num_threads_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-h") == 0)
    {
      if(i + 1 == argc) die("-h requires an argument");
      if(!mem_to_integer(argv[i+1], &args->num_kmers) || args->num_kmers == 0)
        die("Invalid hash size: %s", argv[i+1]);
      args->num_kmers_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-k") == 0)
    {
      if(i + 1 == argc) die("-k requires an argument");
      if(!parse_entire_uint(argv[i+1], &args->kmer_size) ||
         args->kmer_size < 3 || !(args->kmer_size & 0x1))
        die("kmer size (-k) must be an odd int >= 3: %s", argv[i+1]);
      if(args->kmer_size < MIN_KMER_SIZE || args->kmer_size > MAX_KMER_SIZE)
        die("Please recompile with correct kmer size (kmer_size: %u)", args->kmer_size);
      args->kmer_size_set = true;
      i++;
    }
    else {
      args->argv[args->argc++] = argv[i];
    }
  }
}

static void cmd_free(CmdArgs *args)
{
  free(args->argv);
  free(args->cmdline);
}

int main(int argc, char **argv)
{
  if(argc == 1) {
    print_usage(usage, NULL);
    return EXIT_FAILURE;
  }

  int cmd;
  for(cmd = 0; cmd < num_cmds && strcmp(cmds[cmd],argv[1]) != 0; cmd++);
  if(cmd == num_cmds) print_usage(usage, "Unrecognised command: %s", argv[1]);

  int ret;
  CmdArgs args;
  cmd_init(&args, argc, argv);
  ret = funcs[cmd](&args);
  cmd_free(&args);
  return ret;
}
