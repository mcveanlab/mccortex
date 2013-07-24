#include "global.h"
#include "util.h"
#include "cmd.h"

const char *cmds[NUM_CMDS]
  = {"build", "clean", "join", "intersect", "subgraph", "reads", "extend",
     "contigs", "thread", "pview", "pmerge", "call", "diverge", "unique",
     "covg", "place"};

int (*ctx_funcs[NUM_CMDS])(CmdArgs *cmd_args)
  = {ctx_build, ctx_clean, ctx_join, ctx_intersect, ctx_subgraph, ctx_reads,
     ctx_extend, ctx_contigs, ctx_thread, ctx_pview, ctx_pmerge, ctx_call,
     ctx_diverge, ctx_unique, ctx_covg, ctx_place};

// static int ctx_notimpl(CmdArgs *args) {
//   warn("Command not implemented [cmd: %s]", args->cmdline); return EXIT_FAILURE;
// }

void cmd_alloc(CmdArgs *args, int argc, char **argv)
{
  args->genome_size = 3100000000;
  args->num_kmers = 4UL << 20; // 4M
  args->mem_to_use = 1UL << 30; // 1G
  args->num_threads = 2;

  args->genome_size_set = false;
  args->mem_to_use_set = false;
  args->mem_to_use_set = false;
  args->num_threads_set = false;

  // Get command index
  boolean is_ctx_cmd = (strstr(argv[0],"ctx") != NULL);
  if(argc >= 2 && is_ctx_cmd) {
    int cmd;
    for(cmd = 0; cmd < NUM_CMDS && strcmp(cmds[cmd],argv[1]) != 0; cmd++);
    args->cmdidx = cmd < NUM_CMDS ? cmd : -1;
  }
  else args->cmdidx = -1;

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
  for(i = is_ctx_cmd ? 2 : 1; i < argc; i++)
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

void cmd_free(CmdArgs *args)
{
  free(args->argv);
  free(args->cmdline);
}

// Intended for use testing command line parsing etc.
int cmd_run(int argc, char **argv)
{
  CmdArgs args;
  cmd_alloc(&args, argc, argv);
  int ret = -1;

  if(args.cmdidx == -1) warn("Unrecognised command: %s", argc < 2 ? NULL : argv[1]);
  else ret = ctx_funcs[args.cmdidx](&args);

  cmd_free(&args);
  return ret;
}
