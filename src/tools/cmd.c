#include "global.h"
#include "cmd.h"
#include "util.h"
#include "hash_table.h" // for calculating mem usage

const char *cmds[NUM_CMDS]
  = {"build", "view", "clean", "join", "subgraph", "reads", "extend",
     "contigs", "thread", "pview", "pmerge", "call", "diverge", "unique",
     "covg", "place"};

int (*ctx_funcs[NUM_CMDS])(CmdArgs *cmd_args)
  = {ctx_build, ctx_view, ctx_clean, ctx_join, ctx_subgraph,
     ctx_reads, ctx_extend, ctx_contigs, ctx_thread, ctx_pview, ctx_pmerge,
     ctx_call, ctx_diverge, ctx_unique, ctx_covg, ctx_place};

// static int ctx_notimpl(CmdArgs *args) {
//   warn("Command not implemented [cmd: %s]", args->cmdline); return EXIT_FAILURE;
// }

void cmd_accept_options(const CmdArgs *args, const char *accptopts)
{
  if(accptopts == NULL) return;
  if(args->mem_to_use_set && strchr(accptopts,'m') == NULL)
    die("-m <memory> argument not valid for this command");
  if(args->genome_size_set && strchr(accptopts,'g') == NULL)
    die("-g <genomesize> argument not valid for this command");
  if(args->num_threads_set && strchr(accptopts,'t') == NULL)
    die("-t <threads> argument not valid for this command");
  if(args->num_kmers_set && strchr(accptopts,'h') == NULL)
    die("-h <hash-entries> argument not valid for this command");
  if(args->kmer_size_set && strchr(accptopts,'k') == NULL)
    die("-k <kmer-size> argument not valid for this command");
  if(args->file_set && strchr(accptopts,'f') == NULL)
    die("-f <file> argument not valid for this command");
}

void cmd_require_options(const CmdArgs *args, const char *requireopts,
                         const char *usage)
{
  if(requireopts == NULL) return;
  if(args->argc == 0 && *requireopts != '\0') print_usage(usage, NULL);
  for(; *requireopts != '\0'; requireopts++)
  {
    if(*requireopts == 'm') {
      if(!args->mem_to_use_set)
        die("-m <memory> argument required for this command");
    }
    else if(*requireopts == 'g') {
      if(!args->genome_size_set)
        die("-g <genomesize> argument required for this command");
    }
    else if(*requireopts == 't') {
      if(!args->genome_size_set)
        die("-t <threads> argument required for this command");
    }
    else if(*requireopts == 'h') {
      if(!args->num_kmers_set)
        die("-h <hash-entries> argument required for this command");
    }
    else if(*requireopts == 'k') {
      if(!args->kmer_size_set)
        die("-k <kmer-size> argument required for this command");
    }
    else if(*requireopts == 'f') {
      if(!args->file_set)
        die("-f <file> argument required for this command");
    }
    else warn("Ignored required option: %c", *requireopts);
  }
}

void cmd_alloc(CmdArgs *args, int argc, char **argv)
{
  args->genome_size = 3100000000;
  args->num_kmers = 4UL << 20; // 4M
  args->mem_to_use = 1UL << 30; // 1G
  args->num_threads = 2;
  args->kmer_size = 31;
  args->file = NULL;

  args->mem_to_use_set = false;
  args->genome_size_set = false;
  args->num_threads_set = false;
  args->num_kmers_set = false;
  args->kmer_size_set = false;
  args->file_set = false;

  // Get command index
  boolean is_ctx_cmd = (strstr(argv[0],"ctx") != NULL);
  if(argc >= 2 && is_ctx_cmd) {
    int cmd;
    for(cmd = 0; cmd < NUM_CMDS && strcmp(cmds[cmd],argv[1]) != 0; cmd++);
    args->cmdidx = cmd < NUM_CMDS ? cmd : -1;
  }
  else args->cmdidx = -1;

  args->argc = 0;
  args->argv = malloc2(argc * sizeof(char**));

  size_t len = argc - 1; // spaces
  int i;
  for(i = 0; i < argc; i++) len += strlen(argv[i]);
  args->cmdline = malloc2((len+1) * sizeof(char));
  sprintf(args->cmdline, "%s", argv[0]);
  for(i = 1, len = strlen(argv[0]); i < argc; len += strlen(argv[i]) + 1, i++)
    sprintf(args->cmdline+len, " %s", argv[i]);

  // argv[0] = ctx, argv[1] = cmd, argv[2..] = args
  for(i = is_ctx_cmd ? 2 : 1; i < argc; i++)
  {
    if(strcmp(argv[i], "-m") == 0)
    {
      if(i + 1 == argc) die("-m <memory> requires an argument");
      if(!mem_to_integer(argv[i+1], &args->mem_to_use) || args->mem_to_use == 0)
        die("Invalid memory argument: %s", argv[i+1]);
      args->mem_to_use_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0)
    {
      if(i + 1 == argc) die("-g <genome-size> requires an argument");
      if(!bases_to_integer(argv[i+1], &args->genome_size) || args->genome_size == 0)
        die("Invalid genome size argument: %s", argv[i+1]);
      args->genome_size_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0)
    {
      if(i + 1 == argc) die("-t <threads> requires an argument");
      if(!parse_entire_uint(argv[i+1], &args->num_threads) || args->num_threads == 0)
        die("Invalid number of threads: %s", argv[i+1]);
      args->num_threads_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-h") == 0)
    {
      if(i + 1 == argc) die("-h <hash-kmers> requires an argument");
      if(!mem_to_integer(argv[i+1], &args->num_kmers) || args->num_kmers == 0)
        die("Invalid hash size: %s", argv[i+1]);
      args->num_kmers_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-k") == 0)
    {
      if(i + 1 == argc) die("-k <kmer-size> requires an argument");
      if(!parse_entire_uint(argv[i+1], &args->kmer_size) ||
         args->kmer_size < 3 || !(args->kmer_size & 0x1))
        die("kmer size (-k) must be an odd int >= 3: %s", argv[i+1]);
      if(args->kmer_size < MIN_KMER_SIZE || args->kmer_size > MAX_KMER_SIZE)
        die("Please recompile with correct kmer size (kmer_size: %u)", args->kmer_size);
      args->kmer_size_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0)
    {
      if(i + 1 == argc) die("-f <file> requires an argument");
      args->file = argv[i+1];
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

// If your command accepts -h <kmers> and -m <mem> this may be useful
// mem_per_kmer is additional memory per node, above hash table for BinaryKmers
size_t cmd_get_kmers_in_hash(CmdArgs *args, size_t mem_per_kmer)
{
  size_t req_kmers
    = args->num_kmers_set ? args->num_kmers
                          : args->mem_to_use / (sizeof(BinaryKmer) + mem_per_kmer);

  size_t kmers_in_hash, hash_mem, graph_mem;
  hash_mem = hash_table_mem2(req_kmers, &kmers_in_hash);
  graph_mem = hash_mem + kmers_in_hash * mem_per_kmer;

  char graph_mem_str[100], mem_to_use_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(args->mem_to_use, 1, mem_to_use_str);

  if(args->mem_to_use_set && args->num_kmers_set) {
    if(graph_mem > args->mem_to_use) {

      die("-h <kmers> requires more memory than given with -m <mem> [%s > %s]",
          graph_mem_str, mem_to_use_str);
    }
    else
      message("Note: Using less memory than requested (%s < %s), due to: -h <kmer>",
              graph_mem_str, mem_to_use_str);
  }

  message("[memory]  graph: %s\n", graph_mem_str);

  return kmers_in_hash;
}
