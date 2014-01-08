#include "global.h"
#include "cmd.h"
#include "util.h"
#include "hash_table.h" // for calculating mem usage

#include "mem_size.h" // in libs/misc/

const char *cmds[NUM_CMDS]
  = {"build", "view", "healthcheck", "clean", "join", "supernodes", "subgraph",
     "reads", "extend", "contigs", "inferedges", "thread", "pview", "pjoin",
     "call", "diverge", "unique", "covg", "place"};

int (*ctx_funcs[NUM_CMDS])(CmdArgs *cmd_args)
  = {ctx_build, ctx_view, ctx_health_check, ctx_clean, ctx_join, ctx_supernodes,
     ctx_subgraph, ctx_reads, ctx_extend, ctx_contigs, ctx_infer_edges,
     ctx_thread, ctx_pview, ctx_pjoin, ctx_call, ctx_diverge, ctx_unique,
     ctx_covg, ctx_place};

// Not implemented functions
static int ctx_notimpl(CmdArgs *args) {
  warn("Command not implemented [cmd: %s]", args->cmdline); return EXIT_FAILURE;
}

int ctx_covg(CmdArgs *args) { return ctx_notimpl(args); }
int ctx_diverge(CmdArgs *args) { return ctx_notimpl(args); }

void cmd_accept_options(const CmdArgs *args, const char *accptopts,
                        const char *usage)
{
  assert(accptopts != NULL);
  assert(usage != NULL);
  if(args->print_help) print_usage(usage, NULL);
  if(args->num_kmers_set && strchr(accptopts,'n') == NULL)
    print_usage(usage, "-n <hash-entries> argument not valid for this command");
  if(args->mem_to_use_set && strchr(accptopts,'m') == NULL)
    print_usage(usage, "-m <memory> argument not valid for this command");
  if(args->kmer_size_set && strchr(accptopts,'k') == NULL)
    print_usage(usage, "-k <kmer-size> argument not valid for this command");
  if(args->num_threads_set && strchr(accptopts,'t') == NULL)
    print_usage(usage, "-t <threads> argument not valid for this command");
  if(args->use_ncols_set && strchr(accptopts,'c') == NULL)
    print_usage(usage, "-c <ncols> argument not valid for this command");
  if(args->input_file_set && strchr(accptopts,'f') == NULL)
    print_usage(usage, "-f <file> argument not valid for this command");
  if(args->output_file_set && strchr(accptopts,'o') == NULL)
    print_usage(usage, "-o <file> argument not valid for this command");
  if(args->num_ctp_files > 0 && strchr(accptopts,'p') == NULL)
    print_usage(usage, "-p <in.ctp> argument not valid for this command");
  // Check for programmer error - all options should be valid
  const char *p;
  for(p = accptopts; *p != '\0'; p++) {
    if(strchr("nmktcfop", *p) == NULL)
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
    else if(*requireopts == 'k') {
      if(!args->kmer_size_set)
        die("-k <kmer-size> argument required for this command");
    }
    else if(*requireopts == 't') {
      if(!args->num_threads_set)
        die("-t <threads> argument required for this command");
    }
    else if(*requireopts == 'c') {
      if(!args->use_ncols_set)
        die("-c <ncols> argument required for this command");
    }
    else if(*requireopts == 'f') {
      if(!args->input_file_set)
        die("-f <file> argument required for this command");
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

  args->ctp_files = malloc2((size_t)argc * sizeof(char*));
  args->num_ctp_files = 0;

  // Get command index
  boolean is_ctx_cmd = (strstr(argv[0],"ctx") != NULL);
  if(argc >= 2 && is_ctx_cmd) {
    int cmd;
    for(cmd = 0; cmd < NUM_CMDS && strcmp(cmds[cmd],argv[1]) != 0; cmd++);
    args->cmdidx = cmd < NUM_CMDS ? cmd : -1;
  }
  else args->cmdidx = -1;

  args->argc = 0;
  args->argv = malloc2((size_t)argc * sizeof(char**));

  // Get cmdline string
  size_t len = (size_t)argc - 1; // spaces
  int i;
  for(i = 0; i < argc; i++) len += strlen(argv[i]);
  args->cmdline = malloc2((len+1) * sizeof(char));
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
    else if(strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0)
    {
      if(i + 1 == argc) die("%s <kmer-size> requires an argument", argv[i]);
      if(args->num_kmers_set) die("-k <kmer-size> given more than once");
      if(!parse_entire_size(argv[i+1], &args->kmer_size) ||
         args->kmer_size < 3 || !(args->kmer_size & 0x1)) {
        die("kmer size (-k) must be an odd int >= 3: %s", argv[i+1]);
      }
      if(args->kmer_size < MIN_KMER_SIZE || args->kmer_size > MAX_KMER_SIZE) {
        die("Please recompile with correct kmer size (kmer_size: %zu)",
            args->kmer_size);
      }
      args->kmer_size_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--threads") == 0)
    {
      if(i + 1 == argc) die("%s <threads> requires an argument", argv[i]);
      if(args->num_threads_set) die("-t <threads> given more than once");
      if(!parse_entire_size(argv[i+1], &args->num_threads) || args->num_threads == 0)
        die("Invalid number of threads: %s", argv[i+1]);
      args->num_threads_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--ncols") == 0)
    {
      if(i + 1 == argc) die("%s <ncols> requires an argument", argv[i]);
      if(args->num_threads_set) die("-t <ncols> given more than once");
      if(!parse_entire_size(argv[i+1], &args->use_ncols) || args->use_ncols == 0)
        die("Invalid number of colours to use: %s", argv[i+1]);
      args->use_ncols_set = true;
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--file") == 0)
    {
      if(i + 1 == argc) die("%s <file> requires an argument", argv[i]);
      if(args->input_file_set) die("%s <file> given more than once", argv[i]);
      args->input_file = argv[i+1];
      args->input_file_set = true;
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
  free(args->ctp_files);
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

// If your command accepts -n <kmers> and -m <mem> this may be useful
// extra_bits_per_kmer is additional memory per node, above hash table for
// BinaryKmers
size_t cmd_get_kmers_in_hash(const CmdArgs *args, size_t extra_bits_per_kmer,
                             size_t min_num_kmers, boolean use_mem_limit,
                             size_t *graph_mem_ptr)
{
  size_t bits_per_kmer, req_kmers;
  size_t kmers_in_hash, hash_mem, graph_mem, min_kmers_mem;
  char graph_mem_str[100], mem_to_use_str[100];
  char kmers_in_hash_str[100], min_num_kmers_str[100], min_kmers_mem_str[100];
  boolean above_nkmers = false;

  bits_per_kmer = 8*sizeof(BinaryKmer) + extra_bits_per_kmer;

  if(args->num_kmers_set) {
    req_kmers = args->num_kmers;
    above_nkmers = true;
  }
  else if(use_mem_limit && args->mem_to_use_set)
    req_kmers = (8 * args->mem_to_use) / bits_per_kmer;
  else if(min_num_kmers > 0) {
    req_kmers = (size_t)(min_num_kmers / IDEAL_OCCUPANCY);
    above_nkmers = true;
  }
  else {
    req_kmers = args->num_kmers; // take default
    above_nkmers = true;
  }

  hash_mem = hash_table_mem(req_kmers, above_nkmers, &kmers_in_hash);
  graph_mem = hash_mem + (kmers_in_hash * extra_bits_per_kmer) / 8;
  min_kmers_mem = hash_table_mem(min_num_kmers, true, NULL) +
                  (min_num_kmers*bits_per_kmer)/8;

  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(args->mem_to_use, 1, mem_to_use_str);
  ulong_to_str(kmers_in_hash, kmers_in_hash_str);
  bytes_to_str(min_kmers_mem, 1, min_kmers_mem_str);
  ulong_to_str(min_num_kmers, min_num_kmers_str);

  if(kmers_in_hash < min_num_kmers) {
    die("Not enough kmers in hash: require at least %s kmers (min memory: %s)",
        min_num_kmers_str, min_kmers_mem_str);
  }

  // Give a warning about occupancy unless explicitly set to exactly 100%
  if(kmers_in_hash < min_num_kmers/WARN_OCCUPANCY &&
     kmers_in_hash != args->num_kmers)
  {
    warn("Expected hash table occupancy %.2f%% "
         "(you may want to increase -n or -m)",
         (100.0 * min_num_kmers) / kmers_in_hash);
  }

  if(args->mem_to_use_set && args->num_kmers_set) {
    if(args->num_kmers > kmers_in_hash) {
      die("-n <kmers> requires more memory than given with -m <mem> [%s > %s]",
          graph_mem_str, mem_to_use_str);
    }
    else if(use_mem_limit && graph_mem < args->mem_to_use) {
      status("Note: Using less memory than requested (%s < %s); allows for %s kmers",
             graph_mem_str, mem_to_use_str, kmers_in_hash_str);
    }
  }

  if(graph_mem > args->mem_to_use) {
    die("Not enough memory for graph: require at least %s", min_kmers_mem_str);
  }

  status("[memory] graph: %s\n", graph_mem_str);

  if(graph_mem_ptr != NULL) *graph_mem_ptr = graph_mem;

  return kmers_in_hash;
}

// Check memory against memory limit and machine memory
void cmd_check_mem_limit(const CmdArgs *args, size_t mem_requested)
{
  char memstr[100], ramstr[100];
  bytes_to_str(mem_requested, 1, memstr);

  if(mem_requested > args->mem_to_use) {
    if(args->mem_to_use_set)
      die("Need to set higher memory limit [ at least -m %s ]", memstr);
    else
      die("Please set a memory limit [ at least -m %s ]", memstr);
  }

  // Get memory
  size_t ram = getMemorySize();
  bytes_to_str(ram, 1, ramstr);

  if(mem_requested > ram) {
    die("Requesting more memory than is available [ Reqeusted: -m %s RAM: %s ]",
        memstr, ramstr);
  }

  status("[memory] total: %s\n", memstr);
}
