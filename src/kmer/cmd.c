#include "global.h"
#include "cmd.h"
#include "util.h"
#include "hash_table.h" // for calculating mem usage

#include "misc/mem_size.h" // in libs/misc/

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
    print_usage(usage, "-c <ncols> argument not valid for this command");
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
    else if(*requireopts == 'c') {
      if(!args->use_ncols_set)
        die("-c <ncols> argument required for this command");
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
  args->argv = ctx_malloc((size_t)argc * sizeof(char**));

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
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--ncols") == 0)
    {
      if(i + 1 == argc) die("%s <ncols> requires an argument", argv[i]);
      if(args->use_ncols_set) die("-c <ncols> given more than once");
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

void cmd_print_mem(size_t mem_bytes, const char *name)
{
  char mem_str[100];
  bytes_to_str(mem_bytes, 1, mem_str);
  status("[memory] %s: %s", name, mem_str);
}

// If your command accepts -n <kmers> and -m <mem> this may be useful
//  `extra_bits` is additional memory per node, above hash table+BinaryKmers
//  `use_mem_limit` if true, fill args->mem_to_use
size_t cmd_get_kmers_in_hash(const CmdArgs *args, size_t extra_bits,
                             size_t min_num_kmer_req, size_t max_num_kmers_req,
                             bool use_mem_limit, size_t *graph_mem_ptr)
{
  size_t kmers_in_hash, graph_mem = 0, min_kmers_mem;
  char graph_mem_str[100], mem_to_use_str[100];
  char kmers_in_hash_str[100], min_num_kmers_str[100], min_kmers_mem_str[100];

  if(!use_mem_limit && min_num_kmer_req == 0 && !args->num_kmers_set) {
    cmd_print_usage("Cannot read from stream without -n <nkmers> set");
  }

  if(args->num_kmers_set)
    graph_mem = hash_table_mem(args->num_kmers, extra_bits, &kmers_in_hash);
  else if(use_mem_limit)
    graph_mem = hash_table_mem_limit(args->mem_to_use, extra_bits, &kmers_in_hash);
  else if(min_num_kmer_req > 0)
    graph_mem = hash_table_mem((size_t)(min_num_kmer_req/IDEAL_OCCUPANCY),
                               extra_bits, &kmers_in_hash);

  if(max_num_kmers_req > 0 && !args->num_kmers_set)
    kmers_in_hash = MIN2(kmers_in_hash, max_num_kmers_req/IDEAL_OCCUPANCY);

  if(kmers_in_hash < 1024)
    graph_mem = hash_table_mem(1024, extra_bits, &kmers_in_hash);
  // ^ 1024 is a very small default hash table capacity

  min_kmers_mem = hash_table_mem(min_num_kmer_req, extra_bits, NULL);

  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(args->mem_to_use, 1, mem_to_use_str);
  bytes_to_str(min_kmers_mem, 1, min_kmers_mem_str);

  ulong_to_str(kmers_in_hash, kmers_in_hash_str);
  ulong_to_str(min_num_kmer_req, min_num_kmers_str);

  // Give a error/warning about occupancy
  if(kmers_in_hash < min_num_kmer_req)
  {
    die("Not enough kmers in hash: require at least %s kmers (min memory: %s)",
        min_num_kmers_str, min_kmers_mem_str);
  }
  else if(kmers_in_hash < min_num_kmer_req/WARN_OCCUPANCY)
  {
    warn("Expected hash table occupancy %.2f%% "
         "(you may want to increase -n or -m)",
         (100.0 * min_num_kmer_req) / kmers_in_hash);
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
    die("Not enough memory for requested graph: require at least %s [>%s]",
        graph_mem_str, mem_to_use_str);
  }

  cmd_print_mem(graph_mem, "graph");

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

  status("[memory] total: %s of %s RAM\n", memstr, ramstr);
}

const char *cmd_usage = "No usage set";

void cmd_print_usage(const char *errfmt,  ...)
{
  pthread_mutex_lock(&biglock); // lock if never released

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
