#ifndef CMD_H_
#define CMD_H_

#define CMD "ctx"QUOTE_VALUE(MAX_KMER_SIZE)

typedef struct
{
  char *cmdline;
  int cmdidx; // command specified
  boolean print_help;
  // options
  boolean num_kmers_set, mem_to_use_set, kmer_size_set, num_threads_set, use_ncols_set;
  size_t num_kmers, mem_to_use, kmer_size, num_threads, use_ncols;
  boolean input_file_set, output_file_set;
  char *input_file, *output_file;
  size_t num_ctp_files;
  char **ctp_files;
  // arguments not including command:
  int argc;
  char **argv;
} CmdArgs;

#define CMD_ARGS_INIT_MACRO { \
  .cmdline = NULL, .cmdidx = -1, .print_help = false, \
  .num_kmers_set = false, .num_kmers = 4UL<<20 /*4MB*/, \
  .mem_to_use_set = false, .mem_to_use = 1UL<<30 /*1G*/, \
  .kmer_size_set = false, .kmer_size = MAX_KMER_SIZE, \
  .num_threads_set = false, .num_threads = 2, \
  .use_ncols_set = false, .use_ncols = 1, \
  .input_file_set = false, .input_file = NULL, \
  .output_file_set = false, .output_file = NULL, \
  .num_ctp_files = 0, .ctp_files = NULL, \
  .argc = 0, .argv = NULL}

int ctx_build(CmdArgs *args);
int ctx_view(CmdArgs *args);
int ctx_health_check(CmdArgs *args);
int ctx_clean(CmdArgs *args);
int ctx_join(CmdArgs *args);
int ctx_supernodes(CmdArgs *args);
int ctx_subgraph(CmdArgs *args);
int ctx_reads(CmdArgs *args);
int ctx_extend(CmdArgs *args);
int ctx_contigs(CmdArgs *args);
int ctx_infer_edges(CmdArgs *args);
int ctx_thread(CmdArgs *args);
int ctx_pview(CmdArgs *args);
int ctx_pjoin(CmdArgs *args);
int ctx_call(CmdArgs *args);
int ctx_diverge(CmdArgs *args);
int ctx_unique(CmdArgs *args);
int ctx_covg(CmdArgs *args);
int ctx_place(CmdArgs *args);

#define NUM_CMDS 19
extern const char *cmds[NUM_CMDS];
extern int (*ctx_funcs[NUM_CMDS])(CmdArgs *cmd_args);

void cmd_alloc(CmdArgs *args, int argc, char **argv);
void cmd_free(CmdArgs *args);

// accptopts is a string of valid args,
// e.g. "tk" accepts kmer-size and number of threads
// NULL means anything valid, "" means no args valid
void cmd_accept_options(const CmdArgs *args, const char *accptopts,
                        const char *usage);
void cmd_require_options(const CmdArgs *args, const char *requireopts,
                         const char *usage);

// Run a command
int cmd_run(int argc, char **argv);

// If your command accepts -n <kmers> and -m <mem> this may be useful
// extra_bits_per_kmer is additional memory per node, above hash table for
// BinaryKmers
// Resulting graph_mem is always < args->mem_to_use
size_t cmd_get_kmers_in_hash(const CmdArgs *args, size_t extra_bits_per_kmer,
                             size_t min_num_kmers, boolean use_mem_limit,
                             size_t *graph_mem_ptr);

// Check memory against args->mem_to_use and total RAM
void cmd_check_mem_limit(const CmdArgs *args, size_t mem_requested);

#endif
