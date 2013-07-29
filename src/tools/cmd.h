#ifndef CMD_H_
#define CMD_H_

typedef struct
{
  char *cmdline;
  int cmdidx;
  // options
  boolean genome_size_set, num_kmers_set, mem_to_use_set;
  size_t genome_size, num_kmers, mem_to_use;
  boolean kmer_size_set, num_threads_set;
  uint32_t kmer_size, num_threads;
  boolean file_set;
  const char *file;
  // arguments not including command:
  int argc;
  char **argv;
} CmdArgs;

#define CMD "ctx"QUOTE_VALUE(MAX_KMER_SIZE)

int ctx_build(CmdArgs *args);
int ctx_view(CmdArgs *args);
int ctx_clean(CmdArgs *args);
int ctx_join(CmdArgs *args);
int ctx_intersect(CmdArgs *args);
int ctx_subgraph(CmdArgs *args);
int ctx_reads(CmdArgs *args);
int ctx_extend(CmdArgs *args);
int ctx_contigs(CmdArgs *args);
int ctx_thread(CmdArgs *args);
int ctx_pview(CmdArgs *args);
int ctx_pmerge(CmdArgs *args);
int ctx_call(CmdArgs *args);
int ctx_diverge(CmdArgs *args);
int ctx_unique(CmdArgs *args);
int ctx_covg(CmdArgs *args);
int ctx_place(CmdArgs *args);

#define NUM_CMDS 17
extern const char *cmds[NUM_CMDS];
extern int (*ctx_funcs[NUM_CMDS])(CmdArgs *cmd_args);

void cmd_alloc(CmdArgs *args, int argc, char **argv);
void cmd_free(CmdArgs *args);

// accptopts is a string of valid args,
// e.g. "tk" accepts kmer-size and number of threads
// NULL means anything valid, "" means no args valid
void cmd_accept_options(const CmdArgs *args, const char *accptopts);
void cmd_require_options(const CmdArgs *args, const char *requireopts,
                         const char *usage);

// Run a command
int cmd_run(int argc, char **argv);

// If your command accepts -h <kmers> and -m <mem> this may be useful
size_t cmd_get_kmers_in_hash(CmdArgs *args, size_t mem_per_kmer);

#endif
