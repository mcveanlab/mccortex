#ifndef CMD_H_
#define CMD_H_

typedef struct
{
  char *cmdline;
  boolean genome_size_set, num_kmers_set, mem_to_use_set;
  size_t genome_size, num_kmers, mem_to_use;
  boolean kmer_size_set, num_threads_set;
  uint32_t kmer_size, num_threads;
  // arguments not including command:
  int argc;
  char **argv;
} CmdArgs;

#define CMD "ctx"MACRO2STR(MAX_KMER_SIZE)

int ctx_build(CmdArgs *args);
int ctx_clean(CmdArgs *args);
int ctx_join(CmdArgs *args);
int ctx_intersect(CmdArgs *args);
int ctx_subgraph(CmdArgs *args);
int ctx_reads(CmdArgs *args);
int ctx_extend(CmdArgs *args);
int ctx_contigs(CmdArgs *args);
int ctx_thread(CmdArgs *args);
int ctp_view(CmdArgs *args);
int ctp_merge(CmdArgs *args);
int ctx_call(CmdArgs *args);
int ctx_diverge(CmdArgs *args);
int ctx_unique(CmdArgs *args);
int ctx_covg(CmdArgs *args);
int ctx_place(CmdArgs *args);

#endif
