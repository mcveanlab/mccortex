#ifndef ALL_TESTS_H_
#define ALL_TESTS_H_

#include "db_graph.h"
#include "db_alignment.h"
#include "correct_alignment.h"
#include "build_graph.h"
#include "generate_paths.h"

#include "seq_file.h"

//
// Output
//
// Set this to redirect test status output
extern FILE *ctx_tst_out;

// Create our own output function
#define test_status(fmt,...) do {                                              \
  pthread_mutex_lock(&biglock);                                                \
  ftimestamp(ctx_tst_out);                                                     \
  fprintf(ctx_tst_out, "[%s:%i]", BASE_FILE_NAME, __LINE__);                   \
  if(((const char*)(fmt))[0] != '[') fputc(' ', ctx_tst_out);                  \
  fprintf(ctx_tst_out, fmt, ##__VA_ARGS__);                                    \
  if(((const char*)(fmt))[strlen(fmt)-1] != '\n') fputc('\n', ctx_tst_out);    \
  pthread_mutex_unlock(&biglock);                                              \
} while(0)

//
// Testing functions + MACROs
//
extern size_t tests_num_run, tests_num_failed;

// Test MACROs
#define TASSERT2(x,fmt,...) do {                                               \
  tests_num_run++;                                                             \
  if(!(x)) {                                                                   \
    tests_num_failed++;                                                        \
    fprintf(ctx_tst_out,"[%s:%i] Failed: %s\n",__FILE__,__LINE__,QUOTE_VALUE(x));\
    if((fmt) != NULL) test_status(fmt, ##__VA_ARGS__);                         \
  }                                                                            \
} while(0)

#define TASSERT(x,...) do {                                                    \
  tests_num_run++;                                                             \
  if(!(x)) {                                                                   \
    tests_num_failed++;                                                        \
    fprintf(ctx_tst_out,"[%s:%i] Failed: %s\n",__FILE__,__LINE__,QUOTE_VALUE(x));\
  }                                                                            \
} while(0)

//
// Useful functions
//
void rand_bytes(uint8_t *arr, size_t n);
void rand_nucs(Nucleotide *nucs, size_t len);
void rand_bases(char *bases, size_t len);
void bitarr_tostr(const uint8_t *arr, size_t len, char *str);

static inline void seq_read_set(read_t *r, const char *s) {
  size_t len = strlen(s);
  buffer_ensure_capacity(&r->seq, len+1);
  memcpy(r->seq.b, s, len+1);
  r->seq.end = len;
}

//
// Graph setup
//

void _construct_graph_with_paths(dBGraph *graph,
                                 size_t kmer_size, size_t ncols,
                                 char **seqs, size_t nseqs,
                                 CorrectAlnParam path_params);

void _test_add_paths(dBGraph *graph,
                     AsyncIOData *iodata, CorrectAlnInput *task,
                     GenPathWorker *wrkrs, char *seq,
                     size_t exp_npaths, size_t exp_nkmers, size_t exp_pbytes);

static inline void _tests_add_to_graph(dBGraph *graph, const char *str, size_t colour)
{
  build_graph_from_str_mt(graph, colour, str, strlen(str));
}

//
// Graph tests
//

void _test_path_store(const dBGraph *graph);

//
// Functions of tests
//

// util_tests.c
void test_util();

// dna_tests.c
void test_dna_functions();

// bkmer_tests.c
void test_bkmer_functions();

// binary_seq_tests.c
void test_binary_seq_functions();

// hash_table_tests.c
void test_hash_table();

// db_node_tests.c
void test_db_node();

// build_graph_tests.c
void test_build_graph();

// supernode_tests.c
void test_supernode();

// cleaning_tests.c
void test_cleaning();

// subgraph_tests.c
void test_subgraph();

// packed_path_tests.c
void test_packed_path();

// path_tests.c
void test_paths();

// corrected_aln_tests.c
void test_corrected_aln();

// repeat_walker_tests.c
void test_repeat_walker();

// path_set_tests.c
void test_path_sets();

// bubble_caller_tests.c
void test_bubble_caller();

// kmer_occur_tests.c
void test_kmer_occur();

// graph_crawler_tests.c
void test_graph_crawler();

// infer_edges_tests.c
void test_infer_edges_tests();

#endif  /* ALL_TESTS_H_ */
