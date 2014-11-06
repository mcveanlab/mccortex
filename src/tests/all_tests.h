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
  pthread_mutex_lock(&ctx_biglock);                                            \
  timestampf(ctx_tst_out);                                                     \
  fprintf(ctx_tst_out, "[%s:%i]", BASE_FILE_NAME, __LINE__);                   \
  if(((const char*)(fmt))[0] != '[') fputc(' ', ctx_tst_out);                  \
  fprintf(ctx_tst_out, fmt, ##__VA_ARGS__);                                    \
  if(((const char*)(fmt))[strlen(fmt)-1] != '\n') fputc('\n', ctx_tst_out);    \
  pthread_mutex_unlock(&ctx_biglock);                                          \
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

void all_tests_add_paths_multi(dBGraph *graph, const char **seqs, size_t nseqs,
                               CorrectAlnParam params,
                               int exp_npaths, int exp_nkmers);

void all_tests_add_paths(dBGraph *graph, const char *seq,
                         CorrectAlnParam params,
                         int exp_npaths, int exp_nkmers);

void all_tests_construct_graph(dBGraph *graph,
                               size_t kmer_size, size_t ncols,
                               const char **seqs, size_t nseqs,
                               CorrectAlnParam path_params);

static inline void _tests_add_to_graph(dBGraph *graph, const char *str, size_t colour)
{
  build_graph_from_str_mt(graph, colour, str, strlen(str));
}

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

// path_tests.c
void test_paths();

// path_set_tests.c
// void test_path_sets();

// graph_walker_tests.c
void test_graph_walker();

// corrected_aln_tests.c
void test_corrected_aln();

// repeat_walker_tests.c
void test_repeat_walker();

// bubble_caller_tests.c
void test_bubble_caller();

// kmer_occur_tests.c
void test_kmer_occur();

// graph_crawler_tests.c
void test_graph_crawler();

// infer_edges_tests.c
void test_infer_edges_tests();

#endif  /* ALL_TESTS_H_ */
