#ifndef ALL_TESTS_H_
#define ALL_TESTS_H_

#include "db_graph.h"
#include "seq_file.h"

//
// Output
//
// Set this to redirect test status output
extern FILE *ctx_tst_out;

// Create our own output function
#define test_status(fmt,...) fstatus(ctx_tst_out, fmt, ##__VA_ARGS__)

//
// Testing functions + MACROs
//
extern size_t tests_num_run, tests_num_failed;

// Test MACROs
#define TASSERT2(x,fmt,...) do {                                               \
  tests_num_run++;                                                             \
  if(!(x)) { tests_num_failed++;                                               \
    fprintf(ctx_tst_out, "[%s:%i] Failed: %s\n",__FILE__,__LINE__,QUOTE_VALUE(x));\
    if((fmt) != NULL) test_status(fmt, ##__VA_ARGS__);                         \
  } \
} while(0)

#define TASSERT(x) TASSERT2(x,NULL)

//
// Useful functions
//
void fill_rand_bytes(uint8_t *arr, size_t n);
void rand_nucs(Nucleotide *nucs, size_t len);
void bitarr_tostr(const uint8_t *arr, size_t len, char *str);

static inline void seq_read_set(read_t *r, const char *s) {
  size_t len = strlen(s);
  buffer_ensure_capacity(&r->seq, len+1);
  memcpy(r->seq.b, s, len+1);
  r->seq.end = len;
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

// packed_path_tests.c
void test_packed_path();

// path_tests.c
void test_paths();

// corrected_aln_tests.c
void test_corrected_aln();

#endif  /* ALL_TESTS_H_ */
