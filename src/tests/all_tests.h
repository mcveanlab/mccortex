#ifndef ALL_TESTS_H_
#define ALL_TESTS_H_

#include "db_graph.h"

// General

// Set this to redirect test status output
extern FILE *ctx_tst_out;

// Create our own output function
#define test_status(fmt,...) fstatus(ctx_tst_out, fmt, ##__VA_ARGS__)

// bkmer_tests.c
void test_bkmer_functions();

// dna_tests.c
void test_dna_functions();

// hash_table_tests.c
void test_hash_table();

// db_node_tests.c
void test_db_node();

// util_tests.c
void test_util();

// packed_path_tests.c
void test_packed_path();

// Not finished
// path_tests.c
void test_paths();
// void test_path_stores_match(const dBGraph *dbg1, const dBGraph *dbg2);

#endif  /* ALL_TESTS_H_ */
