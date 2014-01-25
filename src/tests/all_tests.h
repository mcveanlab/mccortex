#ifndef ALL_TESTS_H_
#define ALL_TESTS_H_

#include "db_graph.h"

// General
extern FILE *ctx_tst_out;

void test_status(const char *fmt, ...)
__attribute__((format(printf, 1, 2)));

// bkmer_tests.c
void test_bkmer_functions();

// dna_tests.c
void test_dna_functions();

// db_node.c
void test_db_node();

// util.c
void test_util();

// packed_path.c
void test_packed_path();

// path_tests.c
void test_path_stores_match(const dBGraph *dbg1, const dBGraph *dbg2);

#endif  /* ALL_TESTS_H_ */
