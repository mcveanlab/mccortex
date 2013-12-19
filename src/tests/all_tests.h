#ifndef ALL_TESTS_H_
#define ALL_TESTS_H_

#include "db_graph.h"

// test_utils.c: Common utitlities for testing
void rnd_seq(char *seq, size_t len);

// bkmer_tests.c
void test_bkmer_functions();

// dna_tests.c
void test_dna_functions();

// path_tests.c
void test_path_stores_match(const dBGraph *dbg1, const dBGraph *dbg2);

#endif  /* ALL_TESTS_H_ */
