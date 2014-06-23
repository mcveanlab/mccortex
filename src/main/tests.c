#include "global.h"
#include "all_tests.h"
#include "cmd.h"
#include "util.h"

int main(int argc, char **argv)
{
  cortex_init();
  cmd_init(argc, argv);

  ctx_msg_out = NULL;
  ctx_tst_out = stdout;

  test_status("Tests running...");
  test_status("[version] "VERSION_STATUS_STR"\n"); // defined in global.h

  // Binary Kmer tests should work for all values of MAXK
  test_bkmer_functions();

  #if MAX_KMER_SIZE == 31
    // Call remaining tests
    test_util();
    test_dna_functions();
    test_binary_seq_functions();
    test_hash_table();
    test_db_node();
    test_build_graph();
    test_supernode();
    test_subgraph();
    test_cleaning();
    test_paths();
    // test_path_sets(); // DEV: replace with test_path_subset()
    test_corrected_aln();
    test_repeat_walker();
    test_graph_crawler();
    test_bubble_caller();
    test_kmer_occur();
    test_infer_edges_tests();
  #endif

  // Check we free'd all our memory
  size_t still_alloced = alloc_get_num_allocs() - alloc_get_num_frees();
  TASSERT2(still_alloced == 0, "%zu not free'd", still_alloced);

  // Finished
  char num_test_str[100], num_passed_str[100];
  size_t tests_num_passed = tests_num_run - tests_num_failed;
  ulong_to_str(tests_num_run, num_test_str);
  ulong_to_str(tests_num_passed, num_passed_str);

  test_status("Tests passed: %s / %s (%.1f%%)", num_passed_str, num_test_str,
              (100.0*tests_num_passed)/tests_num_run);

  if(tests_num_failed) test_status("%zu tests failed", tests_num_failed);
  else test_status("All tests passed.");

  cortex_destroy();

  // Return 1 if any tests failed, 0 on success
  return tests_num_failed ? 1 : 0;
}
