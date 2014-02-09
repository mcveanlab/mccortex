#include "global.h"
#include "all_tests.h"
#include "util.h"

int main()
{
  cortex_init();
  ctx_msg_out = NULL;
  ctx_tst_out = stdout;
  test_status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  test_status("Tests running...");

  // Call tests
  test_util();
  test_dna_functions();
  test_binary_seq_functions();
  test_bkmer_functions();
  test_hash_table();
  test_db_node();
  test_build_graph();
  test_supernode();
  test_packed_path();
  test_paths();
  test_corrected_aln();

  // Finished
  char num_test_str[100], num_passed_str[100];
  size_t tests_num_passed = tests_num_run - tests_num_failed;
  ulong_to_str(tests_num_run, num_test_str);
  ulong_to_str(tests_num_passed, num_passed_str);

  test_status("Tests passed: %s / %s (%.1f%%)", num_passed_str, num_test_str,
              (100.0*tests_num_passed)/tests_num_run);
  test_status(tests_num_failed ? "Some tests failed." : "All tests passed.");

  cortex_destroy();
  return EXIT_SUCCESS;
}
