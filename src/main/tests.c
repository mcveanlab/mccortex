#include "global.h"
#include "all_tests.h"

int main()
{
  cortex_init();
  ctx_msg_out = NULL;
  ctx_tst_out = stdout;
  test_status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  test_status("Tests running...");
  test_util();
  test_dna_functions();
  test_bkmer_functions();
  test_hash_table();
  test_db_node();
  test_supernode();
  test_packed_path();
  test_paths();
  test_corrected_aln();
  size_t tests_num_passed = tests_num_run - tests_num_failed;
  test_status("Tests passed: %zu / %zu (%.1f%%)", tests_num_passed, tests_num_run,
              (100.0*tests_num_passed)/tests_num_run);
  test_status(tests_num_failed ? "Some tests failed." : "All tests passed.");
  cortex_destroy();
  return EXIT_SUCCESS;
}
