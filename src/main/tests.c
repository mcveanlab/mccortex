#include "global.h"
#include "all_tests.h"

int main()
{
  cortex_init();
  ctx_msg_out = NULL;
  ctx_tst_out = stdout;
  test_status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  test_status("Tests running...");
  test_bkmer_functions();
  test_dna_functions();
  test_db_node();
  test_util();
  test_packed_path();
  test_status("Finished tests.");
  cortex_destroy();
  return EXIT_SUCCESS;
}
