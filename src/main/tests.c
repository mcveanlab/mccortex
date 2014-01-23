#include "global.h"
#include "all_tests.h"

int main()
{
  cortex_init();
  ctx_msg_out = stdout;
  status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  status("Tests running...");
  test_bkmer_functions();
  test_dna_functions();
  test_db_node();
  test_util();
  status("Finished tests.");
  cortex_destroy();
  return EXIT_SUCCESS;
}
