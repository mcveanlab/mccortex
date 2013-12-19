#include "global.h"
#include "all_tests.h"

int main()
{
  ctx_msg_out = stdout;
  seed_random();
  status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  status("Tests running...");
  test_bkmer_to_from_str();
  test_dna_functions();
  status("Finished tests.");
}
