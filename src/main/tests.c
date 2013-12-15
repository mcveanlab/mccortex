#include "global.h"
#include "bkmer_tests.h"

int main()
{
  ctx_msg_out = stdout;
  seed_random();
  status("Min kmer-size: %i, max kmer-size: %i", MIN_KMER_SIZE, MAX_KMER_SIZE);
  status("Tests running...");
  test_bkmer_to_from_str();
  status("Finished tests.");
}
