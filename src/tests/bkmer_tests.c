#include "global.h"
#include "all_tests.h"
#include "binary_kmer.h"

void test_bkmer_to_from_str()
{
  status("Testing binary_kmer_[to|from]_str...");

  char input[MAX_KMER_SIZE+1], result[MAX_KMER_SIZE+1];
  BinaryKmer bkmer;
  size_t k = MAX_KMER_SIZE;

  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k++)
  {
    // Randomise sequence
    rnd_seq(input, k);
    bkmer = binary_kmer_from_str(input, k);
    binary_kmer_to_str(bkmer, k, result);

    if(strcmp(input, result) != 0)
      die("%s vs %s [k=%zu]", input, result, k);
    if(binary_kmer_oversized(bkmer, k))
      die("Oversized kmer %s [k=%zu]", input, k);
  }
}
