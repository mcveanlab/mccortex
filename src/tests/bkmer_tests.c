#include "global.h"
#include "all_tests.h"
#include "binary_kmer.h"

void test_bkmer_str()
{
  test_status("Testing binary_kmer_to_str() binary_kmer_from_str()");

  size_t k;
  BinaryKmer bkmer0, bkmer1;
  char str0[MAX_KMER_SIZE+1], str1[MAX_KMER_SIZE+1];

  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k+=2)
  {
    bkmer0 = binary_kmer_random(k);
    binary_kmer_to_str(bkmer0, k, str0);
    bkmer1 = binary_kmer_from_str(str0, k);
    binary_kmer_to_str(bkmer1, k, str1);

    TASSERT(!binary_kmer_oversized(bkmer0, k));
    TASSERT(!binary_kmer_oversized(bkmer1, k));
    TASSERT(binary_kmers_are_equal(bkmer0,bkmer1));
    TASSERT(strcmp(str0,str1) == 0);
  }
}

void test_bkmer_revcmp()
{
  test_status("Testing binary_kmer_reverse_complement()");

  size_t k;
  BinaryKmer bkmer0, bkmer1, bkmer2;

  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k+=2)
  {
    bkmer0 = binary_kmer_random(k);
    bkmer1 = binary_kmer_reverse_complement(bkmer0, k);
    bkmer2 = binary_kmer_reverse_complement(bkmer1, k);
    TASSERT(!binary_kmer_oversized(bkmer0, k));
    TASSERT(!binary_kmer_oversized(bkmer1, k));
    TASSERT(!binary_kmer_oversized(bkmer2, k));
    // kmer-size is odd, forward != reverse complement
    TASSERT(!binary_kmers_are_equal(bkmer0, bkmer1));
    TASSERT(binary_kmers_are_equal(bkmer0, bkmer2));
  }
}

void test_bkmer_shifts()
{
  test_status("Testing shifting and adding bases");

  size_t i, k;
  BinaryKmer bkmer0, bkmer1, bkmer2;
  Nucleotide nuc;

  // Test shifting and zeroing end base
  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k+=2)
  {
    bkmer0 = binary_kmer_random(k);
    TASSERT(!binary_kmer_oversized(bkmer0, k));

    bkmer1 = bkmer2 = bkmer0;
    bkmer1 = binary_kmer_left_shift_one_base(bkmer1, k);
    bkmer1 = binary_kmer_right_shift_one_base(bkmer1);
    binary_kmer_set_first_nuc(&bkmer2, 0, k);
    TASSERT(binary_kmers_are_equal(bkmer1,bkmer2));
    TASSERT(!binary_kmer_oversized(bkmer1, k));
    TASSERT(!binary_kmer_oversized(bkmer2, k));

    bkmer1 = bkmer2 = bkmer0;
    bkmer1 = binary_kmer_right_shift_one_base(bkmer1);
    bkmer1 = binary_kmer_left_shift_one_base(bkmer1, k);
    binary_kmer_set_last_nuc(&bkmer2, 0);
    TASSERT(binary_kmers_are_equal(bkmer1,bkmer2));
    TASSERT(!binary_kmer_oversized(bkmer1, k));
    TASSERT(!binary_kmer_oversized(bkmer2, k));
  }

  // Shift entirely until zero
  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k+=2)
  {
    bkmer0 = binary_kmer_random(k);
    TASSERT(!binary_kmer_oversized(bkmer0, k));

    bkmer1 = bkmer2 = bkmer0;
    for(i = 0; i < k; i++) {
      bkmer1 = binary_kmer_left_shift_one_base(bkmer1, k);
      bkmer2 = binary_kmer_right_shift_one_base(bkmer2);
    }
    TASSERT(binary_kmers_are_equal(bkmer1,zero_bkmer));
    TASSERT(binary_kmers_are_equal(bkmer2,zero_bkmer));
  }

  // Copy from one bkmer to another by shifting a base at a time
  for(k = MIN_KMER_SIZE; k <= MAX_KMER_SIZE; k+=2)
  {
    bkmer0 = binary_kmer_random(k);
    TASSERT(!binary_kmer_oversized(bkmer0, k));

    // copy from bkmer1 -> bkmer2, shifting right
    bkmer1 = bkmer0;
    bkmer2 = zero_bkmer;
    for(i = 0; i < k; i++) {
      nuc = binary_kmer_last_nuc(bkmer1);
      bkmer1 = binary_kmer_right_shift_one_base(bkmer1);
      bkmer2 = binary_kmer_right_shift_add(bkmer2, k, nuc);
      TASSERT(!binary_kmer_oversized(bkmer1, k));
      TASSERT(!binary_kmer_oversized(bkmer2, k));
    }
    TASSERT(binary_kmers_are_equal(bkmer1,zero_bkmer));
    TASSERT(binary_kmers_are_equal(bkmer2,bkmer0));

    // copy from bkmer1 -> bkmer2, shifting left
    bkmer1 = bkmer0;
    bkmer2 = zero_bkmer;
    for(i = 0; i < k; i++) {
      nuc = binary_kmer_first_nuc(bkmer1, k);
      bkmer1 = binary_kmer_left_shift_one_base(bkmer1, k);
      bkmer2 = binary_kmer_left_shift_add(bkmer2, k, nuc);
      TASSERT(!binary_kmer_oversized(bkmer1, k));
      TASSERT(!binary_kmer_oversized(bkmer2, k));
    }
    TASSERT(binary_kmers_are_equal(bkmer1,zero_bkmer));
    TASSERT(binary_kmers_are_equal(bkmer2,bkmer0));
  }
}

static void test_bkmer_last_nuc()
{
  test_status("Testing binary_kmer_last_nuc()");
  BinaryKmer bkmer = binary_kmer_from_str("CCACGTAAAGC", 11);
  TASSERT(binary_kmer_last_nuc(bkmer) == dna_char_to_nuc('C'));
  bkmer = binary_kmer_right_shift_one_base(bkmer);
  TASSERT(binary_kmer_last_nuc(bkmer) == dna_char_to_nuc('G'));
}

void test_bkmer_functions()
{
  TASSERT(sizeof(BinaryKmer) == NUM_BKMER_WORDS * 8);
  test_bkmer_str();
  test_bkmer_revcmp();
  test_bkmer_shifts();
  test_bkmer_last_nuc();
  // TODO: equal, less than, cmp
}
