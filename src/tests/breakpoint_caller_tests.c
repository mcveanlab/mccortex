#include "global.h"
#include "all_tests.h"

#include "breakpoint_caller.h"
#include "kmer_occur.h"

static void test_kmer_occur_filter()
{
  KOccurBuffer kobuf, kobuftmp;
  kmer_occur_buf_alloc(&kobuf, 16);
  kmer_occur_buf_alloc(&kobuftmp, 16);

  // DEV
  // kograph_filter_stretch();

  kmer_occur_buf_dealloc(&kobuf);
  kmer_occur_buf_dealloc(&kobuftmp);
}

static void test_find_breakpoint()
{
  // DEV
}

// bubble_caller_tests.c
void test_breakpoint_caller()
{
  test_status("Testing breakpoint calling...");

  test_find_breakpoint();
  test_kmer_occur_filter();
}
