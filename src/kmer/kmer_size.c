#include "kmer_size.h"

// Files that are not compiled with MIN_KMER_SIZE and MAX_KMER_SIZE link to
// this object file and discover kmer size limits at run time

int get_min_kmer_size()
{
  return MIN_KMER_SIZE;
}

int get_max_kmer_size()
{
  return MAX_KMER_SIZE;
}
