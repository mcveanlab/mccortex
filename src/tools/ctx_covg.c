#include "global.h"

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"

static const char usage[] =
"usage: ctx_covg <mem> <in.ctx> <calls.raw.vcf> [col1 ...]\n"
"  Clean a cortex binary.\n";;

int main(int argc, char **argv)
{
  // Set up: parse command line args
  // Set up: open binary and VCF file
  // 0) Set up hash table
  // 1) Load sequence from VCF to build hash table
  // 2) Load binary only for kmers already in hash table
  // 3) Parse VCF a second time and print covg information
  return EXIT_SUCCESS;
}
