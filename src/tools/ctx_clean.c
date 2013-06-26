#include "global.h"

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"

// clean by individual
static const char usage[] =
"usage: ctx_clean [OPTIONS] <mem> <genome_size> <in.ctx|ctxlist|colours> <out.ctx>\n"
"  Clean a cortex binary.\n"
"  Options:\n"
"    --tips\n"
"    --supernodes\n"
"  With no options only supernode cleaning is done\n";

int main(int argc, char **argv)
{
  if(argc < 5 || argc > 7) print_usage(usage, NULL);

  // Check cmdline args
  const char *out_path = argv[argc-1];
  const char *in_path = argv[argc-2];
  size_t genome_size, mem_to_use;

  if(!test_file_readable(in_path))
    die("Cannot read input file: %s", in_path);

  if(!bases_to_integer(argv[argc-3], &genome_size))
    die("Cannot parse genome size: %s", argv[argc-3]);

  if(!mem_to_integer(argv[argc-4], &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", argv[argc-4]);

  boolean tip_cleaning = false, supernode_cleaning = false;

  int argi, argend = argc-4;
  for(argi = 1; argi < argend; argi++) {
    if(strcmp(argv[argi],"--tips") == 0) tip_cleaning = true;
    else if(strcmp(argv[argi],"--supernodes") == 0) supernode_cleaning = true;
    else print_usage(usage, "Unknown argument: %s", argv[argi]);
  }

  // default settings
  if(!tip_cleaning && !supernode_cleaning) {
    supernode_cleaning = true;
  }

  if(!test_file_writable(in_path))
    die("Cannot read write output file: %s", out_path);

  // Probe binary file
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(in_path, &is_binary, &kmer_size, &num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read binary file: %s", in_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", in_path);

  // Check if multi-coloured binary - if so warn



  // Initiate graph

  // Load

  // Estimate covg from genome size

  // supernode cleaning:
  // Get supernode covg
  // clean supernodes

  // dump graph

  // free

  return EXIT_SUCCESS;
}
