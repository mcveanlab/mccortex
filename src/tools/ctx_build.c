#include "global.h"

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"

static const char usage[] =
"usage: ctx_build [OPTIONS] <kmer_size> <mem> <out.ctx>\n"
"  Build a cortex binary\n"
"\n"
"  Sequence options:\n"
"    --quality_score_threshold <qual>\n"
"      Filter for quality scores in the input file. [default: 0]\n"
"    --remove_pcr_duplicates\n"
"      Removes PCR duplicate reads by ignoring read pairs if both reads start at\n"
"      the same k-mer as a previous read\n"
"    --cut_homopolymers <bases>\n"
"      Breaks reads at homopolymers of length >= <bases>\n"
"      (i.e. max homopolymer in filtered read == threshold-1) [default: off]\n"
"\n"
"    --nc  indicates a new colour is used\n"
"\n"
"  Loading sequence:\n"
"    --se_list <se.list>               Load a list of FASTA/Q/BAM\n"
"    --pe_list <pe.list1> <pe.list2>   Load paired-end data\n"
"    --seq <in.fa|fq|sam> [in2.fq]     Load sequence data\n"
"      Consecutive sequence options are loaded into the same colour\n"
"\n"
"  Loading binaries\n"
"    --load_binary <data.colours>      Load a binary into new colour(s)\n"
"\n"
"  Loading multiple colours"
"    --colour_list <data.colours>      Load a colour list into new colour(s)\n";

int main(int argc, char **argv)
{
  if(argc < 6) print_usage(usage, NULL);

  const char *out_path = argv[argc-1];
  uint32_t kmer_size, colours_used = 0;
  size_t mem_to_use;

  if(!mem_to_integer(argv[argc-2], &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", argv[argc-2]);

  if(!parse_entire_uint(argv[argc-3], &kmer_size))
    print_usage(usage, "Invalid kmer_size argument: %s", argv[argc-3]);

  if(!(kmer_size & 0x1))
    die("kmer size must be odd");
  if(kmer_size < MIN_KMER_SIZE || kmer_size > MAX_KMER_SIZE)
    die("Please compile for kmer size %u", kmer_size);

  // Validate arguments
  int argi, argend = argc-3;
  uint32_t tmp;
  for(argi = 1; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--quality_score_threshold") == 0) {
      if(argi + 1 >= argend)
        print_usage("--quality_score_threshold requires an arg", NULL);
      if(!parse_entire_uint(argv[argi+1], &tmp))
        die("Invalid --quality_score_threshold argument: %s", argv[argi+1]);
    }
    else if(strcmp(argv[argi],"--remove_pcr_duplicates") == 0) {
    }
    else if(strcmp(argv[argi],"--cut_homopolymers") == 0) {
      if(argi + 1 >= argend)
        print_usage("--cut_homopolymers requires an arg", NULL);
      if(!parse_entire_uint(argv[argi+1], &tmp))
        die("Invalid --cut_homopolymers argument: %s", argv[argi+1]);
    }
    else if(strcmp(argv[argi],"--se_list") == 0) {
      if(argi + 1 >= argend)
        print_usage("--se_list requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --se_list file: %s", argv[argi+1]);
    }
    else if(strcmp(argv[argi],"--pe_list") == 0) {
      if(argi + 2 >= argend)
        print_usage("--pe_list requires two args", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --pe_list file: %s", argv[argi+1]);
      if(!test_file_readable(argv[argi+2]))
        die("Cannot read --pe_list file: %s", argv[argi+2]);
    }
    else if(strcmp(argv[argi],"--seq") == 0) {
      if(argi + 1 >= argend)
        print_usage("--seq requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --seq file: %s", argv[argi+1]);
    }
    else if(strcmp(argv[argi],"--load_binary") == 0) {
      if(argi + 1 >= argend)
        print_usage("--load_binary requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --load_binary file: %s", argv[argi+1]);
      // DEV: get number of colours
    }
    else if(strcmp(argv[argi],"--colour_list") == 0) {
      if(argi + 1 >= argend)
        print_usage("--colour_list requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --colour_list file: %s", argv[argi+1]);
      // DEV: get number of colours
    }
    else if(strcmp(argv[argi],"--nc") == 0) {
      // DEV: increment number of colours if used current colour
    }
    else {
      print_usage(usage, "Unknown command: %s", argv[argi]);
    }
  }

  if(colours_used > NUM_OF_COLOURS)
    die("Please compile for %u colours", NUM_OF_COLOURS);

  if(!test_file_writable(out_path))
    die("Cannot write to file: %s", out_path);

  // Initialise graph, covgs, prefs
  dBGraph *db_graph;

  // Parse arguments, load
  int quality_score_threshold = -1, cut_homopolymers = -1;


  // Dump graph
  binary_dump_graph(out_path, db_graph, CURR_CTX_VERSION, NULL, 0, colours_used);

  message("Dumped version %u cortex binary to: %s", CURR_CTX_VERSION, out_path);

  return EXIT_SUCCESS;
}
