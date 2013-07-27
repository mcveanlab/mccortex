#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"
#include "file_reader.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" build -k <kmer-size> -m <mem> [OPTIONS] <out.ctx>\n"
"  Build a cortex binary.  \n"
"\n"
"  Sequence options:\n"
"    --quality_score_threshold <qual>\n"
"      Filter for quality scores in the input file. [default: 0]\n"
"    --fq_ascii <qual>\n"
"      FASTQ ASCII offset. [default: 0 (auto-detected)]\n"
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
"    --seq <in.fa|fq|sam>              Load sequence data\n"
"    --seq2 <in1> <in2>                Load paired end sequence data\n"
"      Consecutive sequence options are loaded into the same colour\n"
"\n"
"  Loading binaries\n"
"    --load_binary <data.colours>      Load a binary into new colour(s)\n"
"\n"
"  Loading multiple colours\n"
"    --colour_list <data.colours>      Load a colour list into new colour(s)\n";

int ctx_build(CmdArgs *args)
{
  cmd_accept_options(args, "mhk");
  cmd_require_options(args, "k");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  uint32_t kmer_size = args->kmer_size;
  size_t mem_to_use = args->mem_to_use;

  const char *out_path = argv[argc-1];
  uint32_t colours_used = 0;
  boolean current_colour_used = false;

  // Validate arguments
  int argi, argend = argc-1;
  uint32_t tmp;
  for(argi = 0; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--quality_score_threshold") == 0) {
      if(argi + 1 >= argend)
        print_usage("--quality_score_threshold requires an arg", NULL);
      if(!parse_entire_uint(argv[argi+1], &tmp) || tmp > 128)
        die("Invalid --quality_score_threshold argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_ascii") == 0) {
      if(argi + 1 >= argend)
        print_usage("--fq_ascii requires an arg", NULL);
      if(!parse_entire_uint(argv[argi+1], &tmp) || tmp > 128)
        die("Invalid --fq_ascii argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--remove_pcr_duplicates") == 0) {
    }
    else if(strcmp(argv[argi],"--cut_homopolymers") == 0) {
      if(argi + 1 >= argend)
        print_usage("--cut_homopolymers requires an arg", NULL);
      if(!parse_entire_uint(argv[argi+1], &tmp))
        die("Invalid --cut_homopolymers argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--se_list") == 0) {
      if(argi + 1 >= argend)
        print_usage("--se_list requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --se_list file: %s", argv[argi+1]);
      check_colour_or_ctx_list(argv[argi+1], false, true, true, kmer_size, NULL);
      argi += 1;
      current_colour_used = true;
    }
    else if(strcmp(argv[argi],"--pe_list") == 0) {
      if(argi + 2 >= argend)
        print_usage("--pe_list requires two args", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --pe_list file: %s", argv[argi+1]);
      if(!test_file_readable(argv[argi+2]))
        die("Cannot read --pe_list file: %s", argv[argi+2]);
      uint32_t num1, num2;
      num1 = check_colour_or_ctx_list(argv[argi+1], false, false, true, kmer_size, NULL);
      num2 = check_colour_or_ctx_list(argv[argi+2], false, false, true, kmer_size, NULL);
      if(num1 != num2)
        die("--pe_list files diff lengths: %s; %s", argv[argi+1], argv[argi+2]);
      argi += 2;
      current_colour_used = true;
    }
    else if(strcmp(argv[argi],"--seq") == 0) {
      if(argi + 1 >= argend)
        print_usage("--seq requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --seq file: %s", argv[argi+1]);
      argi += 1;
      current_colour_used = true;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      if(argi + 2 >= argend)
        print_usage("--seq2 requires two args", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read first --seq2 file: %s", argv[argi+1]);
      if(!test_file_readable(argv[argi+2]))
        die("Cannot read second --seq2 file: %s", argv[argi+2]);
      argi += 2;
      current_colour_used = true;
    }
    else if(strcmp(argv[argi],"--load_binary") == 0) {
      if(argi + 1 >= argend)
        print_usage("--load_binary requires an arg", NULL);
      // probe binary to get number of colours
      boolean is_binary = false;
      uint32_t bin_kmer_size, bin_num_of_cols, max_col;
      uint64_t bin_num_kmers;
      if(!binary_probe(argv[argi+1], &is_binary, &bin_kmer_size,
                       &bin_num_of_cols, &max_col, &bin_num_kmers)) {
        print_usage(usage, "Cannot read binary file: %s", argv[argi+1]);
      } else if(!is_binary) {
        print_usage(usage, "Input binary file isn't valid: %s", argv[argi+1]);
      } else if(bin_kmer_size != kmer_size) {
        print_usage(usage, "Input binary kmer_size doesn't match [%u vs %u]",
                    bin_kmer_size, kmer_size);
      }
      argi += 1;
      colours_used += bin_num_of_cols;
      current_colour_used = false;
    }
    else if(strcmp(argv[argi],"--colour_list") == 0) {
      if(argi + 1 >= argend)
        print_usage("--colour_list requires an arg", NULL);
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read --colour_list file: %s", argv[argi+1]);
      uint32_t cols_in_list;
      check_colour_or_ctx_list(argv[argi+1], true, true, true, kmer_size, &cols_in_list);
      argi += 1;
      colours_used += cols_in_list;
      current_colour_used = false;
    }
    else if(strcmp(argv[argi],"--nc") == 0) {
      if(!colours_used) print_usage(usage, "--nc not needed, colour is fresh");
      colours_used++;
      current_colour_used = false;
    }
    else {
      print_usage(usage, "Unknown command: %s", argv[argi]);
    }
  }

  if(current_colour_used) colours_used++;

  if(!test_file_writable(out_path))
    die("Cannot write to file: %s", out_path);

  // Initialise graph, covgs, prefs
  dBGraph db_graph;

  size_t mem_per_kmer, req_kmers;

  mem_per_kmer = sizeof(BinaryKmer) + (sizeof(Covg) + sizeof(Edges)) * colours_used;
  req_kmers = args->num_kmers_set ? args->num_kmers : mem_to_use / mem_per_kmer;

  size_t kmers_in_hash, hash_mem, graph_mem;
  hash_mem = hash_table_mem2(req_kmers, &kmers_in_hash);
  graph_mem = kmers_in_hash * mem_per_kmer;

  char graph_mem_str[100], mem_to_use_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(mem_to_use, 1, mem_to_use_str);

  if(args->mem_to_use_set && args->num_kmers_set) {
    if(graph_mem > args->mem_to_use) {

      die("-h <kmers> requires more memory than given with -m <mem> [%s > %s]",
          graph_mem_str, mem_to_use_str);
    }
    else
      message("Note: Using less memory than requested, due to: -h <kmer>");
  }

  db_graph_alloc(&db_graph, kmer_size, colours_used, kmers_in_hash);
  db_graph.col_edges = calloc(kmers_in_hash * colours_used, sizeof(Edges));
  db_graph.col_covgs = calloc(kmers_in_hash * colours_used, sizeof(Covg));

  // Print mem usage
  message("[memory]  graph: %s\n", graph_mem_str);
  hash_table_print_stats(&db_graph.ht);

  // Parse arguments, load
  SeqLoadingStats *stats = seq_loading_stats_create(1000);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1, .empty_colours = false,
                           .update_ginfo = true, .db_graph = &db_graph};

  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

  for(argi = 0; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--quality_score_threshold") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.quality_cutoff = tmp;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_ascii") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.ascii_fq_offset = tmp;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--remove_pcr_duplicates") == 0) {
      prefs.remove_dups_pe = true;
    }
    else if(strcmp(argv[argi],"--cut_homopolymers") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.homopolymer_cutoff = tmp;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--se_list") == 0) {
      parse_filelists(argv[argi+1], NULL, READ_FALIST, &prefs, stats,
                      seq_load_into_db_graph, NULL);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--pe_list") == 0) {
      parse_filelists(argv[argi+1], argv[argi+2], READ_FALIST, &prefs, stats,
                      seq_load_into_db_graph, NULL);
      argi += 2;
    }
    else if(strcmp(argv[argi],"--seq") == 0) {
      seq_parse_se(argv[argi+1], &r1, &r2, &prefs, stats,
                   seq_load_into_db_graph, NULL);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      seq_parse_pe(argv[argi+1], argv[argi+2], &r1, &r2, &prefs, stats,
                   seq_load_into_db_graph, NULL);
      argi += 2;
    }
    else if(strcmp(argv[argi],"--load_binary") == 0) {
      binary_load(argv[argi+1], &db_graph, &prefs, stats);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--colour_list") == 0) {
      parse_filelists(argv[argi+1], NULL, READ_COLOURLIST, &prefs, stats,
                      seq_load_into_db_graph, NULL);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--nc") == 0) {
      prefs.into_colour++;
    }
    else {
      die("Unknown command: %s", argv[argi]);
    }
  }

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);
  seq_loading_stats_free(stats);

  hash_table_print_stats(&db_graph.ht);

  binary_dump_graph(out_path, &db_graph, CURR_CTX_VERSION, NULL, 0, colours_used);

  message("Dumped cortex binary to: %s (format version %u; %u colour%s)\n",
          out_path, CURR_CTX_VERSION, colours_used, colours_used != 1 ? "s" : "");

  message("Done.\n");

  return EXIT_SUCCESS;
}
