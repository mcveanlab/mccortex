#include "global.h"

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"
#include "vcf_parsing.h"

static const char usage[] =
"usage: ctx_covg <mem> <in.ctx> <input.vcf> <out.vcf> [col1 ...]\n"
"  Print coverage on some input file\n";

int main(int argc, char **argv)
{
  if(argc < 5) print_usage(usage, NULL);

  // Parse commandline args
  size_t mem_to_use;
  char *in_ctx_path, *in_vcf_path, *out_vcf_path;

  if(!mem_to_integer(argv[1], &mem_to_use))
    print_usage(usage, "Invalid <mem> arg (try 1GB or 2M): %s", argv[1]);

  in_ctx_path = argv[2];
  in_vcf_path = argv[3];
  out_vcf_path = argv[4];

  // Probe binary
  boolean is_binary = false;
  uint32_t kmer_size, ctx_num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(in_ctx_path, &is_binary, &kmer_size, &ctx_num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read input binary file: %s", in_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", in_ctx_path);

  if(!test_file_readable(in_vcf_path))
    print_usage(usage, "Cannot read input VCF: %s", in_vcf_path);

  size_t i, num_col_given = argc - 5;
  uint32_t colours[num_col_given];

  for(i = 0; i < num_col_given; i++) {
    if(!parse_entire_uint(argv[5+i], colours+i))
      print_usage(usage, "Invalid colour number: %s", argv[5+i]);
  }

  if(!test_file_writable(out_vcf_path))
    print_usage(usage, "Cannot write to output file: %s", out_vcf_path);

  uint32_t cols_used = num_col_given > 0 ? num_col_given : ctx_num_of_cols;

  // Figure out how much mem to use
  size_t mem_per_kmer, kmers_in_hash, hash_mem, graph_mem;

  mem_per_kmer = sizeof(BinaryKmer) + sizeof(Edges) + sizeof(Covg)*cols_used;
  hash_mem = hash_table_mem2(mem_to_use / mem_per_kmer, &kmers_in_hash);
  graph_mem = kmers_in_hash * mem_per_kmer;

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, cols_used, kmers_in_hash);
  db_graph.edges = calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc(db_graph.ht.capacity * cols_used, sizeof(Covg));

  // Print mem usage
  char graph_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  message("[memory]  graph: %s\n", graph_mem_str);
  hash_table_print_stats(&db_graph.ht);

  message("Using kmer size %u with %u colours\n", kmer_size, cols_used);

  // Load sequence from VCF to build hash table
  gzFile vcf;

  StrBuf line;
  strbuf_alloc(&line, 1024);

  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  // DEV: parse header
  uint32_t num_samples = 0;

  vcf_entry_t vcf_entry;
  vcf_entry_alloc(&vcf_entry, num_samples);

  while(strbuf_gzreadline_nonempty(&line, vcf))
  {
    strbuf_chomp(&line);
    vcf_entry_parse(&line, &vcf_entry, num_samples);


  }

  gzclose(vcf);

  // Load binary only for kmers already in the hash table

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = 0,
                           .empty_colours = false,
                           .update_ginfo = true,
                           .db_graph = &db_graph};

  binary_load(in_ctx_path, &db_graph, &prefs, stats);

  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  // DEV: PARSE and PRINT
  while(strbuf_gzreadline_nonempty(&line, vcf))
  {
    strbuf_chomp(&line);
    vcf_entry_parse(&line, &vcf_entry, num_samples);

    
  }

  gzclose(vcf);

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
