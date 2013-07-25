#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_format.h"
#include "vcf_parsing.h"

static const char usage[] =
"usage: "CMD" covg [options] <in.ctx> <input.vcf> <out.vcf>\n"
"  Print coverage on some input file\n";

static int parse_vcf_header(StrBuf *line, gzFile vcf, boolean print,
                            const char *in_vcf_path)
{
  while(strbuf_gzreadline_nonempty(line, vcf) > 0)
  {
    if(print) puts(line->buff);
    if(strncmp(line->buff,"#CHROM",6) == 0) {
      size_t columns = count_char(line->buff, '\t') + 1;
      if(columns < VCFSAMPLES) die("Too few columns in VCF: %s", in_vcf_path);
      return columns - VCFSAMPLES;
    }
    else if(strncmp(line->buff,"##",2) != 0) break;
  }
  die("Missing VCF header: %s", in_vcf_path);
}

void add_str_to_graph(dBGraph *db_graph, const char *contig, size_t contig_len)
{
  BinaryKmer bkmer, tmp_key;
  Nucleotide nuc;
  size_t next_base;
  boolean found;
  uint32_t kmer_size = db_graph->kmer_size;

  binary_kmer_from_str(contig, kmer_size, bkmer);
  binary_kmer_right_shift_one_base(bkmer);

  for(next_base = kmer_size-1; next_base < contig_len; next_base++)
  {
    nuc = binary_nuc_from_char(contig[next_base]);
    binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    db_node_get_key(bkmer, kmer_size, tmp_key);
    hash_table_find_or_insert(&db_graph->ht, tmp_key, &found);
  }
}

int ctx_covg(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  cmd_require_options(args, "m");
  int argc = args->argc;
  char **argv = args->argv;

  size_t mem_to_use = args->mem_to_use;

  if(argc < 3) print_usage(usage, NULL);

  // Parse commandline args
  char *in_ctx_path, *in_vcf_path, *out_vcf_path;

  in_ctx_path = argv[0];
  in_vcf_path = argv[1];
  out_vcf_path = argv[2];

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

  if(num_col_given > 0 && ctx_num_of_cols != num_col_given)
    die("You're using more colours than you need!");

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

  StrBuf line, contig;
  strbuf_alloc(&line, 1024);
  strbuf_alloc(&contig, 1024);

  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  // Read header
  uint32_t num_samples = parse_vcf_header(&line, vcf, false, in_vcf_path);

  vcf_entry_t vcf_entry;
  vcf_entry_alloc(&vcf_entry, num_samples);

  while(strbuf_gzreadline_nonempty(&line, vcf))
  {
    strbuf_chomp(&line);
    vcf_entry_parse(&line, &vcf_entry, num_samples);

    strbuf_reset(&contig);
    strbuf_append_buff(&contig, vcf_entry.lf);

    for(i = 0; i < vcf_entry.num_alts; i++)
    {
      strbuf_shrink(&contig, vcf_entry.lf->len);
      strbuf_append_buff(&contig, vcf_entry.alts[i]);
      strbuf_append_buff(&contig, vcf_entry.rf);

      // Add kmers to the graph
      add_str_to_graph(&db_graph, contig.buff, contig.len);
    }
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

  // Read (and print) header
  parse_vcf_header(&line, vcf, true, in_vcf_path);

  size_t new_samples = num_samples + cols_used;

  size_t j, capacity = 16 * num_samples;
  DeltaArray *covg_array = malloc(sizeof(DeltaArray) * capacity);
  for(i = 0; i < capacity; i++) delta_arr_alloc(&covg_array[i]);

  size_t nodes_cap = 2048, nodes_len = 0;
  hkey_t *nodes = malloc(sizeof(hkey_t) * nodes_cap);

  DeltaArray delta;
  delta_arr_alloc(&delta);

  // PARSE and PRINT
  while(strbuf_gzreadline_nonempty(&line, vcf))
  {
    strbuf_chomp(&line);
    vcf_entry_parse(&line, &vcf_entry, num_samples);

    strbuf_reset(&contig);
    strbuf_append_buff(&contig, vcf_entry.lf);

    if(vcf_entry.num_alts*num_samples > capacity) {
      size_t new_cap = ROUNDUP2POW(vcf_entry.num_alts*num_samples);
      covg_array = realloc(covg_array, sizeof(DeltaArray) * new_cap);
      for(i = capacity; i < new_cap; i++) { delta_arr_alloc(&covg_array[i]); }
      capacity = new_cap;
    }

    for(i = 0; i < vcf_entry.num_alts; i++)
    {
      strbuf_shrink(&contig, vcf_entry.lf->len);
      strbuf_append_buff(&contig, vcf_entry.alts[i]);
      strbuf_append_buff(&contig, vcf_entry.rf);

      if(contig.len > nodes_cap) {
        nodes_cap = ROUNDUP2POW(nodes_cap);
        nodes = realloc(nodes, sizeof(hkey_t) * nodes_cap);
      }

      // DEV check contig length vs kmer_size

      BinaryKmer bkmer;
      Nucleotide nuc;

      binary_kmer_from_str(contig.buff, kmer_size, bkmer);
      nodes[0] = hash_table_find(&db_graph.ht, bkmer);
      uint32_t k, num_nodes = contig.len+1-kmer_size;

      for(j = 1; j <= num_nodes; j++)
      {
        nuc = binary_nuc_from_char(contig.buff[j+kmer_size-1]);
        binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
        nodes[j] = hash_table_find(&db_graph.ht, bkmer);
      }

      // Get contig covg
      for(j = 0; j < num_samples; j++) {
        for(k = 0; k < num_nodes; k++) {

        }
      }
    }

    // Print

  }

  gzclose(vcf);

  for(i = 0; i < num_samples; i++)
      delta_arr_dealloc(&covg_array[i]);

  free(covg_array);

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
