#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
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
  BinaryKmer bkmer, bkey;
  Nucleotide nuc;
  size_t next_base;
  boolean found;
  uint32_t kmer_size = db_graph->kmer_size;

  bkmer = binary_kmer_from_str(contig, kmer_size);
  binary_kmer_right_shift_one_base(&bkmer);

  for(next_base = kmer_size-1; next_base < contig_len; next_base++)
  {
    nuc = binary_nuc_from_char(contig[next_base]);
    binary_kmer_left_shift_add(&bkmer, kmer_size, nuc);
    bkey = db_node_get_key(bkmer, kmer_size);
    hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  }
}

int ctx_covg(CmdArgs *args)
{
  cmd_accept_options(args, "mh", usage);
  cmd_require_options(args, "mh", usage);
  int argc = args->argc;
  char **argv = args->argv;

  if(argc < 3) print_usage(usage, NULL);

  // Parse commandline args
  char *in_ctx_path, *in_vcf_path, *out_vcf_path;

  in_ctx_path = argv[0];
  in_vcf_path = argv[1];
  out_vcf_path = argv[2];

  // Probe binary
  boolean is_binary = false;
  GraphFileHeader gheader = INIT_GRAPH_FILE_HDR;

  if(!graph_file_probe(in_ctx_path, &is_binary, &gheader))
    print_usage(usage, "Cannot read input binary file: %s", in_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", in_ctx_path);

  uint32_t kmer_size = gheader.kmer_size;

  gzFile vcf;
  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  size_t i, num_col_given = argc - 5;
  uint32_t colours[num_col_given];

  for(i = 0; i < num_col_given; i++) {
    if(!parse_entire_uint(argv[5+i], colours+i))
      print_usage(usage, "Invalid colour number: %s", argv[5+i]);
  }

  if(!test_file_writable(out_vcf_path))
    print_usage(usage, "Cannot write to output file: %s", out_vcf_path);

  if(num_col_given > 0 && gheader.num_of_cols != num_col_given)
    die("You're using more colours than you need!");

  uint32_t cols_used = num_col_given > 0 ? num_col_given : gheader.num_of_cols;

  // Figure out how much memory to use
  size_t bits_per_kmer, kmers_in_hash;
  bits_per_kmer = (sizeof(Edges) + sizeof(Covg)*cols_used) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gheader.num_of_kmers, true);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, cols_used, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity * cols_used, sizeof(Covg));

  // Print mem usage
  hash_table_print_stats(&db_graph.ht);

  // Load sequence from VCF to build hash table
  StrBuf line, contig;
  strbuf_alloc(&line, 1024);
  strbuf_alloc(&contig, 1024);

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
      strbuf_append_buff(&contig, &vcf_entry.alts[i]);
      strbuf_append_buff(&contig, vcf_entry.rf);

      // Add kmers to the graph
      add_str_to_graph(&db_graph, contig.buff, contig.len);
    }
  }

  gzclose(vcf);

  // Load binary only for kmers already in the hash table

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .db_graph = &db_graph,
                           .boolean_covgs = false,
                           .must_exist_in_graph = true,
                           .empty_colours = false};

  graph_load(in_ctx_path, &prefs, stats, NULL);

  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  // Read (and print) header
  parse_vcf_header(&line, vcf, true, in_vcf_path);

  size_t new_samples = num_samples + cols_used;

  size_t j;

  size_t nodes_cap = 2048, nodes_len = 0;
  hkey_t *nodes = malloc2(sizeof(hkey_t) * nodes_cap);

  DeltaArray delta;
  delta_arr_alloc(&delta);

  // PARSE and PRINT
  while(strbuf_gzreadline_nonempty(&line, vcf))
  {
    strbuf_chomp(&line);
    vcf_entry_parse(&line, &vcf_entry, num_samples);

    strbuf_reset(&contig);
    strbuf_append_buff(&contig, vcf_entry.lf);

    for(i = 0; i < vcf_entry.num_alts; i++)
    {
      strbuf_shrink(&contig, vcf_entry.lf->len);
      strbuf_append_buff(&contig, &vcf_entry.alts[i]);
      strbuf_append_buff(&contig, vcf_entry.rf);

      if(contig.len > nodes_cap) {
        nodes_cap = ROUNDUP2POW(nodes_cap);
        nodes = realloc2(nodes, sizeof(hkey_t) * nodes_cap);
      }

      // DEV check contig length vs kmer_size

      BinaryKmer bkmer;
      Nucleotide nuc;

      bkmer = binary_kmer_from_str(contig.buff, kmer_size);
      nodes[0] = hash_table_find(&db_graph.ht, bkmer);
      uint32_t k, num_nodes = contig.len+1-kmer_size;

      for(j = 1; j <= num_nodes; j++)
      {
        nuc = binary_nuc_from_char(contig.buff[j+kmer_size-1]);
        binary_kmer_left_shift_add(&bkmer, kmer_size, nuc);
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

  graph_header_dealloc(&gheader);

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
