#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "vcf_parsing.h"

const char ctx_usage[] =
"usage: "CMD" covg [options] <in.ctx> <input.vcf> <out.vcf>\n"
"  Print coverage on some input file\n";

static int parse_vcf_header(StrBuf *line, gzFile vcf, bool print,
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
  bool found;
  size_t kmer_size = db_graph->kmer_size;

  bkmer = binary_kmer_from_str(contig, kmer_size);
  bkmer = binary_kmer_right_shift_one_base(bkmer);

  for(next_base = kmer_size-1; next_base < contig_len; next_base++)
  {
    nuc = dna_char_to_nuc(contig[next_base]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    bkey = bkmer_get_key(bkmer, kmer_size);
    hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  }
}

int ctx_covg(CmdArgs *args)
{
  cmd_accept_options(args, "mh", usage);
  cmd_require_options(args, "mh", usage);
  int argc = args->argc;
  char **argv = args->argv;

  if(argc < 3) cmd_print_usage(NULL);

  // Parse commandline args
  char *in_ctx_path, *in_vcf_path, *out_vcf_path;

  in_ctx_path = argv[0];
  in_vcf_path = argv[1];
  out_vcf_path = argv[2];

  // Probe binary
  bool is_binary = false;
  GraphFileHeader gheader = INIT_GRAPH_FILE_HDR;

  if(!graph_file_probe(in_ctx_path, &is_binary, &gheader))
    cmd_print_usage("Cannot read input binary file: %s", in_ctx_path);
  else if(!is_binary)
    cmd_print_usage("Input binary file isn't valid: %s", in_ctx_path);

  size_t kmer_size = gheader.kmer_size;

  gzFile vcf;
  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  size_t i, num_col_given = argc - 5;
  uint32_t colours[num_col_given];

  for(i = 0; i < num_col_given; i++) {
    if(!parse_entire_uint(argv[5+i], colours+i))
      cmd_print_usage("Invalid colour number: %s", argv[5+i]);
  }

  if(!futil_is_file_writable(out_vcf_path))
    cmd_print_usage("Cannot write to output file: %s", out_vcf_path);

  if(num_col_given > 0 && gheader.num_of_cols != num_col_given)
    die("You're using more colours than you need!");

  size_t cols_used = num_col_given > 0 ? num_col_given : gheader.num_of_cols;

  // Figure out how much memory to use
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  bits_per_kmer = (sizeof(Edges) + sizeof(Covg)*cols_used) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gheader.num_of_kmers, true, &graph_mem);

  cmd_check_mem_limit(args->mem_to_use, graph_mem);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, cols_used, 1, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = ctx_calloc(db_graph.ht.capacity * cols_used, sizeof(Covg));

  // Print mem usage
  hash_table_print_stats(&db_graph.ht);

  // Load sequence from VCF to build hash table
  StrBuf line, contig;
  strbuf_alloc(&line, 1024);
  strbuf_alloc(&contig, 1024);

  // Read header
  size_t num_samples = parse_vcf_header(&line, vcf, false, in_vcf_path);

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

  // Load graph only for kmers already in the hash table

  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs prefs = {.into_colour = 0, .db_graph = &db_graph,
                             .boolean_covgs = false,
                             .must_exist_in_graph = true,
                             .empty_colours = false};

  graph_load(in_ctx_path, prefs, &stats, NULL);

  if((vcf = gzopen(in_vcf_path, "r")) == NULL)
    die("Couldn't open file: %s", in_vcf_path);

  // Read (and print) header
  parse_vcf_header(&line, vcf, true, in_vcf_path);

  size_t new_samples = num_samples + cols_used;

  size_t j;

  size_t nodes_cap = 2048, nodes_len = 0;
  dBNode *nodes = ctx_malloc(sizeof(dBNode) * nodes_cap);

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
        nodes_cap = roundup2pow(nodes_cap);
        nodes = ctx_realloc(nodes, sizeof(dBNode) * nodes_cap);
      }

      // DEV check contig length vs kmer_size

      BinaryKmer bkmer;
      Nucleotide nuc;

      bkmer = binary_kmer_from_str(contig.buff, kmer_size);
      nodes[0].key = hash_table_find(&db_graph.ht, bkmer);
      size_t k, num_nodes = contig.len+1-kmer_size;

      for(j = 1; j <= num_nodes; j++)
      {
        nuc = dna_char_to_nuc(contig.buff[j+kmer_size-1]);
        bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
        nodes[j].key = hash_table_find(&db_graph.ht, bkmer);
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

  ctx_free(covg_array);

  graph_header_dealloc(&gheader);

  ctx_free(db_graph.col_edges);
  ctx_free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
