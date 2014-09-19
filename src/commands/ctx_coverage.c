#include "global.h"
#include "commands.h"
#include "file_util.h"
#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "graph_format.h"
#include "graph_file_reader.h"

const char coverage_usage[] =
"usage: "CMD" coverage [options] <in.ctx> [in2.ctx ..]\n"
"\n"
"  Print contig coverage\n"
"\n"
"  -h, --help           This help message\n"
"  -q, --quiet          Silence status output normally printed to STDERR\n"
"  -f, --force          Overwrite output files\n"
"  -m, --memory <mem>   Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -e, --edges          Print edges as well. Uses hex encoding [TGCA|TGCA].\n"
"  -E, --degree         Print edge degree: 00!,01+,02{, 10-,11=,12<, 20},21>,22*\n"
"  -s, --seq <in>       Sequence file to get coverages for (can specify multiple times)\n"
"  -o, --out <out.txt>  Save output [default: STDOUT]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
// command specific
  {"edges",        no_argument,       NULL, 'e'},
  {"degrees",      no_argument,       NULL, 'E'},
  {"seq",          required_argument, NULL, '1'},
  {"seq",          required_argument, NULL, 's'},
  {NULL, 0, NULL, 0}
};

// Define a vector of Covg
#include "objbuf_macro.h"
create_objbuf(covg_buf,CovgBuffer,Covg);
create_objbuf(edges_buf,EdgesBuffer,Edges);

// [c]AGG[t]
// [a]CCT[g]
static inline void fetch_node_edges(const dBGraph *db_graph, dBNode node,
                                    Edges *dst)
{
  size_t i, ncols = db_graph->num_edge_cols;
  const Edges *edges = &db_node_edges(db_graph, node.key, 0);
  memcpy(dst, edges, ncols * sizeof(Edges));
  if(node.orient == REVERSE) {
    for(i = 0; i < ncols; i++) {
      // dst[i] = rev_nibble_lookup(dst[i]>>4) | (rev_nibble_lookup(dst[i]&0xf)<<4);
      dst[i] = (dst[i]>>4) | (dst[i]<<4);
    }
  }
}

// Print in/outdegree - For debugging mostly
// indegree/outdegree (2 means >=2)
// 00: ! 01: + 02: {
// 10: - 11: = 12: <
// 20: } 21: > 22: *
static inline
void _print_edge_degrees(const Edges *edges, size_t col, size_t ncols,
                         size_t num, FILE *fout)
{
  size_t i, indegree, outdegree;
  const char symbols[3][3] = {"!+{", "-=<", "}>*"};
  for(i = 0; i < num; i++) {
    indegree  = MIN2(edges_get_indegree(edges[i*ncols+col],  FORWARD), 2);
    outdegree = MIN2(edges_get_outdegree(edges[i*ncols+col], FORWARD), 2);
    fputc(symbols[indegree][outdegree], fout);
  }
  fputc('\n', fout);
}

static inline
void _print_edges(const Edges *edges, size_t col, size_t ncols,
                  size_t num, FILE *fout)
{
  size_t i;
  edges_print(fout, edges[col]);
  for(i = 1; i < num; i++) {
    fputc(' ', fout);
    edges_print(fout, edges[i*ncols+col]);
  }
  fputc('\n', fout);
}

static inline void print_read_covg(const dBGraph *db_graph, const read_t *r,
                                   CovgBuffer *covgbuf, EdgesBuffer *edgebuf,
                                   bool print_edges, bool print_edge_degrees,
                                   FILE *fout)
{
  // Find nodes, set covgs
  const size_t kmer_size = db_graph->kmer_size, ncols = db_graph->num_of_cols;
  size_t kmer_length = r->seq.end < kmer_size ? 0 : r->seq.end - kmer_size + 1;

  covg_buf_ensure_capacity(covgbuf, ncols * kmer_length);
  memset(covgbuf->data, 0, ncols * kmer_length * sizeof(Covg));

  if(db_graph->col_edges) {
    edges_buf_ensure_capacity(edgebuf, ncols * kmer_length);
    memset(edgebuf->data, 0, ncols * kmer_length * sizeof(Edges));
  }

  size_t i, j, col, search_start = 0;
  size_t contig_start, contig_end;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode node;
  Covg *covgs;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         0, 0)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size, 0, 0, &search_start);

    bkmer = binary_kmer_from_str(r->seq.b + contig_start, kmer_size);
    bkmer = binary_kmer_right_shift_one_base(bkmer);

    for(i = contig_start, j = contig_start+kmer_size-1; j < contig_end; j++, i++)
    {
      nuc = dna_char_to_nuc(r->seq.b[j]);
      bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
      node = db_graph_find(db_graph, bkmer);
      if(node.key != HASH_NOT_FOUND) {
        covgs = &db_node_covg(db_graph, node.key, 0);
        memcpy(covgbuf->data+i*ncols, covgs, ncols * sizeof(Covg));
        if(db_graph->col_edges) {
          fetch_node_edges(db_graph, node, edgebuf->data+i*ncols);
        }
      }
    }
  }

  // Print one colour per line
  fprintf(fout, ">%s\n%s\n", r->name.b, r->seq.b);
  if(kmer_length == 0) {
    for(i = 0; i < ncols; i++) {
      if(db_graph->col_edges) fputc('\n', fout);
      fputc('\n', fout);
    }
  }
  else {
    for(col = 0; col < ncols; col++)
    {
      if(print_edges)
        _print_edges(edgebuf->data, col, ncols, kmer_length, fout);

      if(print_edge_degrees)
        _print_edge_degrees(edgebuf->data, col, ncols, kmer_length, fout);

      // Print coverages
      fprintf(fout, "%2u", covgbuf->data[col]);
      for(i = 1; i < kmer_length; i++)
        fprintf(fout, " %2u", covgbuf->data[i*ncols+col]);
      fputc('\n', fout);
    }
  }
}

int ctx_coverage(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool print_edges = false, print_edge_degrees = false;
  const char *output_file = NULL;
  SeqFilePtrBuffer sfilebuf;

  seq_file_ptr_buf_alloc(&sfilebuf, 16);

  // tmp args
  size_t i;
  seq_file_t *tmp_sfile;

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!output_file, cmd); output_file = optarg; break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'e': cmd_check(!print_edges,cmd); print_edges = true; break;
      case 'E': cmd_check(!print_edge_degrees,cmd); print_edge_degrees = true; break;
      case '1':
      case 's':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&sfilebuf, tmp_sfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" coverage -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(sfilebuf.len == 0) cmd_print_usage("Require at least one --seq file");
  if(optind == argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  char **graph_paths = argv + optind;
  size_t num_gfiles = argc - optind;
  ctx_assert(num_gfiles > 0);
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // kmer memory = Edges + paths + 1 bit per colour
  bits_per_kmer = sizeof(GPath)*8 +
                  (sizeof(Covg) + print_edges*sizeof(Edges)) * 8 * ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        memargs.mem_to_use_set, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  FILE *fout = futil_open_create(output_file ? output_file : "-", "w");

  //
  // Set up memory
  //
  size_t kmer_size = gfiles[0].hdr.kmer_size;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, print_edges*ncols, kmers_in_hash,
                 DBG_ALLOC_COVGS | (print_edges ? DBG_ALLOC_EDGES : 0));

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  //
  // Load sequence
  //
  CovgBuffer covgbuf;
  EdgesBuffer edgebuf;
  covg_buf_alloc(&covgbuf, 2048);
  edges_buf_alloc(&edgebuf, 2048);

  read_t r;
  seq_read_alloc(&r);

  // Deal with one read at a time
  for(i = 0; i < sfilebuf.len; i++) {
    while(seq_read(sfilebuf.data[i], &r) > 0) {
      print_read_covg(&db_graph, &r, &covgbuf, &edgebuf,
                      print_edges, print_edge_degrees, fout);
    }
    seq_close(sfilebuf.data[i]);
  }

  seq_read_dealloc(&r);
  covg_buf_dealloc(&covgbuf);
  edges_buf_dealloc(&edgebuf);

  seq_file_ptr_buf_dealloc(&sfilebuf);

  fclose(fout);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
