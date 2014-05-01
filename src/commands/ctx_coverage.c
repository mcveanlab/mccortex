#include "global.h"
#include "commands.h"

#include "file_util.h"
#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "graph_format.h"
#include "graph_file_filter.h"

#include "seq_file.h"

const char coverage_usage[] =
"usage: "CMD" coverage [options] <in.ctx> [in2.ctx ..]\n"
"  Print contig coverage\n"
"\n"
"  --edges          Print edges as well\n"
"  --seq <in>       Trusted input (can be used multiple times)\n"
"  --out <out.txt>  Save calls [default: STDOUT]\n";

// Define a vector of Covg
#include "objbuf_macro.h"
create_objbuf(covg_buf,CovgBuffer,Covg);
create_objbuf(edges_buf,EdgesBuffer,Edges);

static inline void fetch_node_edges(const dBGraph *db_graph, dBNode node,
                                    Edges *dst)
{
  size_t i, ncols = db_graph->num_edge_cols;
  const Edges *edges = &db_node_edges(db_graph, node.key, 0);
  memcpy(dst, edges, ncols * sizeof(Edges));
  if(node.orient == REVERSE) {
    for(i = 0; i < ncols; i++) {
      dst[i] = (dst[i]>>4) | (dst[i]<<4);
    }
  }
}

static inline void print_read_covg(const dBGraph *db_graph, const read_t *r,
                                   CovgBuffer *covgbuf, EdgesBuffer *edgebuf,
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
        if(db_graph->col_edges)
          fetch_node_edges(db_graph, node, edgebuf->data+i*ncols);
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
      if(db_graph->col_edges) {
        // Print edges
        edges_print(fout, edgebuf->data[col]);
        for(i = 1; i < kmer_length; i++) {
          fputc(' ', fout);
          edges_print(fout, edgebuf->data[i*ncols+col]);
        }
        fputc('\n', fout);
      }
      // Print coverages
      fprintf(fout, "%2u", covgbuf->data[col]);
      for(i = 1; i < kmer_length; i++)
        fprintf(fout, " %2u", covgbuf->data[i*ncols+col]);
      fputc('\n', fout);
    }
  }
}

int ctx_coverage(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;
  // Already checked that input has at least 3 arguments

  size_t max_seq_inputs = argc / 2;
  seq_file_t *seq_files[max_seq_inputs];
  size_t i, num_seq_files = 0;
  bool print_edges = false;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--seq") == 0) {
      if(argi+1 == argc) cmd_print_usage("--seq <in.fa> requires a file");
      if((seq_files[num_seq_files] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file %s", argv[argi+1]);
      num_seq_files++;
      argi++; // We took an argument
    }
    else if(strcmp(argv[argi],"--edges") == 0) print_edges = true;
    else {
      cmd_print_usage("Unknown arg: %s", argv[argi]);
    }
  }

  if(num_seq_files == 0) cmd_print_usage("Require at least one --seq file");
  if(argi == argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  char **graph_paths = argv + argi;
  size_t num_gfiles = argc - argi;
  GraphFileReader gfiles[num_gfiles];
  size_t ncols = 0, ctx_max_cols = 0, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers, &ctx_max_cols);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // kmer memory = Edges + paths + 1 bit per colour
  bits_per_kmer = (sizeof(Covg) + print_edges*sizeof(Edges)) * 8 * ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        args->mem_to_use_set, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  //
  // Open output file
  //
  FILE *fout;

  if(!args->output_file_set || strcmp("-", args->output_file) == 0)
    fout = stdout;
  else
    fout = fopen(args->output_file, "w");

  if(fout == NULL) die("Cannot open output file: %s", args->output_file);

  //
  // Set up memory
  //
  size_t kmer_size = gfiles[0].hdr.kmer_size;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, print_edges*ncols, kmers_in_hash);
  db_graph.col_covgs = ctx_calloc(db_graph.ht.capacity*ncols, sizeof(Covg));

  if(print_edges)
    db_graph.col_edges = ctx_calloc(db_graph.ht.capacity*ncols, sizeof(Edges));

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++)
    graph_load(&gfiles[i], gprefs, &stats);

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
  for(i = 0; i < num_seq_files; i++) {
    while(seq_read(seq_files[i], &r) > 0) {
      print_read_covg(&db_graph, &r, &covgbuf, &edgebuf, fout);
    }
    seq_close(seq_files[i]);
  }

  seq_read_dealloc(&r);
  covg_buf_dealloc(&covgbuf);
  edges_buf_dealloc(&edgebuf);

  // Finished: do clean up
  if(fout != stdout) fclose(fout);

  if(print_edges) ctx_free(db_graph.col_edges);
  ctx_free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
