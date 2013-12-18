#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "graph_format.h"
#include "path_format.h"
#include "supernode.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"  Pull out contigs from the graph, print statistics\n"
"\n"
"  Options: [ -n <nkmers> | -p <paths.ctp> ]\n"
"   --memory <M>     How much memory to use\n"
"   --colour <c>     Pull out contigs from the given colour [default: 0]\n"
"   --ncontigs <N>   Pull out <N> contigs from random kmers\n"
"   --print          Print contigs in FASTA format\n"
"   --out <out.fa>   Write contigs to a file rather than STDOUT\n"
"   --seed <in.fa>   Use seed kmers from a file. If longer than kmer-size, only\n"
"                    use the first kmer found from each input sequence.\n";

#define MAXPATH 5

typedef struct {
  size_t ncontigs, capacity;
  size_t total_len, total_junc;
  size_t contigs_outdegree[5];
  size_t paths_held[MAXPATH], paths_pickdup[MAXPATH], paths_counter[MAXPATH];
  size_t *lengths, *junctions;
  size_t max_len, max_junc;
  double max_junc_density;
  dBNodeBuffer nodes;
  size_t nprint; // id of next contig printed
} ContigData;

static void contig_print_path_dist(const size_t *hist, size_t n,
                                   const char *name, size_t ncontigs)
{
  char nout_str[90];
  size_t i;

  timestamp(ctx_msg_out);
  message(" %s: ", name);
  for(i = 0; i < n; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(hist[i], nout_str),
            (size_t)((100.0*hist[i])/(2.0*ncontigs)+0.5));
  }
  message("\n");
}

static inline void contig_data_alloc(ContigData *cd, size_t capacity)
{
  size_t *lengths = malloc2(capacity * sizeof(size_t));
  size_t *junctions = malloc2(capacity * sizeof(size_t));
  ContigData tmp = {.ncontigs = 0, .capacity = capacity,
                    .total_len = 0, .total_junc = 0,
                    .contigs_outdegree = {0}, .paths_held = {0},
                    .paths_pickdup = {0}, .paths_counter = {0},
                    .lengths = lengths, .junctions = junctions,
                    .max_len = 0, .max_junc = 0,
                    .max_junc_density = 0, .nprint = 0};
  db_node_buf_alloc(&tmp.nodes, 1024);
  memcpy(cd, &tmp, sizeof(ContigData));
}

static inline void contig_data_ensure_capacity(ContigData *cd, size_t capacity)
{
  if(capacity > cd->capacity) {
    cd->capacity = ROUNDUP2POW(capacity);
    cd->lengths = realloc2(cd->lengths, cd->capacity * sizeof(size_t));
    cd->junctions = realloc2(cd->junctions, cd->capacity * sizeof(size_t));
  }
}

static inline void contig_data_dealloc(ContigData *cd)
{
  free(cd->lengths);
  free(cd->junctions);
  db_node_buf_dealloc(&cd->nodes);
}

static void pulldown_contig(hkey_t node, ContigData *cd,
                            const dBGraph *db_graph, size_t colour,
                            GraphWalker *wlk, RepeatWalker *rptwlk,
                            FILE *fout)
{
  Orientation orient;
  Nucleotide lost_nuc;

  const size_t kmer_size = db_graph->kmer_size;
  size_t njunc = 0;

  db_nodes_reset(&cd->nodes);
  db_node_buf_safe_add(&cd->nodes, node, FORWARD);

  for(orient = 0; orient < 2; orient++)
  {
    if(orient == 1) {
      supernode_reverse(cd->nodes.data, cd->nodes.len);
      node = cd->nodes.data[cd->nodes.len-1].key;
    }

    graph_walker_init(wlk, db_graph, colour, colour, node, orient);
    lost_nuc = binary_kmer_first_nuc(wlk->bkmer, kmer_size);

    while(graph_traverse(wlk) &&
          walker_attempt_traverse(rptwlk, wlk, wlk->node, wlk->orient, wlk->bkmer))
    {
      graph_walker_node_add_counter_paths(wlk, lost_nuc);
      lost_nuc = binary_kmer_first_nuc(wlk->bkmer, kmer_size);
      db_node_buf_safe_add(&cd->nodes, wlk->node, wlk->orient);
    }

    njunc += wlk->fork_count;
    cd->paths_held[MIN2(wlk->num_curr, MAXPATH-1)]++;
    cd->paths_pickdup[MIN2(wlk->num_new, MAXPATH-1)]++;
    cd->paths_counter[MIN2(wlk->num_counter, MAXPATH-1)]++;
    // nloop += rptwlk->nbloom_entries;

    graph_walker_finish(wlk);
    walker_fast_clear(rptwlk, cd->nodes.data, cd->nodes.len);
  }

  if(fout != NULL) {
    fprintf(fout, ">contig%zu\n", cd->nprint);
    db_nodes_print(cd->nodes.data, cd->nodes.len, db_graph, fout);
    putc('\n', fout);
    cd->nprint++;
  }

  // Out degree
  size_t outdegree, len = cd->nodes.len;
  hkey_t firstnode = cd->nodes.data[0].key;
  hkey_t firstorient = opposite_orientation(cd->nodes.data[0].orient);
  hkey_t lastnode = cd->nodes.data[len-1].key;
  hkey_t lastorient = cd->nodes.data[len-1].orient;

  outdegree = edges_get_outdegree(db_graph->col_edges[firstnode], firstorient);
  cd->contigs_outdegree[outdegree]++;
  outdegree = edges_get_outdegree(db_graph->col_edges[lastnode], lastorient);
  cd->contigs_outdegree[outdegree]++;

  cd->lengths[cd->ncontigs] = len;
  cd->junctions[cd->ncontigs] = njunc;
  cd->total_len += len;
  cd->total_junc += njunc;

  cd->max_junc_density = MAX2(cd->max_junc_density, (double)njunc / len);
  cd->max_len = MAX2(cd->max_len, len);
  cd->max_junc = MAX2(cd->max_junc, njunc);

  cd->ncontigs++;
}

struct ParseSeeds {
  ContigData *const cd;
  const dBGraph *const db_graph;
  const size_t colour;
  GraphWalker *const wlk;
  RepeatWalker *const rptwlk;
  FILE *const fout;
};

// Pull down a contig from the first node that is found in the graph
static inline void parse_seed_read(const read_t *r, struct ParseSeeds *ps)
{
  const dBGraph *db_graph = ps->db_graph;

  if(r == NULL) return;
  hkey_t node = seq_reader_first_node(r, -1, -1, ps->colour, db_graph);

  if(node != HASH_NOT_FOUND) {
    pulldown_contig(node, ps->cd, db_graph, ps->colour,
                    ps->wlk, ps->rptwlk, ps->fout);
  }
  else if(ps->fout != NULL) {
    fprintf(ps->fout, ">contig%zu\n\n", ps->cd->nprint);
    ps->cd->nprint++;
  }
}

static inline void parse_seed_reads(read_t *r1, read_t *r2,
                                    int qoffset1, int qoffset2,
                                    const SeqLoadingPrefs *prefs,
                                    SeqLoadingStats *stats,
                                    void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  (void)prefs; (void)stats;

  struct ParseSeeds *ps = (struct ParseSeeds*)ptr;
  contig_data_ensure_capacity(ps->cd, ps->cd->ncontigs+2);
  parse_seed_read(r1, ps);
  parse_seed_read(r2, ps);
}

int ctx_contigs(CmdArgs *args)
{
  cmd_accept_options(args, "mnpo", usage);
  int argc = args->argc;
  char **argv = args->argv;

  // char str[100];
  // printf("Test NaN: %s INF: %s\n", num_to_str(NAN, 2, str),
  //                                  num_to_str(INFINITY, 2, str));

  size_t i, n_rand_contigs = 0, colour = 0;
  boolean print_contigs = false;
  seq_file_t *seed_file = NULL;

  while(argc > 0 && argv[0][0] == '-') {
    if(strcmp(argv[0],"--ncontigs") == 0) {
      if(argc == 1 || !parse_entire_size(argv[1], &n_rand_contigs) ||
         n_rand_contigs == 0) {
        print_usage(usage, "--ncontigs <N> requires an integer argument [>0]");
      }
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--colour") == 0 || strcmp(argv[0],"--color") == 0) {
      if(argc == 1 || !parse_entire_size(argv[1], &colour))
        print_usage(usage, "--colour <c> requires a positive integer argument");
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--seed") == 0 || strcmp(argv[0],"--seeds") == 0) {
      if(argc == 1)
        print_usage(usage, "--seed <in.fa|fq|sam> requires an argument");
      if(seed_file != NULL) {
        print_usage(usage, "Only one --seed allowed, "
                           "try --seed <(cat file1.fa file2.fa ...)");
      }
      if((seed_file = seq_open(argv[1])) == NULL)
        die("Cannot read --seed file: %s", argv[1]);
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--print") == 0) {
      print_contigs = true;
      argv++; argc--;
    }
    else print_usage(usage, "Unknown argument: %s", argv[0]);
  }

  if(argc != 1) print_usage(usage, NULL);
  char *input_ctx_path = argv[0];

  if(seed_file != NULL && n_rand_contigs > 0)
    print_usage(usage, "Please specify one of --seed and --ncontigs");

  if(seed_file == NULL && n_rand_contigs == 0) n_rand_contigs = 1000;

  // Seed random
  seed_random();

  //
  // Open graph file
  //
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open(&file, input_ctx_path, true);

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t path_max_mem = 0, path_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    path_max_mem = MAX2(path_max_mem, pfiles[i].hdr.num_path_bytes);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check colour specified
  if(colour >= graph_file_usedcols(&file)) {
    print_usage(usage, "--colour is too high (%zu > %zu)",
                colour, graph_file_usedcols(&file)-1);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;
  char path_mem_str[100], total_mem_str[100];

  bits_per_kmer = sizeof(Edges)*8 + file.hdr.num_of_cols + sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        file.hdr.num_of_kmers, false);

  graph_mem = hash_table_mem(kmers_in_hash, false, NULL) +
              (bits_per_kmer * kmers_in_hash) / 8;


  // Paths memory
  size_t tmppathsize = paths_merge_needs_tmp(pfiles, num_pfiles) ? path_max_mem : 0;
  path_mem = path_max_mem + tmppathsize;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s", path_mem_str);

  // Total memory
  total_mem = graph_mem + path_mem;
  bytes_to_str(total_mem, 1, total_mem_str);

  if(total_mem > args->mem_to_use)
    die("Requires at least %s memory", total_mem_str);

  //
  // Output file if printing
  //
  FILE *fout = print_contigs ? stdout : NULL;
  if(args->output_file_set) {
    if(!print_contigs)
      warn("Ignoring --out <out> argument (maybe you forgot --print ?)");
    else if((fout = fopen(args->output_file, "r")) == NULL)
      die("Cannot open output file: %s", args->output_file);
  }

  dBGraph db_graph;
  GraphWalker wlk;

  db_graph_alloc(&db_graph, file.hdr.kmer_size, file.hdr.num_of_cols, 1, kmers_in_hash);
  graph_walker_alloc(&wlk);

  size_t node_bit_fields = round_bits_to_words64(db_graph.ht.capacity);

  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = calloc2(node_bit_fields * file.hdr.num_of_cols,
                                  sizeof(uint64_t));
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  RepeatWalker rptwlk;
  walker_alloc(&rptwlk, db_graph.ht.capacity, 22); // 4MB

  uint8_t *path_store = malloc2(path_max_mem);
  path_store_init(&db_graph.pdata, path_store, path_max_mem, path_max_usedcols);

  uint8_t *tmppdata = tmppathsize > 0 ? malloc2(tmppathsize) : NULL;

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           // Graph files
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = true,
                           // Sequence files
                           .into_colour = 0,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false};

  graph_load(&file, &prefs, stats);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  paths_format_merge(pfiles, num_pfiles, false, tmppdata, tmppathsize, &db_graph);

  status("Traversing graph in colour %zu...", colour);

  hkey_t node;
  ContigData cd;

  if(seed_file == NULL)
  {
    contig_data_alloc(&cd, n_rand_contigs);
    while(cd.ncontigs < n_rand_contigs) {
      do { node = db_graph_rand_node(&db_graph); }
      while(!db_node_has_col(&db_graph, node, colour));
      pulldown_contig(node, &cd, &db_graph, colour, &wlk, &rptwlk, fout);
    }
  }
  else
  {
    // Parse seed file
    struct ParseSeeds ps = {.cd = &cd, .db_graph = &db_graph, .colour = colour,
                            .wlk = &wlk, .rptwlk = &rptwlk, .fout = fout};

    contig_data_alloc(&cd, 1024);

    read_t r1, r2;
    seq_read_alloc(&r1);
    seq_read_alloc(&r2);
    seq_parse_se_sf(seed_file, &r1, &r2, &prefs, stats, parse_seed_reads, &ps);
    seq_read_dealloc(&r1);
    seq_read_dealloc(&r2);
    // seq_parse_se_sf() closes seed_file
  }

  seq_loading_stats_free(stats);

  if(args->output_file_set && print_contigs && strcmp(args->output_file, "-") != 0)
    fclose(fout);

  status("\n");

  qsort(cd.lengths, cd.ncontigs, sizeof(size_t), cmp_size);
  qsort(cd.junctions, cd.ncontigs, sizeof(size_t), cmp_size);

  double median_len = MEDIAN(cd.lengths, cd.ncontigs);
  double median_junc = MEDIAN(cd.junctions, cd.ncontigs);

  char ncontigs_str[90], nprint_str[90];
  num_to_str(cd.ncontigs, 1, ncontigs_str);
  num_to_str(cd.nprint, 1, nprint_str);

  status("pulled out %s contigs", ncontigs_str);
  status("printed out %s contigs", nprint_str);

  char len_mean_str[90], len_median_str[90], len_max_str[90], len_total_str[90];
  char jnc_mean_str[90], jnc_median_str[90], jnc_max_str[90], jnc_total_str[90];

  num_to_str((double)cd.total_len / cd.ncontigs, 1, len_mean_str);
  num_to_str((double)cd.total_junc / cd.ncontigs, 1, jnc_mean_str);
  num_to_str(median_len, 1, len_median_str);
  num_to_str(median_junc, 1, jnc_median_str);
  num_to_str(cd.max_len, 1, len_max_str);
  num_to_str(cd.max_junc, 1, jnc_max_str);
  num_to_str(cd.total_len, 1, len_total_str);
  num_to_str(cd.total_junc, 1, jnc_total_str);

  status("Lengths: mean: %s, median: %s, max: %s, total: %s [kmers]",
         len_mean_str, len_median_str, len_max_str, len_total_str);
  status("Junctions: mean: %s, median: %s, max: %s, total: %s [out >1]",
         jnc_mean_str, jnc_median_str, jnc_max_str, jnc_total_str);
  status("Max junction density: %.2f\n", cd.max_junc_density);
  // status("Contigs looping back to a kmer: %zu [%.2f%%]\n", nloop,
  //        (100.0 * nloop) / cd.ncontigs);

  timestamp(ctx_msg_out);
  message(" Outdegree: ");
  char nout_str[90];

  for(i = 0; i <= 4; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(cd.contigs_outdegree[i], nout_str),
            (size_t)((100.0*cd.contigs_outdegree[i])/(2.0*cd.ncontigs)+0.5));
  }
  message("\n");

  contig_print_path_dist(cd.paths_held, MAXPATH, "Paths held", cd.ncontigs);
  contig_print_path_dist(cd.paths_pickdup, MAXPATH, "Paths pickdup", cd.ncontigs);
  contig_print_path_dist(cd.paths_counter, MAXPATH, "Paths counter", cd.ncontigs);

  contig_data_dealloc(&cd);

  // free(visited);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free((void*)db_graph.kmer_paths);
  free(path_store);
  if(tmppdata != NULL) free(tmppdata);

  walker_dealloc(&rptwlk);
  graph_walker_dealloc(&wlk);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
