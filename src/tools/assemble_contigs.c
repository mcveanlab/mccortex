#include "global.h"
#include "assemble_contigs.h"
#include "db_node.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "async_read_io.h"
#include "util.h"
#include "file_util.h"
#include "contig_confidence.h"
#include "gpath_checks.h" // gpath_fetch
#include "gpath_set.h"
#include "gpath_subset.h"

typedef struct
{
  size_t threadid, nthreads;

  GraphWalker wlk;
  RepeatWalker rptwlk;
  dBNodeBuffer nbuf;
  AssembleContigStats stats;

  GPathSet gpset;
  GPathSubset gpsubset;

  // Shared data
  volatile size_t *num_contig_ptr;
  size_t contig_limit;
  uint8_t *visited;
  bool use_missing_info_check;
  size_t *used_paths;
  const dBGraph *db_graph;
  size_t colour;
  const ContigConfidenceTable *conf_table;

  // Output
  FILE *fout;
  pthread_mutex_t *outlock;
} Assembler;

struct ContigStats {
  size_t njunc, ncycles;
  size_t wlk_steps[GRPHWLK_NUM_STATES];
  size_t paths_held[2], paths_cntr[2];
  uint8_t wlk_step_last[2];
  size_t max_step_gap[2];
  double gap_conf[2];
};

static void contig_stats_init(struct ContigStats *stats)
{
  memset(stats, 0, sizeof(struct ContigStats));
  // Initialise doubles
  stats->gap_conf[0] = stats->gap_conf[1] = 1.0;
}

// Returns 0 on success, 1 on failure
static void _assemble_contig(Assembler *assem, hkey_t hkey, const GPath *gpath,
                             struct ContigStats *results)
{
  const dBGraph *db_graph = assem->db_graph;
  GraphWalker *wlk = &assem->wlk;
  RepeatWalker *rptwlk = &assem->rptwlk;
  dBNodeBuffer *nbuf = &assem->nbuf;

  ctx_assert(assem->colour == wlk->ctxcol);
  ctx_assert(assem->colour == wlk->ctpcol);

  GraphStep step;
  size_t i;

  db_node_buf_reset(nbuf);

  if(gpath) {
    // Add gpath to nodes buffer
    dBNode first_node = {.key = hkey, .orient = gpath->orient};
    gpath_fetch(first_node, gpath, nbuf, assem->colour, db_graph);
  } else {
    dBNode first_node = {.key = hkey, .orient = FORWARD};
    db_node_buf_add(nbuf, first_node);
  }

  size_t init_len = nbuf->len;

  struct ContigStats s;
  contig_stats_init(&s);

  for(i = 0; i < 2; i++)
  {
    if(gpath) {
      printf(" Assembling %zu %zu\n", i, (size_t)gpath->orient);
    }

    ctx_assert(nbuf->len >= init_len);
    ctx_assert(!db_graph_check_all_edges(db_graph, nbuf->data, nbuf->len));

    if(i == 1)
      db_nodes_reverse_complement(nbuf->data, nbuf->len);

    // printf("prime %zu -> %zu\n", nbuf->len-init_len, nbuf->len);
    graph_walker_prime(wlk, nbuf->data+nbuf->len-init_len, init_len,
                       init_len, true);

    // graph_walker_prime(wlk, nbuf->data, init_len, init_len, i==0);
    // if(i == 1)
    //   db_nodes_reverse_complement(nbuf->data, nbuf->len);

    size_t init_junc_count = wlk->fork_count;
    bool hit_cycle = false;

    while(graph_walker_next(wlk))
    {
      db_node_buf_add(nbuf, wlk->node);

      // Do some stats
      step = wlk->last_step;
      s.wlk_steps[step.status]++;

      if(step.status == GRPHWLK_USEPATH) {
        ctx_assert(step.path_gap > 0);
        size_t read_length = step.path_gap + db_graph->kmer_size-1 + 2;
        s.max_step_gap[i] = MAX2(s.max_step_gap[i], read_length);
        s.gap_conf[i] *= conf_table_lookup(assem->conf_table, read_length);
      }

      if(!rpt_walker_attempt_traverse(rptwlk, wlk)) { hit_cycle = true; break; }
    }

    s.ncycles += hit_cycle;

    // Grab some stats
    s.njunc += wlk->fork_count - init_junc_count;

    // Record numbers of paths
    s.paths_held[i] = wlk->paths.len;
    s.paths_cntr[i] = wlk->cntr_paths.len;

    // Get failed status
    step = wlk->last_step;
    s.wlk_step_last[i] = step.status;
    if(!hit_cycle) s.wlk_steps[step.status]++;

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, nbuf->data, nbuf->len);
  }

  // If --no-reseed set, mark visited nodes are visited
  // Don't use to seed another contig
  if(gpath == NULL && assem->visited != NULL) {
    for(i = 0; i < nbuf->len; i++)
      (void)bitset_set_mt(assem->visited, nbuf->data[i].key);
  }

  memcpy(results, &s, sizeof(struct ContigStats));
}

// returns 0 on success, 1 otherwise
static int _dump_contig(Assembler *assem, hkey_t hkey, const GPath *gpath,
                        struct ContigStats s)
{
  AssembleContigStats *stats = &assem->stats;
  const dBGraph *db_graph = assem->db_graph;
  dBNodeBuffer *nbuf = &assem->nbuf;
  size_t i, contig_id;

  if(assem->fout != NULL)
  {
    // Print contig
    char kmer_str[MAX_KMER_SIZE+1], left_stat[25], rght_stat[25];
    BinaryKmer seed_bkmer = db_node_get_bkmer(db_graph, hkey);
    binary_kmer_to_str(seed_bkmer, db_graph->kmer_size, kmer_str);
    dna_revcomp_str(kmer_str, kmer_str, db_graph->kmer_size);

    // We have reversed the contig, so left end is now the end we hit when
    // traversing from the seed node forward... FORWARD == 0, REVERSE == 1
    graph_step_status2str(s.wlk_step_last[0], left_stat, sizeof(left_stat));
    graph_step_status2str(s.wlk_step_last[1], rght_stat, sizeof(rght_stat));

    // If end status is healthy, then we got stuck in a repeat
    if(grap_step_status_is_good(s.wlk_step_last[0])) strcpy(left_stat, "HitRepeat");
    if(grap_step_status_is_good(s.wlk_step_last[1])) strcpy(rght_stat, "HitRepeat");

    pthread_mutex_lock(assem->outlock);
    contig_id = assem->num_contig_ptr[0]++;

    if(!assem->contig_limit || contig_id < assem->contig_limit)
    {
      // DEV: add gpath sequence
      // Print in FASTA format with additional info in name
      fprintf(assem->fout, ">contig%zu len=%zu seed=%s %s"
              "lf.status=%s lf.paths.held=%zu lf.paths.cntr=%zu "
              "lf.max_gap=%zu lf.conf=%f "
              "rt.status=%s rt.paths.held=%zu rt.paths.cntr=%zu "
              "rf.max_gap=%zu rf.conf=%f\n",
              contig_id, nbuf->len, kmer_str, gpath ? "PATH " : "",
              left_stat, s.paths_held[0], s.paths_cntr[0], s.max_step_gap[0], s.gap_conf[0],
              rght_stat, s.paths_held[1], s.paths_cntr[1], s.max_step_gap[1], s.gap_conf[1]);

      db_nodes_print(nbuf->data, nbuf->len, db_graph, assem->fout);
      putc('\n', assem->fout);
    }

    pthread_mutex_unlock(assem->outlock);
  }
  else {
    // Lockless update
    contig_id = __sync_fetch_and_add(assem->num_contig_ptr, 1);
  }

  // Generated too many contigs - drop this one without printing or
  // recording any information/statistics on it
  if(assem->contig_limit && contig_id >= assem->contig_limit) {
    return 1; // => stop iterating
  }

  // Update statistics
  for(i = 0; i < 2; i++) {
    stats->paths_held[MIN2(s.paths_held[i], AC_MAX_PATHS-1)]++;
    stats->paths_cntr[MIN2(s.paths_cntr[i], AC_MAX_PATHS-1)]++;
    stats->paths_held_max = MAX2(stats->paths_held_max, s.paths_held[i]);
    stats->paths_cntr_max = MAX2(stats->paths_cntr_max, s.paths_cntr[i]);
  }

  for(i = 0; i < GRPHWLK_NUM_STATES; i++)
    stats->grphwlk_steps[i] += s.wlk_steps[i];

  // Out degree
  dBNode first = db_node_reverse(nbuf->data[0]), last = nbuf->data[nbuf->len-1];

  int outdegree_fw, outdegree_rv;
  outdegree_fw = edges_get_outdegree(db_graph->col_edges[first.key], first.orient);
  outdegree_rv = edges_get_outdegree(db_graph->col_edges[last.key], last.orient);
  stats->contigs_outdegree[outdegree_fw]++;
  stats->contigs_outdegree[outdegree_rv]++;

  size_buf_add(&stats->lengths, nbuf->len);
  size_buf_add(&stats->junctns, s.njunc);

  stats->num_cycles += s.ncycles;
  stats->total_len  += nbuf->len;
  stats->total_junc += s.njunc;

  if(stats->num_contigs == 0) {
    stats->max_junc_density = (double)s.njunc / nbuf->len;
  } else {
    stats->max_junc_density = MAX2(stats->max_junc_density,
                                   (double)s.njunc / nbuf->len);
  }

  stats->num_contigs++;
  return 0; // => keep iterating
}

// Returns 0 if keep iterating, 1 if hit assem->contig_limit
static int _pulldown_contig(hkey_t hkey, Assembler *assem)
{
  struct ContigStats s;

  // Don't use a kmer if it is not in the sample we are assembling
  if(!db_node_has_col(assem->db_graph, hkey, assem->colour)) return 0;

  // Don't use a visited kmer as a seed node if --no-reseed passed
  if(assem->visited != NULL && bitset_get_mt(assem->visited, hkey)) {
    assem->stats.num_reseed_abort++;
    return 0;
  }

  _assemble_contig(assem, hkey, NULL, &s);

  return _dump_contig(assem, hkey, NULL, s);
}

static inline void _seed_rnd_kmers(void *arg)
{
  Assembler *assem = (Assembler*)arg;
  const dBGraph *db_graph = assem->db_graph;

  HASH_ITERATE_PART(&db_graph->ht, assem->threadid, assem->nthreads,
                    _pulldown_contig, assem);
}

static void _seed_from_file(AsyncIOData *data, void *arg)
{
  ctx_assert2(data->r2.seq.end == 0, "Shouldn't have a second read");

  Assembler *assem = (Assembler*)arg;
  const read_t *r1 = &data->r1;
  const size_t kmer_size = assem->db_graph->kmer_size;

  if(r1->seq.end != kmer_size) {
    die("Input read length (%zu) is not kmer_size (%zu): read '%s'",
        r1->seq.end, kmer_size, r1->seq.b);
  }

  // Check all bases are valid
  const char *str;
  for(str = r1->seq.b; *str; str++)
    if(!char_is_acgt(*str))
      die("Invalid read base '%c' (%s): %s", *str, r1->name.b, r1->seq.b);

  // Find node
  dBNode node = db_graph_find_str(assem->db_graph, r1->seq.b);

  if(node.key != HASH_NOT_FOUND)
    _pulldown_contig(node.key, assem);
  else
    assem->stats.num_seeds_not_found++;
}

static int _assemble_from_paths(hkey_t hkey, Assembler *assem)
{
  const GPathStore *gpstore = &assem->db_graph->gpstore;
  const size_t ncols = gpstore->gpset.ncols, colour = assem->colour;
  GPath *gpath = gpath_store_fetch_traverse(gpstore, hkey);
  size_t i, pathid, *used_paths = assem->used_paths;
  struct ContigStats s;

  GPathSet *gpset = &assem->gpset;
  GPathSubset *gpsubset = &assem->gpsubset;

  gpath_set_reset(gpset);

  for(; gpath != NULL; gpath = gpath->next)
  {
    pathid = gpset_get_pkey(&gpstore->gpset, gpath);

    if(!gpath_has_colour(gpath, ncols, colour)) {
      printf("Path doesn't have colour\n");
    }
    else if(bitset_get(used_paths, pathid)) {
      printf("Path used\n");
    }
    else
    {
      GPathNew gpath_cpy = gpath_set_get(&gpstore->gpset, gpath);
      gpath_set_add_mt(gpset, gpath_cpy);
      printf("  Fresh Path!\n");
    }
  }

  gpath_subset_init(gpsubset, gpset);
  gpath_subset_load_set(gpsubset);
  gpath_subset_rmsubstr(gpsubset);

  GPath **list = gpsubset->list.data;

  for(i = 0; i < gpsubset->list.len; i++)
  {
    _assemble_contig(assem, hkey, list[i], &s);
    _dump_contig(assem, hkey, list[i], s);
  }

  return 0; // 0 => keep iterating
}

static void assemble_from_paths(void *arg)
{
  Assembler *assem = (Assembler*)arg;
  const dBGraph *db_graph = assem->db_graph;

  gpath_set_alloc(&assem->gpset, db_graph->gpstore.gpset.ncols,
                  ONE_MEGABYTE, true, false);

  gpath_subset_alloc(&assem->gpsubset);

  HASH_ITERATE_PART(&db_graph->ht, assem->threadid, assem->nthreads,
                    _assemble_from_paths, assem);

  gpath_set_dealloc(&assem->gpset);
  gpath_subset_dealloc(&assem->gpsubset);
}

/**
 * Assemble contig for a given sample.
 *
 * @param seed_files If passed, use seed kmers from sequences. If not given,
 *                   iterate through the hash table.
 * @param contig_limit Stop after printing this many contigs, if zero no limit
 * @param visited Bit array to store visited nodes in. If not NULL, do not use a
 *                seed that has previously been visited. We do not clear this
 *                array.
 * @param seed_with_unused_paths If set, mark paths as used once entirely
 *                               contained in a contig. Unused paths are then
 *                               used to seed contigs.
 */
void assemble_contigs(size_t nthreads,
                      seq_file_t **seed_files, size_t num_seed_files,
                      size_t contig_limit, uint8_t *visited,
                      bool use_missing_info_check, bool seed_with_unused_paths,
                      FILE *fout, const char *out_path,
                      AssembleContigStats *stats,
                      const ContigConfidenceTable *conf_table,
                      const dBGraph *db_graph, size_t colour)
{
  ctx_assert(nthreads > 0);
  ctx_assert(!num_seed_files || seed_files);
  ctx_assert(!seed_with_unused_paths || num_seed_files == 0);

  status("[Assemble] Assembling contigs with %zu threads, walking colour %zu",
         nthreads, colour);
  status("[Assemble] Using missing info check: %s",
         use_missing_info_check ? "yes" : "no");

  if(fout == NULL)
    status("[Assemble]   Not printing contigs");
  else
    status("[Assemble]   Writing contigs to %s", futil_outpath_str(out_path));

  size_t npaths = db_graph->gpstore.num_paths;
  size_t npathwords = (npaths+sizeof(size_t)*8-1)/(sizeof(size_t)*8);
  size_t *used_paths = NULL;
  if(seed_with_unused_paths) used_paths = ctx_calloc(npathwords, sizeof(size_t));

  Assembler *workers = ctx_calloc(nthreads, sizeof(Assembler));
  size_t i, num_contigs = 0;

  pthread_mutex_t outlock;
  if(pthread_mutex_init(&outlock, NULL) != 0) die("Mutex init failed");

  for(i = 0; i < nthreads; i++) {
    Assembler tmp = {.threadid = i, .nthreads = nthreads,
                     .num_contig_ptr = &num_contigs,
                     .contig_limit = contig_limit,
                     .use_missing_info_check = use_missing_info_check,
                     .used_paths = used_paths,
                     .db_graph = db_graph, .colour = colour,
                     .conf_table = conf_table,
                     .visited = visited,
                     .fout = fout, .outlock = &outlock};

    db_node_buf_alloc(&tmp.nbuf, 1024);

    graph_walker_alloc(&tmp.wlk, db_graph);
    graph_walker_setup(&tmp.wlk, use_missing_info_check, colour, colour, db_graph);
    tmp.used_paths = tmp.wlk.used_paths = used_paths;

    rpt_walker_alloc(&tmp.rptwlk, db_graph->ht.capacity, 22); // 4MB
    assemble_contigs_stats_init(&tmp.stats);

    memcpy(&workers[i], &tmp, sizeof(Assembler));
  }

  if(num_seed_files)
  {
    // Start async io reading
    status("[Assemble] Sample seed kmers from:");

    AsyncIOInput *async_tasks = ctx_calloc(num_seed_files, sizeof(AsyncIOInput));

    for(i = 0; i < num_seed_files; i++) {
      async_tasks[i].file1 = seed_files[i];
      status("[Assemble]   %s", futil_outpath_str(seed_files[i]->path));
    }

    asyncio_run_pool(async_tasks, num_seed_files, _seed_from_file,
                     workers, nthreads, sizeof(workers[0]));

    ctx_free(async_tasks);
  }
  else
  {
    // Use random kmers as seeds
    status("[Assemble] Seeding with random kmers...");
    util_run_threads(workers, nthreads, sizeof(workers[0]),
                     nthreads, _seed_rnd_kmers);

    /*
    if(seed_with_unused_paths && npaths > 0)
    {
      // Check if there are unused paths
      size_t top_bits = npaths - (npaths / sizeof(size_t)) * sizeof(size_t);
      for(i = 0; i+1 < npathwords && used_paths[i] != SIZE_MAX; i++) {}

      if(i+1 < npathwords || used_paths[npathwords-1] < bitmask64(top_bits)) {
        status("[Assemble] Seeding with unused paths...");
        util_run_threads(workers, nthreads, sizeof(workers[0]),
                         nthreads, assemble_from_paths);
      } else {
        status("[Assemble] No unused paths to seed with");
      }
    }
    */
  }

  for(i = 0; i < nthreads; i++) {
    db_node_buf_dealloc(&workers[i].nbuf);
    graph_walker_dealloc(&workers[i].wlk);
    rpt_walker_dealloc(&workers[i].rptwlk);
    assemble_contigs_stats_merge(stats, &workers[i].stats);
    assemble_contigs_stats_destroy(&workers[i].stats);
  }

  pthread_mutex_destroy(&outlock);
  ctx_free(workers);
  ctx_free(used_paths);
}
