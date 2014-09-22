#include "global.h"
#include "assemble_contigs.h"
#include "db_node.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "async_read_io.h"
#include "util.h"
#include "file_util.h"

void assemble_contigs_stats_init(AssembleContigStats *stats)
{
  memset(stats, 0, sizeof(*stats));
  // Zero doubles
  stats->max_junc_density = 0;
  size_buf_capacity(&stats->lengths, 1<<20); // ~1 Million
  size_buf_capacity(&stats->junctns, 1<<20); // ~1 Million
}

void assemble_contigs_stats_destroy(AssembleContigStats *stats)
{
  size_buf_dealloc(&stats->lengths);
  size_buf_dealloc(&stats->junctns);
  memset(stats, 0, sizeof(*stats));
}

void assemble_contigs_stats_merge(AssembleContigStats *dst,
                                  const AssembleContigStats *src)
{
  ctx_assert(dst->lengths.len == dst->junctns.len);
  ctx_assert(dst->lengths.len == dst->num_contigs);
  ctx_assert(src->lengths.len == src->junctns.len);
  ctx_assert(src->lengths.len == src->num_contigs);

  size_t i;

  size_buf_append(&dst->lengths, src->lengths.data, src->lengths.len);
  size_buf_append(&dst->junctns, src->junctns.data, src->junctns.len);

  dst->num_contigs += src->num_contigs;
  dst->total_len += src->total_len;
  dst->total_junc += src->total_junc;

  for(i = 0; i < 5; i++)
    dst->contigs_outdegree[i] += src->contigs_outdegree[i];

  for(i = 0; i < AC_MAX_PATHS; i++) {
    dst->paths_held[i] += src->paths_held[i];
    dst->paths_new[i]  += src->paths_new[i];
    dst->paths_cntr[i] += src->paths_cntr[i];
  }

  dst->paths_held_max = MAX2(dst->paths_held_max, src->paths_held_max);
  dst->paths_new_max  = MAX2(dst->paths_new_max,  src->paths_new_max);
  dst->paths_cntr_max = MAX2(dst->paths_cntr_max, src->paths_cntr_max);

  for(i = 0; i < GRPHWLK_NUM_STATES; i++)
    dst->grphwlk_steps[i] += src->grphwlk_steps[i];

  dst->max_junc_density = MAX2(dst->max_junc_density, src->max_junc_density);

  dst->num_reseed_abort += src->num_reseed_abort;
  dst->num_seeds_not_found += src->num_seeds_not_found;
}

#define PREFIX "[Assembled] "

static inline void _print_grphwlk_state(const char *str, uint64_t num,
                                        size_t num_contigs)
{
  char nout_str[50];
  status(PREFIX"  %s: %s\t[ %2zu%% ]", str, ulong_to_str(num, nout_str),
         !num_contigs ? 0 : (size_t)(((100.0*num)/num_contigs)+0.5));
}

static inline void _print_path_dist(const uint64_t *hist, size_t n,
                                    const char *name, size_t num_contigs)
{
  char nout_str[100];
  size_t i;

  timestamp();
  message(PREFIX" %s: ", name);
  for(i = 0; i < n; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(hist[i], nout_str),
            (size_t)((100.0*hist[i])/(2.0*num_contigs)+0.5));
  }
  message("\n");
}


void assemble_contigs_stats_print(const AssembleContigStats *s)
{
  ctx_assert(s->lengths.len == s->junctns.len);
  ctx_assert(s->lengths.len == s->num_contigs);

  size_t i, ncontigs = s->num_contigs;

  if(ncontigs == 0) {
    status("[asm] No contigs assembled");
    return;
  }

  qsort(s->lengths.data, ncontigs, sizeof(s->lengths.data[0]), cmp_size);
  qsort(s->junctns.data, ncontigs, sizeof(s->junctns.data[0]), cmp_size);

  size_t len_n50, jnc_n50;
  size_t len_median, jnc_median, len_mean, jnc_mean;
  size_t len_min, len_max, jnc_min, jnc_max;

  // Calculate N50s
  len_n50 = calc_N50(s->lengths.data, ncontigs, s->total_len);
  jnc_n50 = calc_N50(s->junctns.data, ncontigs, s->total_junc);

  // Calculate medians, means
  len_median = MEDIAN(s->lengths.data, ncontigs);
  jnc_median = MEDIAN(s->junctns.data, ncontigs);
  len_mean = (double)s->total_len / ncontigs;
  jnc_mean = (double)s->total_junc / ncontigs;

  // Calculate min, max
  len_min = s->lengths.data[0];
  jnc_min = s->junctns.data[0];
  len_max = s->lengths.data[ncontigs-1];
  jnc_max = s->junctns.data[ncontigs-1];

  // Print number of contigs
  char num_contigs_str[50], reseed_str[50], seed_not_fnd_str[50];
  long_to_str(ncontigs, num_contigs_str);
  long_to_str(s->num_reseed_abort, reseed_str);
  long_to_str(s->num_seeds_not_found, seed_not_fnd_str);
  status(PREFIX"pulled out %s contigs", num_contigs_str);
  status(PREFIX"no-reseed aborted %s times", reseed_str);
  status(PREFIX"seed kmer not found %s times", seed_not_fnd_str);

  char len_min_str[50], len_max_str[50], len_total_str[50];
  char len_mean_str[50], len_median_str[50], len_n50_str[50];

  char jnc_min_str[50], jnc_max_str[50], jnc_total_str[50];
  char jnc_mean_str[50], jnc_median_str[50], jnc_n50_str[50];

  // Use ulong_to_str instead of num_to_str to get better accuracy
  // e.g. 966 instead of 1K
  ulong_to_str(len_mean, len_mean_str);
  ulong_to_str(jnc_mean, jnc_mean_str);
  ulong_to_str(len_median, len_median_str);
  ulong_to_str(jnc_median, jnc_median_str);
  ulong_to_str(len_n50, len_n50_str);
  ulong_to_str(jnc_n50, jnc_n50_str);
  ulong_to_str(len_min, len_min_str);
  ulong_to_str(jnc_min, jnc_min_str);
  ulong_to_str(len_max, len_max_str);
  ulong_to_str(jnc_max, jnc_max_str);
  ulong_to_str(s->total_len, len_total_str);
  ulong_to_str(s->total_junc, jnc_total_str);

  status(PREFIX"Lengths: mean: %s  median: %s  N50: %s  min: %s  max: %s  total: %s [kmers]",
         len_mean_str, len_median_str, len_n50_str, len_min_str, len_max_str, len_total_str);
  status(PREFIX"Junctions: mean: %s  median: %s  N50: %s  min: %s  max: %s  total: %s [out >1]",
         jnc_mean_str, jnc_median_str, jnc_n50_str, jnc_min_str, jnc_max_str, jnc_total_str);
  status(PREFIX"Max junction density: %.2f\n", s->max_junc_density);

  timestamp();
  message(PREFIX" Outdegree: ");
  char nout_str[50];

  for(i = 0; i <= 4; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(s->contigs_outdegree[i], nout_str),
            (size_t)((100.0*s->contigs_outdegree[i])/(2.0*ncontigs)+0.5));
  }
  message("\n");

  _print_path_dist(s->paths_held, AC_MAX_PATHS, "Paths held",    ncontigs);
  _print_path_dist(s->paths_new,  AC_MAX_PATHS, "Paths pickdup", ncontigs);
  _print_path_dist(s->paths_cntr, AC_MAX_PATHS, "Paths counter", ncontigs);

  const uint64_t *states = s->grphwlk_steps;
  size_t nsteps = s->total_len - s->num_contigs, ncontigends = 2*s->num_contigs;
  status(PREFIX"Traversal succeeded because:");
  _print_grphwlk_state("Go straight   ", states[GRPHWLK_FORWARD], nsteps);
  _print_grphwlk_state("Go colour     ", states[GRPHWLK_COLFWD], nsteps);
  _print_grphwlk_state("Go paths      ", states[GRPHWLK_USEPATH], nsteps);
  status(PREFIX"Traversal halted because:");
  _print_grphwlk_state("No coverage   ", states[GRPHWLK_NOCOVG], ncontigends);
  _print_grphwlk_state("No colour covg", states[GRPHWLK_NOCOLCOVG], ncontigends);
  _print_grphwlk_state("No paths      ", states[GRPHWLK_NOPATHS], ncontigends);
  _print_grphwlk_state("Paths split   ", states[GRPHWLK_SPLIT_PATHS], ncontigends);
  _print_grphwlk_state("Missing paths ", states[GRPHWLK_MISSING_PATHS], ncontigends);

  size_t njunc = states[GRPHWLK_USEPATH] + states[GRPHWLK_NOPATHS] +
                 states[GRPHWLK_SPLIT_PATHS] + states[GRPHWLK_MISSING_PATHS];

  status(PREFIX"Junctions:");
  _print_grphwlk_state("Paths resolved", states[GRPHWLK_USEPATH], njunc);
}

typedef struct
{
  size_t threadid, nthreads;

  GraphWalker wlk;
  RepeatWalker rptwlk;
  dBNodeBuffer nodes;
  AssembleContigStats stats;

  // Shared data
  volatile size_t *num_contig_ptr;
  size_t contig_limit;
  uint8_t *visited;
  const dBGraph *db_graph;
  size_t colour;

  // Output
  FILE *fout;
  pthread_mutex_t *outlock;
} Assembler;

// Returns 0 if keep iterating, 1 if hit assem->contig_limit
static int pulldown_contig(hkey_t hkey, Assembler *assem)
{
  AssembleContigStats *stats = &assem->stats;
  const dBGraph *db_graph = assem->db_graph;

  // Don't use a kmer if it is not in the sample we are assembling
  if(!db_node_has_col(db_graph, hkey, assem->colour)) return 0;

  // Don't use a visited kmer as a seed node if --no-reseed passed
  if(assem->visited != NULL && bitset_get_mt(assem->visited, hkey)) {
    stats->num_reseed_abort++;
    return 0;
  }

  GraphWalker *wlk = &assem->wlk;
  RepeatWalker *rptwlk = &assem->rptwlk;
  dBNodeBuffer *nodes = &assem->nodes;

  dBNode first_node = {.key = hkey, .orient = FORWARD};
  Orientation orient;
  size_t i, njunc = 0;

  db_node_buf_reset(nodes);
  db_node_buf_add(nodes, first_node);

  size_t paths_held[2] = {0}, paths_new[2] = {0}, paths_cntr[2] = {0};
  uint8_t wlk_step_last[2] = {0};

  for(orient = 0; orient < 2; orient++)
  {
    if(orient != FORWARD) {
      db_nodes_reverse_complement(nodes->data, nodes->len);
      first_node = db_node_reverse(first_node);
    }

    graph_walker_init(wlk, db_graph, assem->colour, assem->colour, first_node);

    while(graph_walker_next(wlk) && rpt_walker_attempt_traverse(rptwlk, wlk))
    {
      db_node_buf_add(nodes, wlk->node);
      stats->grphwlk_steps[wlk->last_step.status]++;
    }

    // Grab some stats
    njunc += wlk->fork_count;

    // Record numbers of paths
    paths_held[orient] = wlk->paths.len;
    paths_new[orient]  = wlk->new_paths.len;
    paths_cntr[orient] = wlk->cntr_paths.len;

    // Get failed status
    wlk_step_last[orient] = wlk->last_step.status;

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, nodes->data, nodes->len);
  }

  // If --no-seed set, mark visited nodes are visited
  // Don't use to seed another contig
  if(assem->visited != NULL) {
    for(i = 0; i < nodes->len; i++)
      (void)bitset_set_mt(assem->visited, nodes->data[i].key);
  }

  size_t contig_id;

  if(assem->fout != NULL)
  {
    // Print contig
    char kmer_str[MAX_KMER_SIZE+1], left_stat[25], rght_stat[25];
    BinaryKmer seed_bkmer = db_node_get_bkmer(db_graph, first_node.key);
    binary_kmer_to_str(seed_bkmer, db_graph->kmer_size, kmer_str);
    dna_revcomp_str(kmer_str, kmer_str, db_graph->kmer_size);

    // We have reversed the contig, so left end is now the end we hit when
    // traversing from the seed node forward... FORWARD == 0, REVERSE == 1
    graph_walker_status2str(wlk_step_last[0], left_stat, sizeof(left_stat));
    graph_walker_status2str(wlk_step_last[1], rght_stat, sizeof(rght_stat));

    // If end status is healthy, then we got stuck in a repeat
    if(grphwlk_status_is_good(wlk_step_last[0])) strcpy(left_stat, "HitRepeat");
    if(grphwlk_status_is_good(wlk_step_last[1])) strcpy(rght_stat, "HitRepeat");

    pthread_mutex_lock(assem->outlock);
    contig_id = assem->num_contig_ptr[0]++;

    if(!assem->contig_limit || contig_id < assem->contig_limit)
    {
      // Print in FASTA format with additional info in name
      fprintf(assem->fout, ">contig%zu len=%zu seed=%s "
              "lf.status=%s lf.paths.held=%zu lf.paths.new=%zu lf.paths.cntr=%zu "
              "rt.status=%s rt.paths.held=%zu rt.paths.new=%zu rt.paths.cntr=%zu\n",
              contig_id, nodes->len, kmer_str,
              left_stat, paths_held[0], paths_new[0], paths_cntr[0],
              rght_stat, paths_held[1], paths_new[1], paths_cntr[1]);

      db_nodes_print(nodes->data, nodes->len, db_graph, assem->fout);
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
  for(orient = 0; orient < 2; orient++) {
    stats->grphwlk_steps[wlk_step_last[orient]]++;
    stats->paths_held[MIN2(paths_held[orient], AC_MAX_PATHS-1)]++;
    stats->paths_new [MIN2(paths_new[orient],  AC_MAX_PATHS-1)]++;
    stats->paths_cntr[MIN2(paths_cntr[orient], AC_MAX_PATHS-1)]++;
    stats->paths_held_max = MAX2(stats->paths_held_max, paths_held[orient]);
    stats->paths_new_max  = MAX2(stats->paths_new_max,  paths_new[orient]);
    stats->paths_cntr_max = MAX2(stats->paths_cntr_max, paths_cntr[orient]);
  }

  // Out degree
  dBNode first = db_node_reverse(nodes->data[0]), last = nodes->data[nodes->len-1];

  int outdegree_fw, outdegree_rv;
  outdegree_fw = edges_get_outdegree(db_graph->col_edges[first.key], first.orient);
  outdegree_rv = edges_get_outdegree(db_graph->col_edges[last.key], last.orient);
  stats->contigs_outdegree[outdegree_fw]++;
  stats->contigs_outdegree[outdegree_rv]++;

  size_buf_add(&stats->lengths, nodes->len);
  size_buf_add(&stats->junctns, njunc);

  stats->total_len  += nodes->len;
  stats->total_junc += njunc;

  if(stats->num_contigs == 0) {
    stats->max_junc_density = (double)njunc / nodes->len;
  } else {
    stats->max_junc_density = MAX2(stats->max_junc_density,
                                   (double)njunc / nodes->len);
  }

  stats->num_contigs++;

  return 0; // 0 => keep iterating
}

static inline void _seed_rnd_kmers(void *arg)
{
  Assembler *assem = (Assembler*)arg;
  const dBGraph *db_graph = assem->db_graph;

  HASH_ITERATE_PART(&db_graph->ht, assem->threadid, assem->nthreads,
                    pulldown_contig, assem);
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
  BinaryKmer bkmer = binary_kmer_from_str(r1->seq.b, kmer_size);
  dBNode node = db_graph_find(assem->db_graph, bkmer);

  if(node.key != HASH_NOT_FOUND)
    pulldown_contig(node.key, assem);
  else
    assem->stats.num_seeds_not_found++;
}

/*!
  @param seed_files If passed, use seed kmers from sequences. If not given,
                    iterate through the hash table.
  @param contig_limit Stop after printing this many contigs, if zero no limit
  @param visited Bit array to store visited nodes in. If not NULL, do not use a
                 seed that has previously been visited. We do not clear this
                 array.
 */
void assemble_contigs(size_t nthreads,
                      seq_file_t **seed_files, size_t num_seed_files,
                      size_t contig_limit, uint8_t *visited,
                      FILE *fout, const char *out_path,
                      AssembleContigStats *stats,
                      const dBGraph *db_graph, size_t colour)
{
  ctx_assert(nthreads > 0);
  ctx_assert(!num_seed_files || seed_files);

  status("[Assemble] Assembling contigs with %zu threads, walking colour %zu",
         nthreads, colour);

  if(fout == NULL)
    status("[Assemble]   Not printing contigs");
  else
    status("[Assemble]   Writing contigs to %s", futil_outpath_str(out_path));

  Assembler *workers = ctx_calloc(nthreads, sizeof(Assembler));
  size_t i, num_contigs = 0;

  pthread_mutex_t outlock;
  if(pthread_mutex_init(&outlock, NULL) != 0) die("Mutex init failed");

  for(i = 0; i < nthreads; i++) {
    Assembler tmp = {.threadid = i, .nthreads = nthreads,
                     .num_contig_ptr = &num_contigs,
                     .contig_limit = contig_limit,
                     .db_graph = db_graph, .colour = colour,
                     .visited = visited,
                     .fout = fout, .outlock = &outlock};

    db_node_buf_alloc(&tmp.nodes, 1024);
    graph_walker_alloc(&tmp.wlk);
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
    status("[Assemble]   Sample random seed kmers");
    util_run_threads(workers, nthreads, sizeof(workers[0]),
                     nthreads, _seed_rnd_kmers);
  }

  for(i = 0; i < nthreads; i++) {
    db_node_buf_dealloc(&workers[i].nodes);
    graph_walker_dealloc(&workers[i].wlk);
    rpt_walker_dealloc(&workers[i].rptwlk);
    assemble_contigs_stats_merge(stats, &workers[i].stats);
    assemble_contigs_stats_destroy(&workers[i].stats);
  }

  pthread_mutex_destroy(&outlock);
  ctx_free(workers);
}
