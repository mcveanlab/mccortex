#include "global.h"
#include "assemble_contigs.h"
#include "db_node.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "async_read_io.h"
#include "util.h"
#include "file_util.h"
#include "contig_confidence.h"

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
  const ContigConfidenceTable *conf_table;

  // Output
  FILE *fout;
  pthread_mutex_t *outlock;
} Assembler;

// Returns 0 if keep iterating, 1 if hit assem->contig_limit
static int _pulldown_contig(hkey_t hkey, Assembler *assem)
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
  size_t i, njunc = 0, ncycles = 0;
  GraphStep step;

  db_node_buf_reset(nodes);
  db_node_buf_add(nodes, first_node);

  size_t wlk_steps[GRPHWLK_NUM_STATES] = {0};
  size_t paths_held[2] = {0}, paths_cntr[2] = {0};
  uint8_t wlk_step_last[2] = {0, 0};
  size_t max_step_gap[2] = {0, 0};
  double gap_conf[2] = {1, 1};

  for(orient = 0; orient < 2; orient++)
  {
    if(orient != FORWARD) {
      db_nodes_reverse_complement(nodes->data, nodes->len);
      first_node = db_node_reverse(first_node);
    }

    graph_walker_init(wlk, db_graph, assem->colour, assem->colour, first_node);

    while(graph_walker_next(wlk))
    {
      db_node_buf_add(nodes, wlk->node);

      // Do some stats
      step = wlk->last_step;
      wlk_steps[step.status]++;

      if(step.status == GRPHWLK_USEPATH) {
        ctx_assert(step.path_gap > 0);
        size_t read_length = step.path_gap + db_graph->kmer_size-1 + 2;
        max_step_gap[orient] = MAX2(max_step_gap[orient], read_length);
        gap_conf[orient] *= conf_table_lookup(assem->conf_table, read_length);
      }

      if(!rpt_walker_attempt_traverse(rptwlk, wlk)) { ncycles++; break; }
    }

    // Grab some stats
    njunc += wlk->fork_count;

    // Record numbers of paths
    paths_held[orient] = wlk->paths.len;
    paths_cntr[orient] = wlk->cntr_paths.len;

    // Get failed status
    step = wlk->last_step;
    wlk_step_last[orient] = step.status;
    wlk_steps[step.status]++;

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
    graph_step_status2str(wlk_step_last[0], left_stat, sizeof(left_stat));
    graph_step_status2str(wlk_step_last[1], rght_stat, sizeof(rght_stat));

    // If end status is healthy, then we got stuck in a repeat
    if(grap_step_status_is_good(wlk_step_last[0])) strcpy(left_stat, "HitRepeat");
    if(grap_step_status_is_good(wlk_step_last[1])) strcpy(rght_stat, "HitRepeat");

    pthread_mutex_lock(assem->outlock);
    contig_id = assem->num_contig_ptr[0]++;

    if(!assem->contig_limit || contig_id < assem->contig_limit)
    {
      // Print in FASTA format with additional info in name
      fprintf(assem->fout, ">contig%zu len=%zu seed=%s "
              "lf.status=%s lf.paths.held=%zu lf.paths.cntr=%zu "
              "lf.max_gap=%zu lf.conf=%f "
              "rt.status=%s rt.paths.held=%zu rt.paths.cntr=%zu "
              "rf.max_gap=%zu rf.conf=%f\n",
              contig_id, nodes->len, kmer_str,
              left_stat, paths_held[0], paths_cntr[0], max_step_gap[0], gap_conf[0],
              rght_stat, paths_held[1], paths_cntr[1], max_step_gap[1], gap_conf[1]);

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
    stats->paths_held[MIN2(paths_held[orient], AC_MAX_PATHS-1)]++;
    stats->paths_cntr[MIN2(paths_cntr[orient], AC_MAX_PATHS-1)]++;
    stats->paths_held_max = MAX2(stats->paths_held_max, paths_held[orient]);
    stats->paths_cntr_max = MAX2(stats->paths_cntr_max, paths_cntr[orient]);
  }

  for(i = 0; i < GRPHWLK_NUM_STATES; i++)
    stats->grphwlk_steps[i] += wlk_steps[i];

  // Out degree
  dBNode first = db_node_reverse(nodes->data[0]), last = nodes->data[nodes->len-1];

  int outdegree_fw, outdegree_rv;
  outdegree_fw = edges_get_outdegree(db_graph->col_edges[first.key], first.orient);
  outdegree_rv = edges_get_outdegree(db_graph->col_edges[last.key], last.orient);
  stats->contigs_outdegree[outdegree_fw]++;
  stats->contigs_outdegree[outdegree_rv]++;

  size_buf_add(&stats->lengths, nodes->len);
  size_buf_add(&stats->junctns, njunc);

  stats->num_cycles += ncycles;
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
  BinaryKmer bkmer = binary_kmer_from_str(r1->seq.b, kmer_size);
  dBNode node = db_graph_find(assem->db_graph, bkmer);

  if(node.key != HASH_NOT_FOUND)
    _pulldown_contig(node.key, assem);
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
                      const ContigConfidenceTable *conf_table,
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
                     .conf_table = conf_table,
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
