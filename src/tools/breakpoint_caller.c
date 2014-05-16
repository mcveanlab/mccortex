#include "global.h"
#include "breakpoint_caller.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "kmer_occur.h"
#include "graph_crawler.h"

typedef struct {
  size_t first_runid, num_runs;
} PathRefRun;

typedef struct
{
  // Specific to this instance
  const size_t threadid, nthreads;

  // Temporary memory used by this instance
  GraphCrawler crawlers[2]; // [0] => FORWARD, [1] => REVERSE

  // Paths are built one at a time for each colour,
  // going in each direction (forward for allele+3pflank, reverse for 5pflank)
  // Once all paths are constructed, GraphCrawler collapses them down to merge
  // identical ones.
  // We associate runs along the refence with each path before it is collapsed
  // down, which is stored in allele_refs, allele_run_buf for alleles.
  // If we can collapse down earlier, following Kmer occurance runs along
  // paths could be sped up. This would be fiddly to do though.

  // Flanks and paths
  KOccurRunBuffer koruns_5p, koruns_5p_ended;
  KOccurRunBuffer koruns_3p, koruns_3p_ended;
  dBNodeBuffer flank5pbuf, allelebuf;

  // Where the paths meet the ref
  PathRefRun *allele_refs, *flank5p_refs;
  KOccurRunBuffer allele_run_buf, flank5p_run_buf;

  // Passed to all instances
  const KOGraph kograph;
  const dBGraph *db_graph;
  gzFile gzout;
  pthread_mutex_t *const out_lock;
  size_t *callid;
  const size_t min_ref_nkmers, max_ref_nkmers; // how many kmers of homology req
} BreakpointCaller;

static BreakpointCaller* brkpt_callers_new(size_t num_callers, gzFile gzout,
                                           size_t min_ref_flank,
                                           size_t max_ref_flank,
                                           const KOGraph kograph,
                                           const dBGraph *db_graph)
{
  const size_t ncols = db_graph->num_of_cols;
  BreakpointCaller *callers = ctx_malloc(num_callers * sizeof(BreakpointCaller));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  size_t *callid = ctx_calloc(1, sizeof(size_t));

  // Each colour in each caller can have a GraphCache path at once
  PathRefRun *path_ref_runs = ctx_calloc(num_callers*ncols*2, sizeof(PathRefRun));

  size_t i;
  for(i = 0; i < num_callers; i++)
  {
    BreakpointCaller tmp = {.threadid = i,
                            .nthreads = num_callers,
                            .kograph = kograph,
                            .db_graph = db_graph,
                            .gzout = gzout,
                            .out_lock = out_lock,
                            .callid = callid,
                            .allele_refs = path_ref_runs+i*ncols*2,
                            .flank5p_refs = path_ref_runs+i*ncols*2+ncols,
                            .min_ref_nkmers = min_ref_flank,
                            .max_ref_nkmers = max_ref_flank};

    memcpy(&callers[i], &tmp, sizeof(BreakpointCaller));

    db_node_buf_alloc(&callers[i].allelebuf, 1024);
    db_node_buf_alloc(&callers[i].flank5pbuf, 1024);
    kmer_run_buf_alloc(&callers[i].koruns_5p, 128);
    kmer_run_buf_alloc(&callers[i].koruns_5p_ended, 128);
    kmer_run_buf_alloc(&callers[i].koruns_3p, 128);
    kmer_run_buf_alloc(&callers[i].koruns_3p_ended, 128);
    kmer_run_buf_alloc(&callers[i].allele_run_buf, 128);
    kmer_run_buf_alloc(&callers[i].flank5p_run_buf, 128);
    graph_crawler_alloc(&callers[i].crawlers[0], db_graph);
    graph_crawler_alloc(&callers[i].crawlers[1], db_graph);
  }

  return callers;
}

static void brkpt_callers_destroy(BreakpointCaller *callers, size_t num_callers)
{
  size_t i;
  for(i = 0; i < num_callers; i++) {
    db_node_buf_dealloc(&callers[i].allelebuf);
    db_node_buf_dealloc(&callers[i].flank5pbuf);
    kmer_run_buf_dealloc(&callers[i].koruns_5p);
    kmer_run_buf_dealloc(&callers[i].koruns_5p_ended);
    kmer_run_buf_dealloc(&callers[i].koruns_3p);
    kmer_run_buf_dealloc(&callers[i].koruns_3p_ended);
    kmer_run_buf_dealloc(&callers[i].allele_run_buf);
    kmer_run_buf_dealloc(&callers[i].flank5p_run_buf);
    graph_crawler_dealloc(&callers[i].crawlers[0]);
    graph_crawler_dealloc(&callers[i].crawlers[1]);
  }
  pthread_mutex_destroy(callers[0].out_lock);
  ctx_free(callers[0].out_lock);
  ctx_free(callers[0].callid);
  ctx_free(callers[0].allele_refs);
  ctx_free(callers);
}

static inline void ko_gzprint(gzFile gzout, size_t kmer_size,
                              KOGraph kograph, KOccurRun korun,
                              size_t contig_start)
{
  const char strand[] = {'+','-'};
  const char *chrom = kograph_chrom(kograph,korun).name;
  ctx_assert(korun.first <= korun.last || korun.strand == STRAND_MINUS);
  // get end coord as inclusive coord as start-end
  // (Note: start may be greater than end if strand is minus)
  size_t start, end, qoffset;
  if(korun.strand == STRAND_PLUS) {
    start = korun.first;
    end = korun.last+kmer_size-1;
  } else {
    start = korun.first+kmer_size-1;
    end = korun.last;
  }
  qoffset = korun.qoffset - contig_start;
  // +1 to coords to convert to 1-based
  gzprintf(gzout, "%s:%zu-%zu:%c:%u",
           chrom, start+1, end+1, strand[korun.strand], qoffset+1);
}

static inline void koruns_gzprint(gzFile gzout, size_t kmer_size,
                                  KOGraph kograph,
                                  const KOccurRun *koruns, size_t n,
                                  size_t contig_start)
{
  size_t i;
  if(n == 0) return;
  ko_gzprint(gzout, kmer_size, kograph, koruns[0], contig_start);
  for(i = 1; i < n; i++) {
    gzputc(gzout, ',');
    ko_gzprint(gzout, kmer_size, kograph, koruns[i], contig_start);
  }
}

static void process_contig(BreakpointCaller *caller,
                           const uint32_t *cols, size_t ncols,
                           const dBNodeBuffer *flank5p,
                           const dBNodeBuffer *allelebuf,
                           const KOccurRun *flank5p_runs, size_t num_flank5p_runs,
                           const KOccurRun *flank3p_runs, size_t num_flank3p_runs)
{
  gzFile gzout = caller->gzout;
  KOGraph kograph = caller->kograph;
  const size_t kmer_size = caller->db_graph->kmer_size;

  ctx_assert(ncols > 0);

  // we never re-met the ref
  if(num_flank3p_runs == 0) return;

  // Find first place we meet the ref
  size_t i, end = flank3p_runs[0].qoffset;
  size_t callid = __sync_fetch_and_add((volatile size_t*)caller->callid, 1);

  // Swallow up some of the path into the 3p flank
  size_t shift3p = MIN2(kmer_size-1, end);
  end -= shift3p;

  pthread_mutex_lock(caller->out_lock);

  // 5p flank with list of ref intersections
  gzprintf(gzout, ">call.%zu.5pflank chr=", callid);
  koruns_gzprint(gzout, kmer_size, kograph, flank5p_runs, num_flank5p_runs, 0);
  gzputc(gzout, '\n');
  db_nodes_gzprint(flank5p->data, flank5p->len, caller->db_graph, gzout);
  gzputc(gzout, '\n');

  // 3p flank with list of ref intersections
  gzprintf(gzout, ">call.%zu.3pflank chr=", callid);
  koruns_gzprint(gzout, kmer_size, kograph, flank3p_runs, num_flank3p_runs, end+shift3p);
  gzputc(gzout, '\n');
  db_nodes_gzprint_cont(allelebuf->data+end, allelebuf->len-end, caller->db_graph, gzout);
  gzputc(gzout, '\n');

  // Print path with list of colours
  gzprintf(gzout, ">call.%zu.path cols=%zu", callid, cols[0]);
  for(i = 1; i < ncols; i++) gzprintf(gzout, ",%zu", cols[i]);
  gzputc(gzout, '\n');
  db_nodes_gzprint_cont(allelebuf->data, end, caller->db_graph, gzout);
  gzprintf(gzout, "\n\n");

  pthread_mutex_unlock(caller->out_lock);
}

// If `pickup_new_runs` is true we pick up runs starting at this supernode
static inline bool gcrawler_stop_at_ref_covg(const GraphCache *cache,
                                             const GCacheStep *step,
                                             BreakpointCaller *caller,
                                             KOccurRunBuffer *koruns,
                                             KOccurRunBuffer *koruns_ended,
                                             bool pickup_new_runs)
{
  GCacheSnode *snode = graph_cache_snode(cache, step->supernode);
  GCachePath *path = graph_cache_path(cache, step->pathid);
  const dBNode *nodes = graph_cache_first_node(cache, snode);
  bool forward = (step->orient == FORWARD);

  // Use index of last step as qoffset
  size_t qoffset = path->num_steps-1;

  // Kmer occurance runs are added to koruns_3p_ended only if they end and are
  // longer than the mininum length in kmers (caller->min_ref_nkmers)
  kograph_filter_extend(caller->kograph,
                        nodes, snode->num_nodes, forward,
                        caller->min_ref_nkmers, qoffset,
                        koruns, koruns_ended,
                        pickup_new_runs);

  size_t i, max_ref_run = 0;
  size_t min_run_qoffset = SIZE_MAX, min_ended_run_qoffset = SIZE_MAX;
  for(i = 0; i < koruns->len; i++) {
    max_ref_run = MAX2(max_ref_run, korun_len(koruns->data[i]));
    min_run_qoffset = MIN2(min_run_qoffset, koruns->data[i].qoffset);
  }

  // Stop if all our earliest runs have finished
  for(i = 0; i < koruns_ended->len; i++) {
    min_ended_run_qoffset = MIN2(min_ended_run_qoffset, koruns_ended->data[i].qoffset);
  }

  // Continue if...
  return (koruns_ended->len == 0 && max_ref_run < caller->min_ref_nkmers &&
          min_run_qoffset <= min_ended_run_qoffset);
}

// Try to pick up new runs at each supernode
static bool gcrawler_path_stop_at_ref_covg(GraphCache *cache, GCacheStep *step,
                                           void *arg)
{
  BreakpointCaller *caller = (BreakpointCaller*)arg;

  return gcrawler_stop_at_ref_covg(cache, step, caller,
                                   &caller->koruns_3p,
                                   &caller->koruns_3p_ended,
                                   true);
}

// For 5p flank only pick up new runs starting at the first supernode
static bool gcrawler_flank5p_stop_at_ref_covg(GraphCache *cache, GCacheStep *step,
                                              void *arg)
{
  BreakpointCaller *caller = (BreakpointCaller*)arg;
  const GCachePath *path = graph_cache_path(cache, step->pathid);
  const GCacheStep *first_step = graph_cache_step(cache, path->first_step);
  bool pickup_new_runs = (step == first_step);

  return gcrawler_stop_at_ref_covg(cache, step, caller,
                                   &caller->koruns_5p,
                                   &caller->koruns_5p_ended,
                                   pickup_new_runs) &&
         (caller->koruns_5p.len > 0);
}

static inline void gcrawler_finish_ref_covg(BreakpointCaller *caller,
                                            uint32_t pathid,
                                            KOccurRunBuffer *koruns,
                                            KOccurRunBuffer *koruns_ended,
                                            KOccurRunBuffer *runs_buf,
                                            PathRefRun *ref_runs)
{
  size_t init_len = runs_buf->len;

  // Copy finished runs into array
  kmer_run_buf_ensure_capacity(runs_buf, runs_buf->len+koruns->len+koruns_ended->len);
  kmer_run_buf_append(runs_buf, koruns_ended->data, koruns_ended->len);

  runs_buf->len += koruns_filter(koruns->data, koruns->len,
                                 runs_buf->data+runs_buf->len,
                                 caller->min_ref_nkmers);

  kmer_run_buf_reset(koruns);
  kmer_run_buf_reset(koruns_ended);

  ref_runs[pathid].first_runid = init_len;
  ref_runs[pathid].num_runs = runs_buf->len - init_len;
}

static void gcrawler_path_finish_ref_covg(GraphCache *cache, uint32_t pathid,
                                          void *arg)
{
  (void)cache; // this function passed for callback, don't need cache here
  BreakpointCaller *caller = (BreakpointCaller*)arg;

  gcrawler_finish_ref_covg(caller, pathid,
                           &caller->koruns_3p,
                           &caller->koruns_3p_ended,
                           &caller->allele_run_buf,
                           caller->allele_refs);
}

static void gcrawler_flank5p_finish_ref_covg(GraphCache *cache, uint32_t pathid,
                                             void *arg)
{
  (void)cache; // this function passed for callback, don't need cache here
  BreakpointCaller *caller = (BreakpointCaller*)arg;

  gcrawler_finish_ref_covg(caller, pathid,
                           &caller->koruns_5p,
                           &caller->koruns_5p_ended,
                           &caller->flank5p_run_buf,
                           caller->flank5p_refs);
}

// Traverse from node0 -> node1
static void traverse_5pflank(BreakpointCaller *caller, GraphCrawler *crawler,
                             dBNode node0, dBNode node1)
{
  const dBGraph *db_graph = crawler->cache.db_graph;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  size_t i, num_next;
  BinaryKmer bkmer0 = db_node_get_bkmer(db_graph, node0.key);

  num_next = db_graph_next_nodes(db_graph, bkmer0, node0.orient,
                                 db_node_edges(db_graph, node0.key, 0),
                                 next_nodes, next_nucs);

  // Find index of previous node
  for(i = 0; i < num_next && !db_nodes_are_equal(next_nodes[i],node1); i++) {}

  ctx_assert(i < num_next && db_nodes_are_equal(next_nodes[i],node1));

  kmer_run_buf_reset(&caller->koruns_5p);
  kmer_run_buf_reset(&caller->koruns_5p_ended);
  kmer_run_buf_reset(&caller->flank5p_run_buf);

  // Go backwards to get 5p flank
  // NULL means loop from 0..(ncols-1)
  graph_crawler_fetch(crawler, node0,
                      next_nodes, next_nucs, i, num_next,
                      NULL, db_graph->num_of_cols,
                      gcrawler_flank5p_stop_at_ref_covg,
                      gcrawler_flank5p_finish_ref_covg,
                      caller);
}

static inline
KOccurRun* fetch_ref_contact(const GraphCache *cache, uint32_t pathid,
                             const PathRefRun *ref_runs,
                             KOccurRunBuffer *runbuf)
{
  // Get path
  const GCachePath *path = graph_cache_path(cache, pathid);
  const GCacheStep *steps = graph_cache_step(cache, path->first_step);
  size_t num_steps = path->num_steps;

  // Get runs along the ref
  PathRefRun ref_run = ref_runs[pathid];
  size_t num_runs = ref_run.num_runs;
  KOccurRun *koruns = runbuf->data + ref_run.first_runid;
  koruns_sort_by_qoffset(koruns, num_runs);

  // Set qoffset to be the kmer offset in the path
  size_t r, s, offset = 0;

  for(s = r = 0; s < num_steps; s++) {
    for(; r < num_runs && koruns[r].qoffset == s; r++) {
      koruns[r].qoffset = offset;
    }

    if(r == num_runs) break;

    const GCacheSnode *snode = graph_cache_snode(cache, steps[s].supernode);
    offset += snode->num_nodes;
  }

  return koruns;
}

// Walk the graph remembering the last time we met the ref
// When traversal fails, dump sequence up to last meeting with the ref
static void follow_break(BreakpointCaller *caller, dBNode node)
{
  size_t i, j, k, num_next;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  const dBGraph *db_graph = caller->db_graph;

  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);
  Edges edges = db_node_get_edges(db_graph, node.key, 0);

  num_next = db_graph_next_nodes(db_graph, bkey, node.orient, edges,
                                 next_nodes, next_nucs);

  // Filter out next nodes in the reference
  for(i = 0, j = 0; i < num_next; i++) {
    if(kograph_num(caller->kograph, next_nodes[i].key) == 0) {
      next_nodes[j] = next_nodes[i];
      next_nucs[j] = next_nucs[i];
      j++;
    }
  }

  // Abandon if all options are in ref or none are
  if(j == num_next || j == 0) return;
  num_next = j;

  // Follow all paths not in ref, in all colours
  GraphCrawler *fw_crawler = &caller->crawlers[node.orient];
  GraphCrawler *rv_crawler = &caller->crawlers[!node.orient];
  dBNodeBuffer *allelebuf = &caller->allelebuf, *flank5pbuf = &caller->flank5pbuf;
  GCMultiColPath *flank5p_multicolpath, *allele_multicolpath;
  KOccurRun *flank5p_runs, *flank3p_runs;
  size_t flank5p_pathid, allele_pathid;
  size_t num_flank5p_runs, num_flank3p_runs;

  // We fetch 5' flanks in all colours then merge matching paths
  // we stop fetching a single path if it stops tracking the reference
  // Alternatively, we could fetch the 5' flank in everyone and stop after a
  // given distance, then check for that set of paths how much it tracks the
  // reference. This has the advantage of scaling much better with number of
  // samples, but not so well as min_ref_nkmers increases (since we fetch
  // many flanks that can't be used) - I think this is less of a worry.

  // Loop over possible next nodes at this junction
  for(i = 0; i < num_next; i++)
  {
    // Go backwards to get 5p flank
    traverse_5pflank(caller, rv_crawler, db_node_reverse(next_nodes[i]),
                     db_node_reverse(node));

    // Loop over the flanks we got
    for(j = 0; j < rv_crawler->num_paths; j++)
    {
      // Get 5p flank
      db_node_buf_reset(flank5pbuf);
      graph_crawler_get_path_nodes(rv_crawler, j, flank5pbuf);
      flank5p_multicolpath = &rv_crawler->multicol_paths[j];
      flank5p_pathid = flank5p_multicolpath->pathid;

      // Fetch 3pflank ref position
      num_flank5p_runs = caller->flank5p_refs[flank5p_pathid].num_runs;
      flank5p_runs = fetch_ref_contact(&rv_crawler->cache, flank5p_pathid,
                                       caller->flank5p_refs,
                                       &caller->flank5p_run_buf);

      koruns_reverse(flank5p_runs, num_flank5p_runs);
      koruns_sort_by_qoffset(flank5p_runs, num_flank5p_runs);
      db_nodes_reverse_complement(flank5pbuf->data, flank5pbuf->len);

      if(num_flank5p_runs > 0)
      {
        // Reset caller
        kmer_run_buf_reset(&caller->koruns_3p);
        kmer_run_buf_reset(&caller->koruns_3p_ended);
        kmer_run_buf_reset(&caller->allele_run_buf);

        // functions gcrawler_path_stop_at_ref_covg(),
        //           gcrawler_path_finish_ref_covg()
        // both fill koruns_3p, koruns_3p_ended and allele_run_buf

        // Only traverse in the colours we have a flank for
        graph_crawler_fetch(fw_crawler, node,
                            next_nodes, next_nucs, j, num_next,
                            flank5p_multicolpath->cols,
                            flank5p_multicolpath->num_cols,
                            gcrawler_path_stop_at_ref_covg,
                            gcrawler_path_finish_ref_covg,
                            caller);

        // Assemble contigs - fetch forwards for each path for given 5p flank
        for(k = 0; k < fw_crawler->num_paths; k++)
        {
          // Fetch nodes
          db_node_buf_reset(allelebuf);
          graph_crawler_get_path_nodes(fw_crawler, k, allelebuf);

          // DEV: this should really be troo
          // ctx_assert(allelebuf->len > 0);
          if(allelebuf->len > 0)
          {
            allele_multicolpath = &fw_crawler->multicol_paths[k];
            allele_pathid = allele_multicolpath->pathid;

            // Fetch 3pflank ref position
            num_flank3p_runs = caller->allele_refs[allele_pathid].num_runs;
            flank3p_runs = fetch_ref_contact(&fw_crawler->cache, allele_pathid,
                                             caller->allele_refs,
                                             &caller->allele_run_buf);

            process_contig(caller,
                           allele_multicolpath->cols,
                           allele_multicolpath->num_cols,
                           flank5pbuf, allelebuf,
                           flank5p_runs, num_flank5p_runs,
                           flank3p_runs, num_flank3p_runs);
          }
        }
      }
    }
  }
}

void breakpoint_caller_node(hkey_t hkey, BreakpointCaller *caller)
{
  Edges edges;

  graph_crawler_reset(&caller->crawlers[0]);
  graph_crawler_reset(&caller->crawlers[1]);

  // check node is in the ref
  if(kograph_num(caller->kograph, hkey) > 0) {
    edges = db_node_get_edges(caller->db_graph, hkey, 0);
    if(edges_get_outdegree(edges, FORWARD) > 1) {
      follow_break(caller, (dBNode){.key = hkey, .orient = FORWARD});
    }
    if(edges_get_outdegree(edges, REVERSE) > 1) {
      follow_break(caller, (dBNode){.key = hkey, .orient = REVERSE});
    }
  }
}

static void breakpoint_caller(void *ptr)
{
  BreakpointCaller *caller = (BreakpointCaller*)ptr;
  ctx_assert(caller->db_graph->num_edge_cols == 1);

  HASH_ITERATE_PART(&caller->db_graph->ht, caller->threadid, caller->nthreads,
                    breakpoint_caller_node, caller);
}

static void breakpoints_print_header(gzFile gzout, const CmdArgs *args,
                                     const char **seq_paths, size_t nseq_paths)
{
  char datestr[9], cwd[PATH_MAX + 1];
  size_t i;

  ctx_assert(nseq_paths > 0);

  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  gzprintf(gzout, "##fileFormat=CtxBreakpointsv0.1\n##fileDate=%s\n", datestr);
  gzprintf(gzout, "##cmd=\"%s\"\n", args->cmdline);
  if(futil_get_current_dir(cwd) != NULL) gzprintf(gzout, "##wkdir=%s\n", cwd);

  gzprintf(gzout, "##reference=%s", seq_paths[0]);
  for(i = 1; i < nseq_paths; i++) gzprintf(gzout, ":%s", seq_paths[i]);
  gzputc(gzout, '\n');

  gzprintf(gzout, "##ctxVersion=\""VERSION_STATUS_STR"\"\n");
  gzprintf(gzout, "##ctxKmerSize=%i\n", MAX_KMER_SIZE);
}

void breakpoints_call(size_t num_of_threads,
                      const read_t *reads, size_t num_reads,
                      gzFile gzout, const char *out_path,
                      const char **seq_paths, size_t num_seq_paths,
                      size_t min_ref_flank, size_t max_ref_flank,
                      const CmdArgs *args, dBGraph *db_graph)
{
  KOGraph kograph = kograph_create(reads, num_reads, true,
                                   num_of_threads, db_graph);

  BreakpointCaller *callers = brkpt_callers_new(num_of_threads, gzout,
                                                min_ref_flank, max_ref_flank,
                                                kograph, db_graph);

  status("Running BreakpointCaller with %zu thread%s, output to: %s",
         num_of_threads, util_plural_str(num_of_threads),
         futil_outpath_str(out_path));

  status("  Finding breakpoints after at least %zu kmers (%zubp) of homology",
         min_ref_flank, min_ref_flank+db_graph->kmer_size-1);

  breakpoints_print_header(gzout, args, seq_paths, num_seq_paths);

  util_run_threads(callers, num_of_threads, sizeof(callers[0]),
                   num_of_threads, breakpoint_caller);

  char call_num_str[100];
  ulong_to_str(callers[0].callid[0], call_num_str);
  status("  %s calls printed to %s", call_num_str, futil_outpath_str(out_path));

  brkpt_callers_destroy(callers, num_of_threads);
  kograph_free(kograph);
}
