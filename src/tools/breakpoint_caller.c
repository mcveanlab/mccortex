#include "global.h"
#include "breakpoint_caller.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "kmer_occur.h"
#include "graph_crawler.h"
#include "json_hdr.h"

typedef struct {
  uint32_t first_runid, num_runs;
} PathRefRun;

typedef struct
{
  // Specific to this instance
  const size_t nthreads;

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
  KOccurRunBuffer koruns_tmp;
  dBNodeBuffer flank5pbuf, allelebuf;

  // Where the paths meet the ref
  PathRefRun *allele_refs, *flank5p_refs;
  KOccurRunBuffer allele_run_buf, flank5p_run_buf;

  // Passed to all instances
  const KOGraph *kograph;
  const dBGraph *db_graph;
  gzFile gzout;
  pthread_mutex_t *const out_lock;
  size_t *callid;
  const size_t min_ref_nkmers, max_ref_nkmers; // how many kmers of homology req
} BreakpointCaller;

// We clear the graph cache after each fork is dealt with, so there should
// be an upper bound on how many paths the graph crawler generates.
// This upper bound is used when allocating the memory to store ref runs
// (stretches where the ref runs along with the sample)
#define MAX_REFRUNS_PER_ORIENT(ncols) ((ncols)*4)
#define MAX_REFRUNS_PER_CALLER(ncols) MAX_REFRUNS_PER_ORIENT(ncols)*2

static BreakpointCaller* brkpt_callers_new(size_t num_callers,
                                           gzFile gzout,
                                           size_t min_ref_nkmers,
                                           size_t max_ref_nkmers,
                                           const KOGraph *kograph,
                                           const dBGraph *db_graph)
{
  ctx_assert(num_callers > 0);

  const size_t ncols = db_graph->num_of_cols;
  BreakpointCaller *callers = ctx_malloc(num_callers * sizeof(BreakpointCaller));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  size_t *callid = ctx_calloc(1, sizeof(size_t));

  // Each colour in each caller can have a GraphCache path at once
  PathRefRun *path_ref_runs = ctx_calloc(num_callers*MAX_REFRUNS_PER_CALLER(ncols),
                                         sizeof(PathRefRun));

  size_t i;
  for(i = 0; i < num_callers; i++)
  {
    BreakpointCaller tmp = {.nthreads = num_callers,
                            .kograph = kograph,
                            .db_graph = db_graph,
                            .gzout = gzout,
                            .out_lock = out_lock,
                            .callid = callid,
                            .allele_refs = path_ref_runs,
                            .flank5p_refs = path_ref_runs+MAX_REFRUNS_PER_ORIENT(ncols),
                            .min_ref_nkmers = min_ref_nkmers,
                            .max_ref_nkmers = max_ref_nkmers};

    memcpy(&callers[i], &tmp, sizeof(BreakpointCaller));

    path_ref_runs += MAX_REFRUNS_PER_CALLER(ncols);

    db_node_buf_alloc(&callers[i].allelebuf, 1024);
    db_node_buf_alloc(&callers[i].flank5pbuf, 1024);
    korun_buf_alloc(&callers[i].koruns_5p, 128);
    korun_buf_alloc(&callers[i].koruns_5p_ended, 128);
    korun_buf_alloc(&callers[i].koruns_3p, 128);
    korun_buf_alloc(&callers[i].koruns_3p_ended, 128);
    korun_buf_alloc(&callers[i].koruns_tmp, 128);
    korun_buf_alloc(&callers[i].allele_run_buf, 128);
    korun_buf_alloc(&callers[i].flank5p_run_buf, 128);
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
    korun_buf_dealloc(&callers[i].koruns_5p);
    korun_buf_dealloc(&callers[i].koruns_5p_ended);
    korun_buf_dealloc(&callers[i].koruns_3p);
    korun_buf_dealloc(&callers[i].koruns_3p_ended);
    korun_buf_dealloc(&callers[i].koruns_tmp);
    korun_buf_dealloc(&callers[i].allele_run_buf);
    korun_buf_dealloc(&callers[i].flank5p_run_buf);
    graph_crawler_dealloc(&callers[i].crawlers[0]);
    graph_crawler_dealloc(&callers[i].crawlers[1]);
  }
  pthread_mutex_destroy(callers[0].out_lock);
  ctx_free(callers[0].out_lock);
  ctx_free(callers[0].callid);
  ctx_free(callers[0].allele_refs);
  ctx_free(callers);
}

static void process_contig(BreakpointCaller *caller,
                           const uint32_t *cols, size_t ncols,
                           const dBNodeBuffer *flank5p,
                           const dBNodeBuffer *allelebuf,
                           const KOccurRun *flank5p_runs, size_t nflank5p_runs,
                           const KOccurRun *flank3p_runs, size_t nflank3p_runs)
{
  gzFile gzout = caller->gzout;
  const KOGraph *kograph = caller->kograph;
  const size_t kmer_size = caller->db_graph->kmer_size;

  ctx_assert(ncols > 0);

  // we never re-met the ref
  if(nflank3p_runs == 0) { /*status("  No 3p");*/ return; }

  // status("  got a call");

  // Find first place we meet the ref
  size_t callid = __sync_fetch_and_add((volatile size_t*)caller->callid, 1);

  // Swallow up some of the path into the 3p flank
  size_t i, flank3pidx = flank3p_runs[0].qoffset;
  size_t extra3pbases = MIN2(kmer_size-1, flank3pidx);
  size_t num_path_kmers = flank3pidx - extra3pbases;
  size_t kmer3poffset = kmer_size-1-extra3pbases;

  pthread_mutex_lock(caller->out_lock);

  // This can be set to anything without a '.' in it
  const char prefix[] = "call";

  // 5p flank with list of ref intersections
  gzprintf(gzout, ">brkpnt.%s%zu.5pflank chr=", prefix, callid);
  koruns_gzprint(gzout, kmer_size, kograph, flank5p_runs, nflank5p_runs, 0, 0);
  gzputc(gzout, '\n');
  db_nodes_gzprint(flank5p->b, flank5p->len, caller->db_graph, gzout);
  gzputc(gzout, '\n');

  // 3p flank with list of ref intersections
  gzprintf(gzout, ">brkpnt.%s%zu.3pflank chr=", prefix, callid);
  koruns_gzprint(gzout, kmer_size, kograph, flank3p_runs, nflank3p_runs,
                 flank3pidx, kmer3poffset);
  gzputc(gzout, '\n');
  db_nodes_gzprint_cont(allelebuf->b+num_path_kmers,
                        allelebuf->len-num_path_kmers,
                        caller->db_graph, gzout);
  gzputc(gzout, '\n');

  // Print path with list of colours
  gzprintf(gzout, ">brkpnt.%s%zu.path cols=%zu", prefix, callid, cols[0]);
  for(i = 1; i < ncols; i++) gzprintf(gzout, ",%zu", cols[i]);
  gzputc(gzout, '\n');
  db_nodes_gzprint_cont(allelebuf->b, num_path_kmers, caller->db_graph, gzout);
  gzprintf(gzout, "\n\n");

  pthread_mutex_unlock(caller->out_lock);
}


static inline bool gcrawler_stop_at_ref_covg(const GraphCache *cache,
                                             const GCacheStep *step,
                                             BreakpointCaller *caller,
                                             KOccurRunBuffer *koruns,
                                             KOccurRunBuffer *koruns_tmp,
                                             KOccurRunBuffer *koruns_ended)
{
  const GCacheUnitig *unitig = gc_step_get_unitig(cache, step);
  const GCachePath *path = gc_step_get_path(cache, step);
  const dBNode *nodes = gc_unitig_get_nodes(cache, unitig);
  bool forward = (step->orient == FORWARD);

  // Get node index of first node in the last step
  size_t qoffset = gc_path_get_nkmers(cache, path) -
                   gc_step_get_nkmers(cache, gc_path_last_step(cache, path));

  // Kmer occurance runs are added to koruns_3p_ended only if they end and are
  // longer than the mininum length in kmers (caller->min_ref_nkmers)
  kograph_filter_extend(caller->kograph,
                        nodes, unitig->num_nodes, forward,
                        caller->min_ref_nkmers, qoffset,
                        koruns, koruns_tmp, koruns_ended);

  size_t i, min_run_qoffset = SIZE_MAX, min_ended_run_qoffset = SIZE_MAX;
  size_t len, maxlen = 0;

  for(i = 0; i < koruns->len; i++) {
    min_run_qoffset = MIN2(min_run_qoffset, koruns->b[i].qoffset);
    len = korun_len(koruns->b[i]);
    maxlen = MAX2(maxlen, len);
  }

  // Stop if all our earliest runs have finished
  for(i = 0; i < koruns_ended->len; i++) {
    min_ended_run_qoffset = MIN2(min_ended_run_qoffset, koruns_ended->b[i].qoffset);
    len = korun_len(koruns->b[i]);
    maxlen = MAX2(maxlen, len);
  }

  // Continue if...
  return (min_run_qoffset <= min_ended_run_qoffset) &&
         (!caller->max_ref_nkmers || maxlen < caller->max_ref_nkmers);
}

// Try to pick up new runs at each unitig
static bool gcrawler_stop_at_ref_covg_path(const GraphCache *cache,
                                           const GCacheStep *step,
                                           void *arg)
{
  BreakpointCaller *caller = (BreakpointCaller*)arg;

  return gcrawler_stop_at_ref_covg(cache, step, caller,
                                   &caller->koruns_3p,
                                   &caller->koruns_tmp,
                                   &caller->koruns_3p_ended);
}

// For 5p flank only pick up new runs starting at the first unitig
static bool gcrawler_stop_at_ref_covg_flank5p(const GraphCache *cache,
                                              const GCacheStep *step,
                                              void *arg)
{
  BreakpointCaller *caller = (BreakpointCaller*)arg;

  bool cont = gcrawler_stop_at_ref_covg(cache, step, caller,
                                        &caller->koruns_5p,
                                        &caller->koruns_tmp,
                                        &caller->koruns_5p_ended);

  // size_t i;
  // KOccurRunBuffer *koruns = &caller->koruns_5p;
  // for(i = 0; i < koruns->len; i++)
  //   status("qoffset: %u", koruns->b[i].qoffset);

  return cont && (caller->koruns_5p.len > 0);
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
  korun_buf_capacity(runs_buf, runs_buf->len+koruns->len+koruns_ended->len);
  korun_buf_push(runs_buf, koruns_ended->b, koruns_ended->len);

  runs_buf->len += koruns_filter(runs_buf->b + runs_buf->len,
                                 caller->min_ref_nkmers,
                                 koruns->b, koruns->len);

  // status("Added %zu at end", runs_buf->len - init_len - koruns_ended->len);

  korun_buf_reset(koruns);
  korun_buf_reset(koruns_ended);

  ctx_assert(pathid < MAX_REFRUNS_PER_ORIENT(caller->db_graph->num_of_cols));

  ref_runs[pathid].first_runid = init_len;
  ref_runs[pathid].num_runs = runs_buf->len - init_len;
}

static void gcrawler_finish_ref_covg_path(const GraphCache *cache,
                                          uint32_t pathid,
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

static void gcrawler_finish_ref_covg_flank5p(const GraphCache *cache,
                                             uint32_t pathid,
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


static inline KOccurRun* fetch_ref_contact(const PathRefRun *ref_run,
                                           KOccurRunBuffer *runbuf,
                                           size_t *nruns)
{
  // Get runs along the ref and sort them
  *nruns = ref_run->num_runs;
  KOccurRun *koruns = runbuf->b + ref_run->first_runid;
  koruns_sort_by_qoffset(koruns, *nruns);
  return koruns;
}

// Traverse from node0 -> node1
static void traverse_5pflank(BreakpointCaller *caller, GraphCrawler *crawler,
                             dBNode node0, dBNode node1)
{
  const dBGraph *db_graph = crawler->cache.db_graph;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  size_t i, num_next;
  BinaryKmer bkmer0 = db_node_get_bkey(db_graph, node0.key);

  num_next = db_graph_next_nodes(db_graph, bkmer0, node0.orient,
                                 db_node_edges(db_graph, node0.key, 0),
                                 next_nodes, next_nucs);

  // Find index of previous node
  for(i = 0; i < num_next && !db_nodes_are_equal(next_nodes[i],node1); i++) {}

  ctx_assert(i < num_next && db_nodes_are_equal(next_nodes[i],node1));

  korun_buf_reset(&caller->koruns_5p);
  korun_buf_reset(&caller->koruns_5p_ended);
  korun_buf_reset(&caller->flank5p_run_buf);

  // Go backwards to get 5p flank
  // NULL means loop from 0..(ncols-1)
  graph_crawler_fetch(crawler, node0,
                      next_nodes, i, num_next,
                      NULL, db_graph->num_of_cols,
                      gcrawler_stop_at_ref_covg_flank5p,
                      gcrawler_finish_ref_covg_flank5p,
                      caller);
}

// Walk the graph remembering the last time we met the ref
// When traversal fails, dump sequence up to last meeting with the ref
static void follow_break(BreakpointCaller *caller, dBNode node)
{
  size_t i, j, k, num_next;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  size_t nonref_idx[4], num_nonref_next = 0;
  const dBGraph *db_graph = caller->db_graph;

  BinaryKmer bkey = db_node_get_bkey(db_graph, node.key);
  Edges edges = db_node_get_edges(db_graph, node.key, 0);

  num_next = db_graph_next_nodes(db_graph, bkey, node.orient, edges,
                                 next_nodes, next_nucs);

  // Filter out next nodes in the reference
  for(i = 0; i < num_next; i++) {
    if(!kograph_occurs(caller->kograph, next_nodes[i].key)) {
      nonref_idx[num_nonref_next++] = i;
    }
  }

  // Abandon if no non-ref kmers next
  if(num_nonref_next == 0) return;

  // debug
  // char nstr[MAX_KMER_SIZE+3];
  // db_node_to_str(db_graph, node, nstr);
  // status("follow %s", nstr);

  // Follow all paths not in ref, in all colours
  GraphCrawler *fw_crawler = &caller->crawlers[node.orient];
  GraphCrawler *rv_crawler = &caller->crawlers[!node.orient];
  dBNodeBuffer *allelebuf = &caller->allelebuf, *flank5pbuf = &caller->flank5pbuf;
  GCMultiColPath *flank5p_multicolpath, *allele_multicolpath;
  KOccurRun *flank5p_runs, *flank3p_runs;
  size_t flank5p_pathid, allele_pathid;
  size_t nflank5p_runs, nflank3p_runs;

  // We fetch 5' flanks in all colours then merge matching paths
  // we stop fetching a single path if it stops tracking the reference
  // Alternatively, we could fetch the 5' flank in everyone and stop after a
  // given distance, then check for that set of paths how much it tracks the
  // reference. This has the advantage of scaling much better with number of
  // samples, but not so well as min_ref_nkmers increases (since we fetch
  // many flanks that can't be used) - I think this is less of a worry.

  // Loop over possible next nodes at this junction
  for(i = 0; i < num_nonref_next; i++)
  {
    size_t next_idx = nonref_idx[i];

    // Go backwards to get 5p flank
    traverse_5pflank(caller, rv_crawler, db_node_reverse(next_nodes[next_idx]),
                     db_node_reverse(node));

    // if(!rv_crawler->num_paths) { status("No 5p"); }

    // Loop over the flanks that we got
    for(j = 0; j < rv_crawler->num_paths; j++)
    {
      // Get 5p flank
      db_node_buf_reset(flank5pbuf);
      graph_crawler_get_path_nodes(rv_crawler, j, flank5pbuf);
      flank5p_multicolpath = &rv_crawler->multicol_paths[j];
      flank5p_pathid = flank5p_multicolpath->pathid;

      // Fetch 5pflank ref position
      flank5p_runs = fetch_ref_contact(&caller->flank5p_refs[flank5p_pathid],
                                       &caller->flank5p_run_buf,
                                       &nflank5p_runs);

      koruns_reverse(flank5p_runs, nflank5p_runs, flank5pbuf->len);
      koruns_sort_by_qoffset(flank5p_runs, nflank5p_runs);
      db_nodes_reverse_complement(flank5pbuf->b, flank5pbuf->len);

      // if(!nflank5p_runs) { status("  No 5p (b)"); }

      if(nflank5p_runs > 0)
      {
        // Reset caller
        korun_buf_reset(&caller->koruns_3p);
        korun_buf_reset(&caller->koruns_3p_ended);
        korun_buf_reset(&caller->allele_run_buf);

        // functions gcrawler_stop_at_ref_covg_path(),
        //           gcrawler_finish_ref_covg_path()
        // both fill koruns_3p, koruns_3p_ended and allele_run_buf

        // Only traverse in the colours we have a flank for
        graph_crawler_fetch(fw_crawler, node,
                            next_nodes, next_idx, num_next,
                            flank5p_multicolpath->cols,
                            flank5p_multicolpath->num_cols,
                            gcrawler_stop_at_ref_covg_path,
                            gcrawler_finish_ref_covg_path,
                            caller);

        // Assemble contigs - fetch forwards for each path for given 5p flank
        for(k = 0; k < fw_crawler->num_paths; k++)
        {
          // Fetch nodes
          db_node_buf_reset(allelebuf);
          graph_crawler_get_path_nodes(fw_crawler, k, allelebuf);
          ctx_assert(allelebuf->len > 0);

          allele_multicolpath = &fw_crawler->multicol_paths[k];
          allele_pathid = allele_multicolpath->pathid;

          // Fetch 3pflank ref position
          flank3p_runs = fetch_ref_contact(&caller->allele_refs[allele_pathid],
                                           &caller->allele_run_buf,
                                           &nflank3p_runs);

          process_contig(caller,
                         allele_multicolpath->cols,
                         allele_multicolpath->num_cols,
                         flank5pbuf, allelebuf,
                         flank5p_runs, nflank5p_runs,
                         flank3p_runs, nflank3p_runs);
        }
      }
    }
  }
}

static inline int breakpoint_caller_node(hkey_t hkey, BreakpointCaller *caller)
{
  // DEBUG
  // const dBGraph *db_graph = caller->db_graph;
  // char kstr[MAX_KMER_SIZE+1];
  // BinaryKmer bkmer = db_node_get_bkey(db_graph, hkey);
  // binary_kmer_to_str(bkmer, db_graph->kmer_size, kstr);
  // if(strcmp(kstr,"GTTGCTCATGA")) return 0; // skip all but given kmer
  // status("brk %s\n", kstr);
  //

  // check node is in the ref
  if(kograph_occurs(caller->kograph, hkey))
  {
    graph_crawler_reset(&caller->crawlers[0]);
    graph_crawler_reset(&caller->crawlers[1]);
    follow_break(caller, (dBNode){.key = hkey, .orient = FORWARD});
    follow_break(caller, (dBNode){.key = hkey, .orient = REVERSE});
  }

  return 0; // => keep iterating
}

static void breakpoint_caller(void *ptr, size_t threadid)
{
  BreakpointCaller *caller = (BreakpointCaller*)ptr;
  ctx_assert(caller->db_graph->num_edge_cols == 1);

  HASH_ITERATE_PART(&caller->db_graph->ht, threadid, caller->nthreads,
                    breakpoint_caller_node, caller);
}

// Print JSON header to gzout
static void breakpoints_print_header(gzFile gzout, const char *out_path,
                                     char **seq_paths, size_t nseq_paths,
                                     const read_t *reads, size_t nreads,
                                     bool load_ref_edges,
                                     size_t min_ref_nkmers,
                                     size_t max_ref_nkmers,
                                     cJSON **hdrs, size_t nhdrs,
                                     size_t ref_col,
                                     const dBGraph *db_graph)
{
  size_t i;
  ctx_assert(nseq_paths > 0);

  // Construct cJSON
  cJSON *json = cJSON_CreateObject();

  cJSON_AddStringToObject(json, "file_format", "CtxBreakpoints");
  cJSON_AddNumberToObject(json, "format_version", BREAKPOINT_FORMAT_VERSION);

  // Add standard cortex headers
  json_hdr_make_std(json, out_path, hdrs, nhdrs, db_graph,
                    hash_table_nkmers(&db_graph->ht));

  // Update reference colour
  // json.graph.colours[ref_col]
  cJSON *hdr_graph  = json_hdr_get(json,      "graph",   cJSON_Object, out_path);
  cJSON *hdr_colour = json_hdr_get(hdr_graph, "colours", cJSON_Array,  out_path);
  cJSON *ref = cJSON_GetArrayItem(hdr_colour, ref_col);
  ctx_assert(ref != NULL);
  cJSON_AddTrueToObject(ref, "is_ref");

  // Add parameters used in bubble calling to the header
  json_hdr_augment_cmd(json, "breakpoints", "min_ref_flank_kmers",
                                            cJSON_CreateInt(min_ref_nkmers));
  json_hdr_augment_cmd(json, "breakpoints", "max_ref_flank_kmers",
                                            cJSON_CreateInt(max_ref_nkmers));
  json_hdr_augment_cmd(json, "breakpoints", "load_ref_edges",
                                            cJSON_CreateBool(load_ref_edges));

  // Add paths to reference files
  cJSON *ref_files = cJSON_CreateArray();
  for(i = 0; i < nseq_paths; i++)
  {
    // Get absolute path to output file if possible
    char abspath[PATH_MAX + 1];
    char *ref_path = realpath(seq_paths[i], abspath) ? abspath : seq_paths[i];
    cJSON_AddItemToArray(ref_files, cJSON_CreateString(ref_path));
  }
  json_hdr_augment_cmd(json, "breakpoints", "ref_files", ref_files);

  // List contigs
  cJSON *contigs = cJSON_CreateArray();
  for(i = 0; i < nreads; i++) {
    cJSON *contig = cJSON_CreateObject();
    cJSON_AddStringToObject(contig, "id", reads[i].name.b);
    cJSON_AddNumberToObject(contig, "length", reads[i].seq.end);
    cJSON_AddItemToArray(contigs, contig);
  }
  json_hdr_augment_cmd(json, "breakpoints", "contigs", contigs);

  // Write header to file
  json_hdr_gzprint(json, gzout);

  // Print comments about the format
  gzputs(gzout, "\n");
  gzputs(gzout, "# This file was generated with McCortex\n");
  gzputs(gzout, "#   written by Isaac Turner <turner.isaac@gmail.com>\n");
  gzputs(gzout, "#   url: "MCCORTEX_URL"\n");
  gzputs(gzout, "# \n");
  gzputs(gzout, "# Comment lines begin with a # and are ignored, but must come after the header\n");
  gzputs(gzout, "# Format is:\n");
  gzputs(gzout, "#   chr=seq:start-end:strand:offset\n");
  gzputs(gzout, "#   all coordinates are 1-based\n");
  gzputs(gzout, "#   <strand> is + or -. If +, start <= end. If -, start >= end.\n");
  gzputs(gzout, "#   <offset> is the position in the sequence where ref starts agreeing\n");
  gzputs(gzout, "\n");

  cJSON_Delete(json);
}

void breakpoints_call(size_t nthreads, size_t ref_col,
                      gzFile gzout, const char *out_path,
                      const read_t *reads, size_t num_reads,
                      char **seq_paths, size_t num_seq_paths,
                      bool load_ref_edges,
                      size_t min_ref_nkmers, size_t max_ref_nkmers,
                      cJSON **hdrs, size_t nhdrs,
                      dBGraph *db_graph)
{
  ctx_assert(!max_ref_nkmers || min_ref_nkmers <= max_ref_nkmers);
  // Temporarily hide edges from kograph_create if we don't want to load edges
  Edges *tmp_edges = db_graph->col_edges;
  if(!load_ref_edges) db_graph->col_edges = NULL;

  KOGraph kograph = kograph_create(reads, num_reads, true, ref_col,
                                   nthreads, db_graph);

  // Restore graph edges
  db_graph->col_edges = tmp_edges;

  BreakpointCaller *callers = brkpt_callers_new(nthreads, gzout,
                                                min_ref_nkmers, max_ref_nkmers,
                                                &kograph, db_graph);

  status("Running BreakpointCaller with %zu thread%s, output to: %s",
         nthreads, util_plural_str(nthreads),
         futil_outpath_str(out_path));

  status("  Finding breakpoints after at least %zu kmers (%zubp) of homology",
         min_ref_nkmers, min_ref_nkmers+db_graph->kmer_size-1);

  breakpoints_print_header(gzout, out_path,
                           seq_paths, num_seq_paths,
                           reads, num_reads,
                           load_ref_edges,
                           min_ref_nkmers, min_ref_nkmers,
                           hdrs, nhdrs,
                           ref_col, db_graph);

  util_run_threads(callers, nthreads, sizeof(callers[0]),
                   nthreads, breakpoint_caller);

  char call_num_str[100];
  ulong_to_str(callers[0].callid[0], call_num_str);
  status("  %s calls printed to %s", call_num_str, futil_outpath_str(out_path));

  brkpt_callers_destroy(callers, nthreads);
  kograph_dealloc(&kograph);
}
