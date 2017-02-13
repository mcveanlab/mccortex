#include "global.h"
#include "bubble_caller.h"
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"
#include "cmd.h"
#include "seq_reader.h"
#include "db_unitig.h"
#include "binary_seq.h"
#include "graph_crawler.h"
#include "json_hdr.h"

#include <pthread.h> // multithreading

BubbleCaller* bubble_callers_new(size_t num_callers,
                                 const BubbleCallingPrefs *prefs,
                                 gzFile gzout,
                                 const dBGraph *db_graph)
{
  ctx_assert(num_callers > 0);

  // Max usage is 4 * max_allele_len * cols
  size_t i;
  size_t max_path_len = MAX2(prefs->max_flank_len, prefs->max_allele_len);

  BubbleCaller *callers = ctx_malloc(num_callers * sizeof(BubbleCaller));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  uint64_t *nbubbles_ptr = ctx_calloc(1, sizeof(uint64_t));

  for(i = 0; i < num_callers; i++)
  {
    bool *haploid_seen = ctx_calloc(prefs->nhaploid_cols, sizeof(bool));

    BubbleCaller tmp = {.nthreads = num_callers,
                        .haploid_seen = haploid_seen,
                        .num_haploid_bubbles = 0,
                        .num_serial_bubbles = 0,
                        .nbubbles_ptr = nbubbles_ptr,
                        .prefs = prefs,
                        .db_graph = db_graph, .gzout = gzout,
                        .out_lock = out_lock};

    memcpy(&callers[i], &tmp, sizeof(BubbleCaller));

    callers[i].unitig_map = kh_init(uint32to32);

    // First two buffers don't actually need to grow
    db_node_buf_alloc(&callers[i].flank5p, prefs->max_flank_len);
    db_node_buf_alloc(&callers[i].pathbuf, max_path_len);

    graph_walker_alloc(&callers[i].wlk, db_graph);
    rpt_walker_alloc(&callers[i].rptwlk, db_graph->ht.capacity, 22); // 4MB

    graph_cache_alloc(&callers[i].cache, db_graph);
    cache_stepptr_buf_alloc(&callers[i].spp_forward, 1024);
    cache_stepptr_buf_alloc(&callers[i].spp_reverse, 1024);
    strbuf_alloc(&callers[i].output_buf, 2048);
  }

  return callers;
}

void bubble_callers_destroy(BubbleCaller *callers, size_t num_callers)
{
  ctx_assert(num_callers > 0);

  size_t i;
  for(i = 0; i < num_callers; i++)
  {
    ctx_free(callers[i].haploid_seen);

    kh_destroy(uint32to32, callers[i].unitig_map);

    db_node_buf_dealloc(&callers[i].flank5p);
    db_node_buf_dealloc(&callers[i].pathbuf);

    rpt_walker_dealloc(&callers[i].rptwlk);
    graph_walker_dealloc(&callers[i].wlk);

    graph_cache_dealloc(&callers[i].cache);
    cache_stepptr_buf_dealloc(&callers[i].spp_forward);
    cache_stepptr_buf_dealloc(&callers[i].spp_reverse);
    strbuf_dealloc(&callers[i].output_buf);
  }
  pthread_mutex_destroy(callers[0].out_lock);
  ctx_free(callers[0].out_lock);
  ctx_free(callers[0].nbubbles_ptr);
  ctx_free(callers);
}

// Print JSON header to gzout
static void bubble_caller_print_header(gzFile gzout, const char* out_path,
                                       const BubbleCallingPrefs *prefs,
                                       cJSON **hdrs, size_t nhdrs,
                                       const dBGraph *db_graph)
{
  size_t i;

  // Construct cJSON
  cJSON *jsonhdr = cJSON_CreateObject();

  cJSON_AddStringToObject(jsonhdr, "file_format", "CtxBubbles");
  cJSON_AddNumberToObject(jsonhdr, "format_version", BUBBLE_FORMAT_VERSION);

  // Add standard cortex headers
  json_hdr_make_std(jsonhdr, out_path, hdrs, nhdrs, db_graph,
                    hash_table_nkmers(&db_graph->ht));

  // Add parameters used in bubble calling to the header
  json_hdr_augment_cmd(jsonhdr, "bubbles", "max_flank_kmers",  cJSON_CreateInt(prefs->max_flank_len));
  json_hdr_augment_cmd(jsonhdr, "bubbles", "max_allele_kmers", cJSON_CreateInt(prefs->max_allele_len));
  cJSON *haploids = cJSON_CreateArray();
  for(i = 0; i < prefs->nhaploid_cols; i++)
    cJSON_AddItemToArray(haploids, cJSON_CreateInt(prefs->haploid_cols[i]));
  json_hdr_augment_cmd(jsonhdr, "bubbles", "haploid_colours", haploids);

  // Write header to file
  json_hdr_gzprint(jsonhdr, gzout);

  // Print comments about the format
  gzputs(gzout, "\n");
  gzputs(gzout, "# This file was generated with McCortex\n");
  gzputs(gzout, "#   written by Isaac Turner <turner.isaac@gmail.com>\n");
  gzputs(gzout, "#   url: "MCCORTEX_URL"\n");
  gzputs(gzout, "# \n");
  gzputs(gzout, "# Comment lines begin with a # and are ignored, but must come after the header\n");
  gzputs(gzout, "\n");

  cJSON_Delete(jsonhdr);
}

static void branch_to_str(const dBNode *nodes, size_t len, bool print_first_kmer,
                          StrBuf *sbuf, const dBGraph *db_graph)
{
  size_t i = print_first_kmer, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;

  if(print_first_kmer) {
    strbuf_ensure_capacity(sbuf, sbuf->end + kmer_size);
    bkmer = db_node_oriented_bkmer(db_graph, nodes[0]);
    binary_kmer_to_str(bkmer, kmer_size, sbuf->b+sbuf->end);
    sbuf->end += kmer_size;
  }

  // i == 1 if print_first_kmer, otherwise 0
  strbuf_ensure_capacity(sbuf, sbuf->end + len + 1); // +1 for '\n'
  for(; i < len; i++) {
    nuc = db_node_get_last_nuc(nodes[i], db_graph);
    sbuf->b[sbuf->end++] = dna_nuc_to_char(nuc);
  }

  sbuf->b[sbuf->end++] = '\n';
  sbuf->b[sbuf->end] = '\0';
}


// Potential bubble - filter ref and duplicate alleles
static void print_bubble(BubbleCaller *caller,
                         GCacheStep **steps, size_t num_paths)
{
  const BubbleCallingPrefs *prefs = caller->prefs;
  const dBGraph *db_graph = caller->db_graph;
  size_t i;

  dBNodeBuffer *flank5p = &caller->flank5p;
  if(flank5p->len == 0)
  {
    // Haven't fetched 5p flank yet
    // flank5p[0] already contains the first node
    flank5p->len = 1;
    db_unitig_extend(flank5p, prefs->max_flank_len, db_graph);
    db_nodes_reverse_complement(flank5p->b, flank5p->len);
  }

  //
  // Print Bubble
  //

  // write to string buffer then flush to gzFile
  StrBuf *sbuf = &caller->output_buf;
  strbuf_reset(sbuf);

  // Temporary node buffer to use
  dBNodeBuffer *pathbuf = &caller->pathbuf;
  db_node_buf_reset(pathbuf);

  // Get bubble number (threadsafe nbubbles_ptr++)
  size_t id = __sync_fetch_and_add((volatile uint64_t*)caller->nbubbles_ptr, 1);

  // This can be set to anything without a '.' in it
  const char prefix[] = "call";

  // 5p flank
  // strbuf_sprintf(sbuf, ">bubble.%s%zu.5pflank kmers=%zu\n", prefix, id, flank5p->len);
  strbuf_append_str(sbuf, ">bubble.");
  strbuf_append_str(sbuf, prefix);
  strbuf_append_ulong(sbuf, id);
  strbuf_append_str(sbuf, ".5pflank kmers=");
  strbuf_append_ulong(sbuf, flank5p->len);
  strbuf_append_char(sbuf, '\n');
  branch_to_str(flank5p->b, flank5p->len, true, sbuf, db_graph);

  // 3p flank
  db_node_buf_reset(pathbuf);
  const GCacheUnitig *unitig = gc_step_get_unitig(&caller->cache, steps[0]);
  gc_db_unitig_fetch_nodes(&caller->cache, unitig, steps[0]->orient, pathbuf);

  // strbuf_sprintf(sbuf, ">bubble.%s%zu.3pflank kmers=%zu\n", prefix, id, pathbuf->len);
  strbuf_append_str(sbuf, ">bubble.");
  strbuf_append_str(sbuf, prefix);
  strbuf_append_ulong(sbuf, id);
  strbuf_append_str(sbuf, ".3pflank kmers=");
  strbuf_append_ulong(sbuf, pathbuf->len);
  strbuf_append_char(sbuf, '\n');
  branch_to_str(pathbuf->b, pathbuf->len, false, sbuf, db_graph);

  // Print alleles
  for(i = 0; i < num_paths; i++)
  {
    db_node_buf_reset(pathbuf);
    gc_step_fetch_nodes(&caller->cache, steps[i], pathbuf);

    // strbuf_sprintf(sbuf, ">bubble.%s%zu.branch.%zu kmers=%zu\n",
    //                prefix, id, i, pathbuf->len);
    strbuf_append_str(sbuf, ">bubble.");
    strbuf_append_str(sbuf, prefix);
    strbuf_append_ulong(sbuf, id);
    strbuf_append_str(sbuf, ".branch.");
    strbuf_append_ulong(sbuf, i);
    strbuf_append_str(sbuf, " kmers=");
    strbuf_append_ulong(sbuf, pathbuf->len);
    strbuf_append_char(sbuf, '\n');

    branch_to_str(pathbuf->b, pathbuf->len, false, sbuf, db_graph);
  }

  strbuf_append_char(sbuf, '\n');

  ctx_assert(strlen(sbuf->b) == sbuf->end);

  // lock, print, unlock
  pthread_mutex_lock(caller->out_lock);
  gzwrite(caller->gzout, sbuf->b, sbuf->end);
  pthread_mutex_unlock(caller->out_lock);
}

// `fork_node` is a node with outdegree > 1
void find_bubbles(BubbleCaller *caller, dBNode fork_node)
{
  graph_cache_reset(&caller->cache);

  const dBGraph *db_graph = caller->db_graph;
  GraphCache *cache = &caller->cache;
  GraphWalker *wlk = &caller->wlk;
  RepeatWalker *rptwlk = &caller->rptwlk;

  // char tmpstr[MAX_KMER_SIZE+3];
  // db_node_to_str(db_graph, fork_node, tmpstr);
  // status("Calling from %s", tmpstr);

  dBNode nodes[4];
  Nucleotide bases[4];
  size_t i, num_next, num_edges_in_col;
  BinaryKmer fork_bkmer = db_node_get_bkey(db_graph, fork_node.key);

  num_next = db_graph_next_nodes(db_graph, fork_bkmer, fork_node.orient,
                                 db_node_edges(db_graph, fork_node.key, 0),
                                 nodes, bases);

  // loop over alleles, then colours
  Colour colour, colours_loaded = db_graph->num_of_cols;
  bool node_has_col[4];

  uint32_t pathid;

  for(colour = 0; colour < colours_loaded; colour++)
  {
    if(!db_node_has_col(db_graph, fork_node.key, colour)) continue;

    // Determine if this fork is a fork in the current colour
    num_edges_in_col = 0;
    for(i = 0; i < num_next; i++) {
      node_has_col[i] = (db_node_has_col(db_graph, nodes[i].key, colour) > 0);
      num_edges_in_col += node_has_col[i];
    }

    graph_walker_setup(wlk, true, colour, colour, db_graph);

    for(i = 0; i < num_next; i++)
    {
      if(node_has_col[i])
      {
        graph_walker_start(wlk, fork_node);
        graph_walker_force(wlk, nodes[i], num_edges_in_col > 1);

        pathid = graph_crawler_load_path_limit(cache, nodes[i], wlk, rptwlk,
                                               caller->prefs->max_allele_len);

        graph_walker_finish(wlk);
        graph_crawler_reset_rpt_walker(rptwlk, cache, pathid);
      }
    }
  }

  // Set up 5p flank
  caller->flank5p.b[0] = db_node_reverse(fork_node);
  caller->flank5p.len = 0; // set to one to signify we haven't fetched flank yet
}

static bool paths_all_share_unitig(const GraphCache *cache,
                                   khash_t(uint32to32) *unitig_map,
                                   GCacheStep const*const* steps,
                                   size_t num_paths)
{
  uint32_t unitig;
  khiter_t kiter;
  kh_clear(uint32to32, unitig_map);
  const GCachePath *path;
  const GCacheStep *step;
  size_t i;
  int hret = 0;

  for(i = 0; i < num_paths; i++) {
    path = gc_step_get_path(cache, steps[i]);
    step = gc_path_first_step(cache, path);
    for(; step < steps[i]; step++) {
      unitig = gc_step_encode_uint32(step);
      kiter = kh_put(uint32to32, unitig_map, unitig, &hret);
      if(hret < 0) die("khash table failed: out of memory?");
      if(hret > 0) kh_value(unitig_map, kiter) = 0; // init if not in table
      kh_value(unitig_map, kiter)++;
    }
  }

  // Look for hits and wipe at the same time
  bool shared_unitig = false;
  for(kiter = kh_begin(unitig_map); kiter != kh_end(unitig_map); ++kiter) {
    if(kh_exist(unitig_map, kiter)) {
      shared_unitig |= (kh_value(unitig_map, kiter) == num_paths);
      kh_del(uint32to32, unitig_map, kiter);
    }
  }

  return shared_unitig;
}

// Remove paths that are both seen in a haploid sample (e.g. repeat)
// Returns number of paths
static size_t remove_haploid_paths(const GraphCache *cache,
                                   GCacheStep **steps, size_t num_paths,
                                   bool *haploid_seen,
                                   const size_t *haploid_cols,
                                   size_t nhaploid_cols)
{
  size_t r, p;
  memset(haploid_seen, 0, sizeof(bool)*nhaploid_cols);

  for(p = 0; p < num_paths; )
  {
    for(r = 0; r < nhaploid_cols; r++)
    {
      if(graph_cache_step_has_colour(cache, steps[p], haploid_cols[r]))
      {
        // Drop path if already haploid_seen
        if(haploid_seen[r]) break;
        haploid_seen[r] = 1;
      }
    }

    // Drop path
    if(r < nhaploid_cols) {
      SWAP(steps[p], steps[num_paths-1]);
      num_paths--;
    }
    else p++;
  }

  return num_paths;
}

// returns true if paths contain a bubble after filtering, otherwise false
static bool filter_bubbles(BubbleCaller *bc, GCacheStepPtrBuf *ends)
{
  if(!graph_cache_is_3p_flank(&bc->cache, ends->b, ends->len)) {
    // status("fail: Not 3p");
    return false;
  }

  if((ends->len = graph_cache_remove_dupes(&bc->cache, ends->b, ends->len)) < 2) {
    // status("fail: all dupes");
    return false;
  }

  if((ends->len = remove_haploid_paths(&bc->cache,
                                       ends->b, ends->len,
                                       bc->haploid_seen,
                                       bc->prefs->haploid_cols,
                                       bc->prefs->nhaploid_cols)) < 2)
  {
    // status("fail: haploid");
    bc->num_haploid_bubbles++; // haploid bubble removed
    return false;
  }

  // remove serial bubbles by dropping all paths if they all share a unitig
  if(bc->prefs->remove_serial_bubbles &&
     paths_all_share_unitig(&bc->cache, bc->unitig_map,
                            (GCacheStep const*const*)ends->b, ends->len))
  {
    // status("fail: serial");
    bc->num_serial_bubbles++;
    return false;
  }

  return true;
}

// Load GCacheSteps into caller->spp_forward (if they traverse the unitig forward)
// or caller->spp_reverse (if they traverse the unitig in reverse)
void find_bubbles_ending_with(BubbleCaller *bc, GCacheUnitig *unitig)
{
  // possible 3p flank (i.e. bubble end)
  // record paths that go through here forwards, and in reverse
  cache_stepptr_buf_reset(&bc->spp_forward);
  cache_stepptr_buf_reset(&bc->spp_reverse);

  uint32_t stepid = unitig->stepid;
  GCacheStep *step;

  while(stepid != UINT32_MAX)
  {
    step = graph_cache_step(&bc->cache, stepid);
    if(step->orient == FORWARD) {
      cache_stepptr_buf_push(&bc->spp_forward, &step, 1);
    } else {
      cache_stepptr_buf_push(&bc->spp_reverse, &step, 1);
    }
    stepid = step->next_step;
  }

  // Filter out non-bubbles
  if(!filter_bubbles(bc, &bc->spp_forward)) cache_stepptr_buf_reset(&bc->spp_forward);
  if(!filter_bubbles(bc, &bc->spp_reverse)) cache_stepptr_buf_reset(&bc->spp_reverse);
}

static void write_bubbles_to_file(BubbleCaller *caller)
{
  // Loop over unitigs checking if they are 3p flanks
  size_t nunitigs = graph_cache_num_unitigs(&caller->cache);
  GCacheUnitig *unitig;
  size_t i;

  // status("Got %zu nunitigs", nunitigs);

  for(i = 0; i < nunitigs; i++)
  {
    unitig = graph_cache_unitig(&caller->cache, i);
    find_bubbles_ending_with(caller, unitig);

    // status("ends: %zu %zu", caller->spp_forward.len, caller->spp_reverse.len);

    if(caller->spp_forward.len > 1)
      print_bubble(caller, caller->spp_forward.b, caller->spp_forward.len);
    if(caller->spp_reverse.len > 1)
      print_bubble(caller, caller->spp_reverse.b, caller->spp_reverse.len);
  }
}

static inline int bubble_caller_node(hkey_t hkey, BubbleCaller *caller)
{
  Edges edges = db_node_get_edges(caller->db_graph, hkey, 0);
  if(edges_get_outdegree(edges, FORWARD) > 1) {
    find_bubbles(caller, (dBNode){.key = hkey, .orient = FORWARD});
    write_bubbles_to_file(caller);
  }
  if(edges_get_outdegree(edges, REVERSE) > 1) {
    find_bubbles(caller, (dBNode){.key = hkey, .orient = REVERSE});
    write_bubbles_to_file(caller);
  }

  return 0; // => keep iterating
}

void bubble_caller(void *args, size_t threadid)
{
  BubbleCaller *caller = (BubbleCaller*)args;

  HASH_ITERATE_PART(&caller->db_graph->ht, threadid, caller->nthreads,
                    bubble_caller_node, caller);
}

void invoke_bubble_caller(size_t num_of_threads,
                          const BubbleCallingPrefs *prefs,
                          gzFile gzout, const char *out_path,
                          cJSON **hdrs, size_t nhdrs,
                          const dBGraph *db_graph)
{
  ctx_assert(db_graph->num_edge_cols == 1);
  ctx_assert(db_graph->node_in_cols != NULL);
  size_t i;

  status("Calling bubbles with %zu threads, output: %s", num_of_threads, out_path);

  StrBuf tmpstr = {0,0,0};
  for(i = 0; i < prefs->nhaploid_cols; i++)
    strbuf_sprintf(&tmpstr, "\t%zu", prefs->haploid_cols[i]);
  status("Haploid colours:%s", tmpstr.b);
  strbuf_dealloc(&tmpstr);

  // Print header
  bubble_caller_print_header(gzout, out_path, prefs, hdrs, nhdrs, db_graph);

  BubbleCaller *callers = bubble_callers_new(num_of_threads, prefs,
                                             gzout, db_graph);

  // Run
  util_run_threads(callers, num_of_threads, sizeof(callers[0]),
                   num_of_threads, bubble_caller);

  // Report number of bubble called+printed
  uint64_t nhaploid = 0, nserial = 0, nbubbles = callers[0].nbubbles_ptr[0];

  for(i = 0; i < num_of_threads; i++) {
    nhaploid += callers[i].num_haploid_bubbles;
    nserial += callers[i].num_serial_bubbles;
  }

  char n0[ULONGSTRLEN];
  status("Bubble Caller called %s bubbles\n", ulong_to_str(nbubbles, n0));
  status("Haploid bubbles dropped: %s", ulong_to_str(nhaploid, n0));
  status("Serial bubbles dropped: %s", ulong_to_str(nserial, n0));

  status("Turn bubble file into VCF with:");
  status("   bwa index ref.fa");
  status("   scripts/cortex_print_flanks.sh %s > %s.flanks", out_path, out_path);
  status("   bwa mem ref.fa %s.flanks > %s.sam", out_path, out_path);
  status("   ctx%i calls2vcf -F %s.sam -o output.vcf %s ref.fa",
         (int)get_max_kmer_size(), out_path, out_path);
  timestamp(); message("\n");

  // Clean up
  bubble_callers_destroy(callers, num_of_threads);
}
