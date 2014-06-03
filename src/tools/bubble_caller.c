#include "global.h"
#include "bubble_caller.h"
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"
#include "cmd.h"
#include "seq_reader.h"
#include "path_store.h"
#include "supernode.h"
#include "binary_seq.h"
#include "graph_crawler.h"
#include "caller_output.h"

#include <time.h> // printing datetime
#include <pthread.h> // multithreading

BubbleCaller* bubble_callers_new(size_t num_callers,
                                 BubbleCallingPrefs prefs,
                                 gzFile gzout,
                                 const dBGraph *db_graph)
{
  ctx_assert(num_callers > 0);

  // Max usage is 4 * max_allele_len * cols
  size_t i;
  size_t max_path_len = MAX2(prefs.max_flank_len, prefs.max_allele_len);

  BubbleCaller *callers = ctx_malloc(num_callers * sizeof(BubbleCaller));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  size_t *num_bubbles_ptr = ctx_calloc(1, sizeof(size_t));

  for(i = 0; i < num_callers; i++)
  {
    BubbleCaller tmp = {.threadid = i, .nthreads = num_callers,
                        .haploid_seen = ctx_calloc(1+prefs.num_haploid, sizeof(bool)),
                        .num_bubbles_ptr = num_bubbles_ptr,
                        .prefs = prefs,
                        .db_graph = db_graph, .gzout = gzout,
                        .out_lock = out_lock};
  
    memcpy(&callers[i], &tmp, sizeof(BubbleCaller));

    // First two buffers don't actually need to grow
    db_node_buf_alloc(&callers[i].flank5p, prefs.max_flank_len);
    db_node_buf_alloc(&callers[i].pathbuf, max_path_len);

    graph_walker_alloc(&callers[i].wlk);
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
  ctx_free(callers[0].num_bubbles_ptr);
  ctx_free(callers);
}


static void bubble_caller_print_header(gzFile gzout, const char* out_path,
                                       const dBGraph *db_graph)
{
  caller_gzprint_header(gzout, out_path, "CtxBubblesv0.1", db_graph);
}

static void branch_to_str(const dBNode *nodes, size_t len, bool print_first_kmer,
                          StrBuf *sbuf, const dBGraph *db_graph)
{
  size_t i = print_first_kmer, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;

  if(print_first_kmer) {
    strbuf_ensure_capacity(sbuf, sbuf->len + kmer_size);
    bkmer = db_node_oriented_bkmer(db_graph, nodes[0]);
    binary_kmer_to_str(bkmer, kmer_size, sbuf->buff+sbuf->len);
    sbuf->len += kmer_size;
  }

  // i == 1 if print_first_kmer, otherwise 0
  strbuf_ensure_capacity(sbuf, sbuf->len + len + 1); // +1 for '\n'
  for(; i < len; i++) {
    nuc = db_node_get_last_nuc(nodes[i], db_graph);
    sbuf->buff[sbuf->len++] = dna_nuc_to_char(nuc);
  }

  sbuf->buff[sbuf->len++] = '\n';
  sbuf->buff[sbuf->len] = '\0';
}

// Remove paths that are both seen in a haploid sample (e.g. repeat)
// Returns number of paths
static size_t remove_haploid_paths(const GraphCache *cache,
                                   GCacheStep **steps, size_t num_paths,
                                   bool *haploid_seen,
                                   const size_t *haploid_cols,
                                   size_t num_haploid)
{
  size_t r, p;
  memset(haploid_seen, 0, sizeof(bool)*num_haploid);

  for(p = 0; p < num_paths; )
  {
    for(r = 0; r < num_haploid; r++)
    {
      if(graph_cache_step_has_colour(cache, steps[p], haploid_cols[r]))
      {
        // Drop path if already haploid_seen
        if(haploid_seen[r]) break;
        haploid_seen[r] = 1;
      }
    }

    // Drop path
    if(r < num_haploid) {
      SWAP(steps[p], steps[num_paths-1]);
      num_paths--;
    }
    else p++;
  }

  return num_paths;
}


// Potential bubble - filter ref and duplicate alleles
static void print_bubble(BubbleCaller *caller,
                         GCacheStep **steps, size_t num_paths)
{
  const BubbleCallingPrefs prefs = caller->prefs;
  const dBGraph *db_graph = caller->db_graph;
  GCacheSnode *snode;
  size_t i;

  dBNodeBuffer *flank5p = &caller->flank5p;
  if(flank5p->len == 0)
  {
    // Haven't fetched 5p flank yet
    // flank5p[0] already contains the first node
    flank5p->len = 1;
    supernode_extend(flank5p, prefs.max_flank_len, db_graph);
    db_nodes_reverse_complement(flank5p->data, flank5p->len);
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

  // Get bubble number (threadsafe num_bubbles_ptr++)
  size_t id = __sync_fetch_and_add((volatile size_t*)caller->num_bubbles_ptr, 1);

  // 5p flank
  strbuf_sprintf(sbuf, ">bubble.%zu.5pflank kmers=%zu\n", id, flank5p->len);
  branch_to_str(flank5p->data, flank5p->len, true, sbuf, db_graph);

  // 3p flank
  db_node_buf_reset(pathbuf);
  snode = graph_cache_snode(&caller->cache, steps[0]->supernode);
  graph_cache_snode_fetch_nodes(&caller->cache, snode, steps[0]->orient, pathbuf);

  strbuf_sprintf(sbuf, ">bubble.%zu.3pflank kmers=%zu\n", id, pathbuf->len);
  branch_to_str(pathbuf->data, pathbuf->len, false, sbuf, db_graph);

  // Print alleles
  for(i = 0; i < num_paths; i++)
  {
    db_node_buf_reset(pathbuf);
    graph_cache_step_fetch_nodes(&caller->cache, steps[i], pathbuf);
    strbuf_sprintf(sbuf, ">bubble.%zu.branch.%zu kmers=%zu\n", id, i, pathbuf->len);
    branch_to_str(pathbuf->data, pathbuf->len, false, sbuf, db_graph);
  }

  strbuf_append_char(sbuf, '\n');

  ctx_assert(strlen(sbuf->buff) == sbuf->len);

  // lock, print, unlock
  pthread_mutex_lock(caller->out_lock);
  gzputs(caller->gzout, sbuf->buff);
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
  BinaryKmer fork_bkmer = db_node_get_bkmer(db_graph, fork_node.key);

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

    for(i = 0; i < num_next; i++)
    {
      if(node_has_col[i])
      {
        graph_walker_init(wlk, db_graph, colour, colour, fork_node);
        graph_walker_force(wlk, nodes[i].key, bases[i], num_edges_in_col > 1);

        pathid = graph_crawler_load_path_limit(cache, nodes[i], wlk, rptwlk,
                                               caller->prefs.max_allele_len);

        graph_walker_finish(wlk);
        graph_crawler_reset_rpt_walker(rptwlk, cache, pathid);
      }
    }
  }

  // Set up 5p flank
  caller->flank5p.data[0] = db_node_reverse(fork_node);
  caller->flank5p.len = 0; // set to one to signify we haven't fetched flank yet
}

static void remove_non_bubbles(BubbleCaller *caller, GCacheStepPtrBuf *endsteps)
{
  if(!graph_cache_is_3p_flank(&caller->cache, endsteps->data, endsteps->len)) {
    cache_stepptr_buf_reset(endsteps);
  }
  else
  {
    endsteps->len = graph_cache_remove_dupes(&caller->cache,
                                             endsteps->data, endsteps->len);

    endsteps->len = remove_haploid_paths(&caller->cache,
                                         endsteps->data, endsteps->len,
                                         caller->haploid_seen,
                                         caller->prefs.haploid_cols,
                                         caller->prefs.num_haploid);

    if(endsteps->len < 2) cache_stepptr_buf_reset(endsteps);
  }
}

// Load GCacheSteps into caller->spp_forward (if they traverse the snode forward)
// or caller->spp_reverse (if they traverse the snode in reverse)
void find_bubbles_ending_with(BubbleCaller *caller, GCacheSnode *snode)
{
  // possible 3p flank (i.e. bubble end)
  // record paths that go through here forwards, and in reverse
  cache_stepptr_buf_reset(&caller->spp_forward);
  cache_stepptr_buf_reset(&caller->spp_reverse);

  uint32_t stepid = snode->first_step;
  GCacheStep *step;

  while(stepid != UINT32_MAX)
  {
    step = graph_cache_step(&caller->cache, stepid);
    if(step->orient == FORWARD) {
      cache_stepptr_buf_add(&caller->spp_forward, step);
    } else {
      cache_stepptr_buf_add(&caller->spp_reverse, step);
    }
    stepid = step->next_step;
  }

  // Filter out non-bubbles
  remove_non_bubbles(caller, &caller->spp_forward);
  remove_non_bubbles(caller, &caller->spp_reverse);
}

static void write_bubbles_to_file(BubbleCaller *caller)
{
  // Loop over supernodes checking if they are 3p flanks
  size_t snode_count = graph_cache_num_snodes(&caller->cache);
  GCacheSnode *snode;
  size_t i;

  for(i = 0; i < snode_count; i++)
  {
    snode = graph_cache_snode(&caller->cache, i);
    find_bubbles_ending_with(caller, snode);

    if(caller->spp_forward.len > 1)
      print_bubble(caller, caller->spp_forward.data, caller->spp_forward.len);
    if(caller->spp_reverse.len > 1)
      print_bubble(caller, caller->spp_reverse.data, caller->spp_reverse.len);
  }
}

void bubble_caller_node(hkey_t hkey, BubbleCaller *caller)
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
}

void bubble_caller(void *args)
{
  BubbleCaller *caller = (BubbleCaller*)args;

  HASH_ITERATE_PART(&caller->db_graph->ht, caller->threadid, caller->nthreads,
                    bubble_caller_node, caller);
}

void invoke_bubble_caller(size_t num_of_threads, BubbleCallingPrefs prefs,
                          gzFile gzout, const char *out_path,
                          const dBGraph *db_graph)
{
  ctx_assert(db_graph->num_edge_cols == 1);
  ctx_assert(db_graph->node_in_cols != NULL);

  status("Calling bubbles with %zu threads, output: %s", num_of_threads, out_path);

  // Print header
  bubble_caller_print_header(gzout, out_path, db_graph);

  BubbleCaller *callers = bubble_callers_new(num_of_threads, prefs,
                                             gzout, db_graph);

  // Run
  util_run_threads(callers, num_of_threads, sizeof(callers[0]),
                   num_of_threads, bubble_caller);

  // Report number of bubble called+printed
  size_t num_of_bubbles = callers[0].num_bubbles_ptr[0];
  char num_bubbles_str[100];
  ulong_to_str(num_of_bubbles, num_bubbles_str);
  status("%s bubbles called with Paths-Bubble-Caller\n", num_bubbles_str);

  // Clean up
  bubble_callers_destroy(callers, num_of_threads);
}
