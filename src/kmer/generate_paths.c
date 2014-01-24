#include "global.h"

#include <pthread.h>
#include "msgpool.h" // pool for getting jobs
#include <unistd.h> // usleep

#include "generate_paths.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "db_alignment.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "path_store.h"
#include "path_format.h"
#include "path_store_thread_safe.h"
#include "async_read_io.h"

//
// Multithreaded code to add paths to the graph from sequence data
// Uses async_read_io to allow multiple readers, multiple workers
//


// #define CTXVERBOSE 1

struct GenPathWorker
{
  pthread_t thread;

  // We take jobs from the pool
  MsgPool *pool;
  AsyncIOData data; // current data
  GeneratePathsTask task; // current task

  dBGraph *const db_graph;
  GraphWalker wlk, wlk2;
  RepeatWalker rptwlk, rptwlk2;

  // Nodes from the reads
  dBAlignment alignment; // raw nodes from reads with gaps
  // contig with gaps filled
  // we use revcontig when walking backwards
  dBNodeBuffer contig, revcontig;

  // Statistics on gap traversal
  SeqLoadingStats *stats;
  uint64_t *gap_ins_histgrm, *gap_err_histgrm, histgrm_len;

  // Nucleotides and positions of junctions
  Nucleotide *nuc_fw, *nuc_rv;
  size_t *pos_fw, *pos_rv;
  size_t num_fw, num_rv, junc_arrsize;

  // Packed representation of path
  uint8_t *packed;
  size_t packed_memcap;
};

#define INIT_BUFLEN 1024

// seq gap of N bases can be filled by MAX2(0, N±(N*GAP_VARIANCE+GAP_WIGGLE))
#define GAP_VARIANCE 0.1
#define GAP_WIGGLE 5

#define MAX_CONTEXT 1000

// Should we print all paths?
boolean gen_paths_print_inserts = false;
volatile size_t print_contig_id = 0;


size_t gen_paths_worker_est_mem(const dBGraph *db_graph)
{
  size_t job_mem = 1024*4; // Assume 1024 bytes per read, 2 reads, seq+qual
  size_t graph_walker_mem = 2*graph_walker_est_mem();
  size_t rpt_wlker_mem = 2*rpt_walker_est_mem(db_graph->ht.capacity, 22);
  size_t alignment_mem = db_alignment_est_mem();
  size_t contig_mem = 2*INIT_BUFLEN*sizeof(dBNode);
  size_t gap_hist_mem = 2*INIT_BUFLEN*sizeof(size_t);
  size_t junc_mem = INIT_BUFLEN * sizeof(Nucleotide) +
                    INIT_BUFLEN * sizeof(Nucleotide) +
                    INIT_BUFLEN * sizeof(size_t) +
                    INIT_BUFLEN * sizeof(size_t);
  size_t packed_mem = INIT_BUFLEN;

  return job_mem + graph_walker_mem + rpt_wlker_mem + alignment_mem + contig_mem +
         gap_hist_mem + junc_mem + packed_mem + sizeof(GenPathWorker);
}

static void _gen_paths_worker_alloc(GenPathWorker *wrkr, dBGraph *db_graph)
{
  GenPathWorker tmp = {.db_graph = db_graph, .pool = NULL};

  // Current job
  asynciodata_alloc(&tmp.data);

  // Graph traversal
  graph_walker_alloc(&tmp.wlk);
  graph_walker_alloc(&tmp.wlk2);
  rpt_walker_alloc(&tmp.rptwlk, db_graph->ht.capacity, 22); // 4MB
  rpt_walker_alloc(&tmp.rptwlk2, db_graph->ht.capacity, 22); // 4MB

  // Node buffers
  db_alignment_alloc(&tmp.alignment);
  db_node_buf_alloc(&tmp.contig, INIT_BUFLEN);
  db_node_buf_alloc(&tmp.revcontig, INIT_BUFLEN);
  // Gap Histogram
  tmp.histgrm_len = INIT_BUFLEN;
  tmp.gap_ins_histgrm = calloc2(tmp.histgrm_len, sizeof(*tmp.gap_ins_histgrm));
  tmp.gap_err_histgrm = calloc2(tmp.histgrm_len, sizeof(*tmp.gap_err_histgrm));
  // Stats
  tmp.stats = seq_loading_stats_create(0);
  // Junction data
  tmp.junc_arrsize = INIT_BUFLEN;
  tmp.nuc_fw = malloc2(tmp.junc_arrsize * sizeof(*tmp.nuc_fw));
  tmp.nuc_rv = malloc2(tmp.junc_arrsize * sizeof(*tmp.nuc_rv));
  tmp.pos_fw = malloc2(tmp.junc_arrsize * sizeof(*tmp.pos_fw));
  tmp.pos_rv = malloc2(tmp.junc_arrsize * sizeof(*tmp.pos_rv));
  tmp.num_fw = tmp.num_rv = 0;
  // Packed path temporary space
  tmp.packed_memcap = INIT_BUFLEN;
  tmp.packed = malloc2(tmp.packed_memcap);

  memcpy(wrkr, &tmp, sizeof(GenPathWorker));
}

static void _gen_paths_worker_dealloc(GenPathWorker *wrkr)
{
  asynciodata_dealloc(&wrkr->data);
  graph_walker_dealloc(&wrkr->wlk);
  graph_walker_dealloc(&wrkr->wlk2);
  rpt_walker_dealloc(&wrkr->rptwlk);
  rpt_walker_dealloc(&wrkr->rptwlk2);
  db_alignment_dealloc(&wrkr->alignment);
  db_node_buf_dealloc(&wrkr->contig);
  db_node_buf_dealloc(&wrkr->revcontig);
  free(wrkr->gap_ins_histgrm);
  free(wrkr->gap_err_histgrm);
  seq_loading_stats_free(wrkr->stats);
  free(wrkr->nuc_fw); free(wrkr->nuc_rv);
  free(wrkr->pos_fw); free(wrkr->pos_rv);
  free(wrkr->packed);
}


GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph)
{
  size_t i;
  GenPathWorker *workers = malloc2(n * sizeof(GenPathWorker));
  for(i = 0; i < n; i++) _gen_paths_worker_alloc(&workers[i], graph);
  return workers;
}

void gen_paths_workers_dealloc(GenPathWorker *workers, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) _gen_paths_worker_dealloc(&workers[i]);
  free(workers);
}

static inline void worker_gap_cap(GenPathWorker *wrkr, size_t max_gap)
{
  if(wrkr->histgrm_len < max_gap) {
    max_gap = roundup2pow(max_gap);
    wrkr->gap_ins_histgrm = realloc2(wrkr->gap_ins_histgrm, max_gap*sizeof(uint64_t));
    wrkr->gap_err_histgrm = realloc2(wrkr->gap_err_histgrm, max_gap*sizeof(uint64_t));
    // Zero new memory
    size_t newmem = (max_gap-wrkr->histgrm_len)*sizeof(uint64_t);
    memset(wrkr->gap_ins_histgrm+wrkr->histgrm_len, 0, newmem);
    memset(wrkr->gap_err_histgrm+wrkr->histgrm_len, 0, newmem);
    wrkr->histgrm_len = max_gap;
  }
}

static inline void worker_nuc_cap(GenPathWorker *wrkr, size_t cap)
{
  if(cap > wrkr->junc_arrsize) {
    wrkr->junc_arrsize = roundup2pow(cap);
    wrkr->nuc_fw = realloc2(wrkr->nuc_fw, wrkr->junc_arrsize * sizeof(Nucleotide));
    wrkr->nuc_rv = realloc2(wrkr->nuc_rv, wrkr->junc_arrsize * sizeof(Nucleotide));
    wrkr->pos_fw = realloc2(wrkr->pos_fw, wrkr->junc_arrsize * sizeof(size_t));
    wrkr->pos_rv = realloc2(wrkr->pos_rv, wrkr->junc_arrsize * sizeof(size_t));
  }
}

static inline void worker_packed_cap(GenPathWorker *wrkr, size_t nbytes)
{
  if(nbytes > wrkr->packed_memcap) {
    wrkr->packed_memcap = roundup2pow(nbytes);
    wrkr->packed = realloc2(wrkr->packed, wrkr->packed_memcap * sizeof(Nucleotide));
  }
}

static inline void merge_arrays(uint64_t *restrict a,
                                const uint64_t *restrict b,
                                size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) a[i] += b[i];
}

// Merge stats into workers[0]
void generate_paths_merge_stats(GenPathWorker *wrkrs, size_t num_workers)
{
  size_t i, n;

  if(num_workers <= 1) return;

  for(i = 1; i < num_workers; i++)
  {
    n = wrkrs[i].histgrm_len;
    worker_gap_cap(&wrkrs[0], n);
    merge_arrays(wrkrs[0].gap_ins_histgrm, wrkrs[i].gap_ins_histgrm, n);
    merge_arrays(wrkrs[0].gap_err_histgrm, wrkrs[i].gap_err_histgrm, n);

    // Zero arrays
    size_t arr_bytes = sizeof(*wrkrs[i].gap_ins_histgrm)*wrkrs[i].histgrm_len;
    memset(wrkrs[i].gap_ins_histgrm, 0, arr_bytes);
    memset(wrkrs[i].gap_err_histgrm, 0, arr_bytes);
  }
}


// assume nbits > 0
#define bits_in_top_word(nbits) ((((nbits) - 1) & 3) + 1)

// Returns number of paths added
static inline size_t _juncs_to_paths(const size_t *restrict pos_pl,
                                     const size_t *restrict pos_mn,
                                     size_t num_pl, size_t num_mn,
                                     const Nucleotide *nuc_pl,
                                     const boolean pl_is_fw,
                                     GenPathWorker *wrkr)
{
  size_t i, num_added = 0;
  size_t ctpcol = wrkr->task.ctpcol;

  hkey_t node;
  Orientation orient;
  size_t start_mn, start_pl, pos;
  PathLen plen, plen_orient;
  boolean added;

  const dBNode *nodes = wrkr->contig.data;
  dBGraph *db_graph = wrkr->db_graph;

  // Create packed path with four diff offsets (0..3), point to correct one
  size_t pckd_memsize = (num_pl+3)/4;
  uint8_t *packed_ptrs[4], *packed_ptr;

  worker_packed_cap(wrkr, pckd_memsize*4);
  for(i = 0; i < 4; i++) packed_ptrs[i] = wrkr->packed + i*pckd_memsize;

  pack_bases(packed_ptrs[0]+sizeof(PathLen), nuc_pl, num_pl);
  packed_cpy(packed_ptrs[1]+sizeof(PathLen), packed_ptrs[0], 1, num_pl);
  packed_cpy(packed_ptrs[2]+sizeof(PathLen), packed_ptrs[0], 2, num_pl);
  packed_cpy(packed_ptrs[3]+sizeof(PathLen), packed_ptrs[0], 3, num_pl);

  // pl => plus in direction
  // mn => minus against direction

  // 0,1,2,3
  // 3,2,1

  // DEV: should always add longest to shortest path
  // this is not the case with rev currently?
  for(start_pl = 0, start_mn = num_mn-1; start_mn != SIZE_MAX; start_mn--)
  {
    if(pl_is_fw)
      while(start_pl < num_pl && pos_pl[start_pl] < pos_mn[start_mn]) start_pl++;
    else
      while(start_pl < num_pl && pos_pl[start_pl] > pos_mn[start_mn]) start_pl++;

    if(start_pl == num_pl) break;

    // nodes[pos] is the kmer we are adding this path to
    pos = (pl_is_fw ? pos_mn[start_mn] - 1 : pos_mn[start_mn] + 1);
    start_pl -= (start_pl > 0 && pos_pl[start_pl-1] == pos);

    // Check path is not too long (MAX_PATHLEN is the limit)
    plen = (PathLen)MIN2(num_pl - start_pl, MAX_PATHLEN);
    node = nodes[pos].key;
    orient = pl_is_fw ? nodes[pos].orient : rev_orient(nodes[pos].orient);

    dBNode tmp = {.key = node, .orient = orient};
    path_format_is_path_valid(db_graph, tmp, wrkr->task.ctxcol,
                              nuc_pl+start_pl, plen);

    #ifdef CTXVERBOSE
      char kmerstr[MAX_KMER_SIZE+1];
      binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, kmerstr);
      printf(" %s:%i) start_pl: %zu start_mn: %zu {%zu}\n", kmerstr, orient,
             start_pl, start_mn, pos_mn[start_mn]);
    #endif

    // Write orient and length to packed representation
    plen_orient = plen | ((PathLen)orient << PATH_LEN_BITS);
    packed_ptr = packed_ptrs[start_pl&3] + start_pl/4;
    memcpy(packed_ptr, &plen_orient, sizeof(PathLen));

    // mask top byte!
    size_t top_idx = sizeof(PathLen)+(plen-1)/4;
    uint8_t top_byte = packed_ptr[top_idx];
    packed_ptr[top_idx] &= (uint8_t)(255 >> (8 - bits_in_top_word(plen)));

    added = path_store_mt_find_or_add(node, db_graph, ctpcol, packed_ptr, plen);
    packed_ptr[top_idx] = top_byte; // restore top byte

    #ifdef CTXVERBOSE
      printf("We %s\n", added ? "added" : "abandoned");
    #endif

    // If the path already exists, all of its subpaths also already exist
    if(!added && plen < MAX_PATHLEN) break;
    num_added++;
  }

  return num_added;
}

#undef bits_in_top_word

static void worker_junctions_to_paths(GenPathWorker *wrkr)
{
  size_t num_fw = wrkr->num_fw, num_rv = wrkr->num_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide *nuc_fw = wrkr->nuc_fw, *nuc_rv = wrkr->nuc_rv;

  worker_packed_cap(wrkr, (MAX2(num_fw, num_rv)+3)/4);

  // status("num_fw: %zu num_rv: %zu", num_fw, num_rv);

  // Reverse pos_rv, nuc_rv
  size_t i, j, tmp_pos;
  Nucleotide tmp_nuc;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp_pos);
    SWAP(nuc_rv[i], nuc_rv[j], tmp_nuc);
  }

  // _juncs_to_paths returns the number of paths added
  size_t n;

  n = _juncs_to_paths(pos_fw, pos_rv, num_fw, num_rv, nuc_fw, true, wrkr);
  if(n) _juncs_to_paths(pos_rv, pos_fw, num_rv, num_fw, nuc_rv, false, wrkr);
}

static void worker_contig_to_junctions(GenPathWorker *wrkr)
{
  // status("nodebuf: %zu", wrkr->contig.len+MAX_KMER_SIZE+1);
  worker_nuc_cap(wrkr, wrkr->contig.len);

  if(gen_paths_print_inserts) {
    pthread_mutex_lock(&biglock);
    printf(">contig%zu\n", print_contig_id++);
    db_nodes_print(wrkr->contig.data, wrkr->contig.len, wrkr->db_graph, stdout);
    printf("\n");
    pthread_mutex_unlock(&biglock);
  }

  // Find forks in this colour
  Edges edges;
  size_t  i, num_fw = 0, num_rv = 0;
  int indegree, outdegree;
  Nucleotide *nuc_fw = wrkr->nuc_fw, *nuc_rv = wrkr->nuc_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide nuc;
  BinaryKmer bkmer;

  dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = wrkr->db_graph->kmer_size;

  const dBNode *nodes = wrkr->contig.data;
  const size_t contig_len = wrkr->contig.len;
  const size_t ctxcol = wrkr->task.ctxcol;

  for(i = 0; i < contig_len; i++)
  {
    edges = db_node_edges(db_graph, ctxcol, nodes[i].key);
    outdegree = edges_get_outdegree(edges, nodes[i].orient);
    indegree = edges_get_indegree(edges, nodes[i].orient);

    if(indegree > 1 && i > 0)
    {
      bkmer = db_node_bkmer(db_graph, nodes[i-1].key);
      nuc = nodes[i-1].orient == FORWARD
              ? dna_nuc_complement(binary_kmer_first_nuc(bkmer, kmer_size))
              : binary_kmer_last_nuc(bkmer);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
    }

    if(outdegree > 1 && i+1 < contig_len && num_rv > 0)
    {
      bkmer = db_node_bkmer(db_graph, nodes[i+1].key);
      nuc = db_node_last_nuc(bkmer, nodes[i+1].orient, kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
    }
  }

  wrkr->num_fw = num_fw;
  wrkr->num_rv = num_rv;

  if(num_fw > 0 && num_rv > 0)
    worker_junctions_to_paths(wrkr);
}


static void prime_for_traversal(GraphWalker *wlk,
                                const dBNode *block, size_t n, boolean forward,
                                size_t ctxcol, size_t ctpcol,
                                const dBGraph *db_graph)
{
  assert(n > 0);
  dBNode node0;

  if(n > MAX_CONTEXT) {
    if(forward) block = block + n - MAX_CONTEXT;
    n = MAX_CONTEXT;
  }

  if(forward) { node0 = block[0]; block++; }
  else { node0 = db_node_reverse(block[n-1]); }

  // status("prime_for_traversal() %zu:%i n=%zu %s", (size_t)node0.key, node0.orient,
  //        n, forward ? "forward" : "reverse");

  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node0.key, node0.orient);
  // graph_walker_fast_traverse(wlk, block, n-1, forward);
  graph_walker_slow_traverse(wlk, block, n-1, forward);

  // char tmpbkmer[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(wlk->bkmer, db_graph->kmer_size, tmpbkmer);
  // status("  primed: %s:%i\n", tmpbkmer, wlk->orient);
}

// Resets node buffer, GraphWalker and RepeatWalker after use
static boolean worker_traverse_gap(dBNode end_node, dBNodeBuffer *contig,
                                   size_t gap_min, size_t gap_max,
                                   size_t *gap_len_ptr,
                                   GraphWalker *wlk, RepeatWalker *rptwlk)
{
  boolean traversed = false;

  size_t init_len = contig->len, max_len = contig->len + gap_max;
  db_node_buf_ensure_capacity(contig, contig->len + gap_max + 1);

  while(contig->len < max_len && graph_traverse(wlk) &&
        rpt_walker_attempt_traverse(rptwlk, wlk, wlk->node, wlk->orient, wlk->bkmer))
  {
    if(wlk->node == end_node.key && wlk->orient == end_node.orient) {
      traversed = true;
      break;
    }
    dBNode node = {.key = wlk->node, .orient = wlk->orient};
    contig->data[contig->len++] = node;
  }

  size_t gap_len = contig->len - init_len;

  // Clean up GraphWalker, RepeatWalker
  graph_walker_finish(wlk);
  rpt_walker_fast_clear(rptwlk, contig->data+init_len, gap_len);

  // Fail if we bridged but too short
  if(gap_len < gap_min) traversed = false;
  if(!traversed) contig->len = init_len;

  *gap_len_ptr = gap_len;
  return traversed;
}

// Traversed from both sides of the gap
// nodes in contig1 will be the *reverse*
// Resets node buffers, GraphWalkers and RepeatWalkers after use
static boolean worker_traverse_gap2(dBNodeBuffer *contig0, dBNodeBuffer *contig1,
                                    size_t gap_min, size_t gap_max,
                                    size_t *gap_len_ptr,
                                    GraphWalker *wlk0, RepeatWalker *rptwlk0,
                                    GraphWalker *wlk1, RepeatWalker *rptwlk1)
{
  const size_t init_len0 = contig0->len, init_len1 = contig1->len;

  db_node_buf_ensure_capacity(contig0, contig0->len + gap_max + 1);
  db_node_buf_ensure_capacity(contig1, contig1->len + gap_max + 1);

  boolean traversed = false;

  boolean use[2] = {true,true};
  GraphWalker *wlk[2] = {wlk0, wlk1};
  RepeatWalker *rptwlk[2] = {rptwlk0, rptwlk1};
  dBNodeBuffer *contig[2] = {contig0, contig1};
  size_t i, gap_len = 0;

  //
  // char str[MAX_KMER_SIZE+1];
  // BinaryKmer bkmer;
  // for(i = 0; i < 2; i++) {
  //   bkmer = db_node_bkmer(wlk0->db_graph, wlk[i]->node);
  //   binary_kmer_to_str(bkmer, wlk0->db_graph->kmer_size, str);
  //   status("primed %zu: %zu:%i %s\n", i, (size_t)wlk[i]->node, (int)wlk[i]->orient, str);
  // }
  //

  while(gap_len < gap_max && (use[0] || use[1])) {
    for(i = 0; i < 2; i++) {
      if(use[i]) {
        use[i] = (graph_traverse(wlk[i]) &&
                  rpt_walker_attempt_traverse(rptwlk[i], wlk[i],
                                              wlk[i]->node, wlk[i]->orient,
                                              wlk[i]->bkmer));

        if(use[i]) {
          //
          // bkmer = db_node_bkmer(wlk0->db_graph, wlk[i]->node);
          // binary_kmer_to_str(bkmer, wlk0->db_graph->kmer_size, str);
          // status("%zu: %zu:%i %s\n", i, (size_t)wlk[i]->node, (int)wlk[i]->orient, str);
          //

          if(wlk[0]->node == wlk[1]->node && wlk[0]->orient != wlk[1]->orient) {
            traversed = true;
            use[0] = use[1] = false; // set both to false to exit loop
            break;
          }

          dBNode node = {.key = wlk[i]->node, .orient = wlk[i]->orient};
          contig[i]->data[contig[i]->len++] = node;
          gap_len++;
        }
      }
    }
  }

  // status("worker_traverse_gap2() gap_len:%zu gap_min:%zu gap_max:%zu %s",
  //        gap_len, gap_min, gap_max, traversed ? "Good" : "Bad");

  // Clean up GraphWalker, RepeatWalker
  graph_walker_finish(wlk0);
  rpt_walker_fast_clear(rptwlk0, contig0->data+init_len0, contig0->len-init_len0);
  graph_walker_finish(wlk1);
  rpt_walker_fast_clear(rptwlk1, contig1->data+init_len1, contig1->len-init_len1);

  // Fail if we bridged but too short
  if(gap_len < gap_min) traversed = false;
  if(!traversed) { contig0->len = init_len0; contig1->len = init_len1; }

  *gap_len_ptr = gap_len;
  return traversed;
}

static boolean traverse_one_way(GenPathWorker *wrkr,
                                size_t start_idx, size_t gap_idx, size_t end_idx,
                                size_t gap_min, size_t gap_max,
                                size_t *gap_len_ptr)
{
  boolean traversed;
  const size_t block0len = gap_idx-start_idx, block1len = end_idx-gap_idx;
  const size_t ctxcol = wrkr->task.ctxcol, ctpcol = wrkr->task.ctpcol;
  const dBNodeBuffer *nodes = &wrkr->alignment.nodes;
  dBGraph *db_graph = wrkr->db_graph;

  // Start traversing forward
  prime_for_traversal(&wrkr->wlk, nodes->data+start_idx, block0len, true,
                      ctxcol, ctpcol, db_graph);

  traversed
    = worker_traverse_gap(nodes->data[gap_idx],
                          &wrkr->contig,
                          gap_min, gap_max, gap_len_ptr,
                          &wrkr->wlk, &wrkr->rptwlk);

  if(traversed) return true;

  // Start traversing backwards
  prime_for_traversal(&wrkr->wlk, nodes->data+gap_idx, block1len, false,
                      ctxcol, ctpcol, db_graph);

  traversed
    = worker_traverse_gap(db_node_reverse(nodes->data[gap_idx-1]),
                          &wrkr->revcontig,
                          gap_min, gap_max, gap_len_ptr,
                          &wrkr->wlk, &wrkr->rptwlk);

  return traversed;
}

static boolean traverse_two_way(GenPathWorker *wrkr,
                                size_t start_idx, size_t gap_idx, size_t end_idx,
                                size_t gap_min, size_t gap_max,
                                size_t *gap_len_ptr)
{
  const size_t block0len = gap_idx-start_idx, block1len = end_idx-gap_idx;
  const size_t ctxcol = wrkr->task.ctxcol, ctpcol = wrkr->task.ctpcol;
  const dBNodeBuffer *nodes = &wrkr->alignment.nodes;
  dBGraph *db_graph = wrkr->db_graph;

  prime_for_traversal(&wrkr->wlk, nodes->data+start_idx, block0len, true,
                      ctxcol, ctpcol, db_graph);
  prime_for_traversal(&wrkr->wlk2, nodes->data+gap_idx, block1len, false,
                      ctxcol, ctpcol, db_graph);

  // status("Traversing two-way... %zu[len:%zu]:%zu[len:%zu]",
  //        start_idx, block0len, gap_idx, block1len);

  return worker_traverse_gap2(&wrkr->contig, &wrkr->revcontig,
                              gap_min, gap_max, gap_len_ptr,
                              &wrkr->wlk, &wrkr->rptwlk,
                              &wrkr->wlk2, &wrkr->rptwlk2);
}

// Construct a contig from the gapped alignment
// Returns the next start index of the alignment
static void worker_generate_contigs(GenPathWorker *wrkr)
{
  const dBAlignment *algnmnt = &wrkr->alignment;
  const dBNodeBuffer *nodes = &algnmnt->nodes;
  const uint32Buffer *gaps = &algnmnt->gaps;

  const size_t num_align_nodes = nodes->len;

  dBNodeBuffer *contig = &wrkr->contig, *revcontig = &wrkr->revcontig;

  if(nodes->len == 0) return;

  // worker_generate_contigs ensures contig is at least nodes->len long
  boolean both_reads = (algnmnt->used_r1 && algnmnt->used_r2);
  size_t i, j, start_idx = 0, gap_idx, end_idx, block0len, block1len;
  size_t gap_est, gap_len, gap_min, gap_max;
  boolean traversed = false;

  // Assemble a contig in... contig
  db_node_buf_reset(contig);
  end_idx = db_alignment_next_gap(algnmnt, 0);
  block0len = end_idx;

  while(end_idx < num_align_nodes)
  {
    // block0 [start_idx..gap_idx-1], block1 [gap_idx..end_idx]
    gap_idx = end_idx;
    end_idx = db_alignment_next_gap(algnmnt, end_idx);
    block1len = end_idx - gap_idx;

    // status("start_idx: %zu gap_idx: %zu end_idx: %zu", start_idx, gap_idx, end_idx);

    // We've got a gap to traverse
    long gap_min_long, gap_max_long;
    boolean is_mp;

    // Get bound for acceptable bridge length (min,max length values)
    is_mp = (both_reads && gap_idx == algnmnt->r2strtidx);
    gap_est = gaps->data[gap_idx];
    if(is_mp) gap_est += algnmnt->r1enderr;

    gap_min_long = (long)gap_est - (long)(gap_est * GAP_VARIANCE + GAP_WIGGLE);
    gap_max_long = (long)gap_est + (long)(gap_est * GAP_VARIANCE + GAP_WIGGLE);

    if(is_mp) {
      gap_min_long += wrkr->task.ins_gap_min;
      gap_max_long += wrkr->task.ins_gap_max;
    }

    gap_min = (size_t)MAX2(0, gap_min_long);
    gap_max = (size_t)MAX2(0, gap_max_long);

    // Copy block0 (start_idx..end_idx) into contig buffer
    if(!traversed) {
      db_node_buf_reset(contig);
      db_node_buf_append(contig, nodes->data+start_idx, block0len);
    }

    db_node_buf_reset(revcontig);

    // Alternative traversing from both sides
    if(wrkr->task.one_way_gap_traverse)
      traversed = traverse_one_way(wrkr, start_idx, gap_idx, end_idx,
                                   gap_min, gap_max, &gap_len);
    else
      traversed = traverse_two_way(wrkr, start_idx, gap_idx, end_idx,
                                   gap_min, gap_max, &gap_len);

    // status("traversal: %s!\n", traversed ? "worked" : "failed");

    if(traversed)
    {
      // reverse and copy from revcontig -> contig
      db_node_buf_ensure_capacity(contig, contig->len + revcontig->len);
      for(i = 0, j = contig->len+revcontig->len-1; i < revcontig->len; i++, j--)
        contig->data[j] = db_node_reverse(revcontig->data[i]);

      contig->len += revcontig->len;

      // Copy block1
      db_node_buf_append(contig, nodes->data+gap_idx, block1len);

      // Update gap stats
      worker_gap_cap(wrkr, gap_len);

      // gap_est is the sequence gap (number of missing kmers)
      if(is_mp) {
        size_t ins_gap = (size_t)MAX2((long)gap_len - (long)gap_est, 0);
        wrkr->gap_ins_histgrm[ins_gap]++;
        // size_t err_gap = gap_len - ins_gap;
        // if(err_gap) wrkr->gap_err_histgrm[err_gap]++;
      } else {
        wrkr->gap_err_histgrm[gap_len]++;
      }
    }
    else
    {
      // Process the contig we have constructed
      worker_contig_to_junctions(wrkr);

      // Move start up the the gap we are stuck on
      start_idx = gap_idx;
      block0len = block1len;
    }
  }

  if(!traversed) {
    db_node_buf_reset(contig);
    db_node_buf_append(contig, nodes->data+start_idx, block0len);
  }

  // Use last contig constructed
  worker_contig_to_junctions(wrkr);
}

static void worker_reads_to_nodebuf(GenPathWorker *wrkr)
{
  AsyncIOData *data = &wrkr->data;
  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  // if(r1->seq.end < 100) {
  //   printf("1>%s\n", r1->seq.b);
  //   if(r2 != NULL) printf("2>%s\n", r2->seq.b);
  // }

  uint8_t fq_cutoff1, fq_cutoff2;
  fq_cutoff1 = fq_cutoff2 = wrkr->task.fq_cutoff;

  if(fq_cutoff1 > 0)
  {
    fq_cutoff1 += wrkr->data.fq_offset1;
    fq_cutoff2 += wrkr->data.fq_offset2;
  }

  uint8_t hp_cutoff = wrkr->task.hp_cutoff;

  // Second read is in reverse orientation - need in forward
  if(r2 != NULL && wrkr->task.read_pair_FR)
    seq_read_reverse_complement(r2);

  // Update stats
  if(r2 == NULL) wrkr->stats->num_se_reads++;
  else wrkr->stats->num_pe_reads += 2;

  db_alignment_from_reads(&wrkr->alignment, r1, r2,
                          fq_cutoff1, fq_cutoff2, hp_cutoff,
                          wrkr->db_graph);

  // For debugging
  // db_alignment_print(&wrkr->alignment, wrkr->db_graph);
}

static inline void worker_do_job(GenPathWorker *wrkr)
{
  worker_reads_to_nodebuf(wrkr);
  worker_generate_contigs(wrkr);
}

// pthread method, loop: grabs job, does processing
static void* generate_paths_worker(void *ptr)
{
  GenPathWorker *wrkr = (GenPathWorker*)ptr;

  AsyncIOData data;
  while(msgpool_read(wrkr->pool, &data, &wrkr->data)) {
    wrkr->data = data;
    memcpy(&wrkr->task, data.ptr, sizeof(GeneratePathsTask));
    worker_do_job(wrkr);
  }

  pthread_exit(NULL);
}

static void start_gen_path_workers(GenPathWorker *workers, size_t nworkers)
{
  size_t i;
  int rc;

  // Thread attribute joinable
  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  // Start threads
  for(i = 0; i < nworkers; i++) {
    rc = pthread_create(&workers[i].thread, &thread_attr,
                        generate_paths_worker, (void*)&workers[i]);
    if(rc != 0) die("Creating thread failed");
  }

  // Finished with thread attribute
  pthread_attr_destroy(&thread_attr);

  // Join threads
  for(i = 0; i < nworkers; i++) {
    rc = pthread_join(workers[i].thread, NULL);
    if(rc != 0) die("Joining thread failed");
  }
}

// workers array must be at least as long as tasks
// Stats are merged into workers[0]
void generate_paths(GeneratePathsTask *tasks, size_t num_tasks,
                    GenPathWorker *workers, size_t num_workers)
{
  status("[MkPaths] %zu input%s being handled by %zu worker%s",
         num_tasks, num_tasks != 1 ? "s" : "",
         num_workers, num_workers != 1 ? "s" : "");

  if(!num_tasks) return;
  assert(num_workers > 0);
  assert(workers[0].db_graph->path_kmer_locks != NULL);

  size_t i, max_gap = 0;
  MsgPool pool;
  // msgpool_alloc_yield(&pool, MSGPOOLRSIZE, sizeof(AsyncIOData));
  msgpool_alloc_spinlock(&pool, MSGPOOLRSIZE, sizeof(AsyncIOData));

  for(i = 0; i < num_tasks; i++) max_gap = MAX2(max_gap, tasks[i].ins_gap_max);

  for(i = 0; i < num_workers; i++) {
    workers[i].pool = &pool;
    worker_gap_cap(&workers[i], max_gap*2+100);
  }

  // Start async io reading
  AsyncIOWorker *asyncio_workers;
  AsyncIOReadTask *asyncio_tasks = malloc(num_tasks * sizeof(AsyncIOReadTask));

  for(i = 0; i < num_tasks; i++) {
    AsyncIOReadTask aio_task = {.file1 = tasks[i].file1,
                                .file2 = tasks[i].file2,
                                .fq_offset = tasks[i].fq_offset,
                                .ptr = &tasks[i]};

    memcpy(&asyncio_tasks[i], &aio_task, sizeof(AsyncIOReadTask));
  }

  asyncio_workers = asyncio_read_start(&pool, asyncio_tasks, num_tasks);

  // Run the workers until the pool is closed
  start_gen_path_workers(workers, num_workers);

  // status("Pool open: %i", pool.open);

  // Finish with the async io (waits until queue is empty)
  asyncio_read_finish(asyncio_workers, num_tasks);
  free(asyncio_tasks);
  msgpool_dealloc(&pool);

  // Merge gap counts into worker[0]
  generate_paths_merge_stats(workers, num_workers);
}


// Save gap size distribution
// base_fmt is the beginning of the file name - the reset is <num>.csv or something
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void gen_paths_dump_gap_sizes(const char *base_fmt,
                              const uint64_t *arr, size_t arrlen,
                              size_t kmer_size, boolean insert_sizes,
                              size_t nreads)
{
  assert(arrlen > 0);

  // Print summary statistics: min, mean, median, mode, max
  size_t i, min, max, total, ngaps = 0, mode = 0;
  max = total = arr[0];

  for(min = 0; min < arrlen && arr[min] == 0; min++) {}

  if(min == arrlen) {
    if(insert_sizes) status("No insert gaps traversed");
    else status("No seq error gaps traversed");
    return;
  }

  for(i = 1; i < arrlen; i++) {
    if(arr[i] > 0) max = i;
    if(arr[i] > arr[mode]) mode = i;
    ngaps += arr[i];
    total += arr[i] * i;
  }

  double mean = (double)total / ngaps;
  float median = find_hist_median(arr, arrlen, ngaps);

  size_t ninputs = insert_sizes ? nreads/2 : nreads;
  char ngaps_str[100], ninputs_str[100];
  ulong_to_str(ngaps, ngaps_str);
  ulong_to_str(ninputs, ninputs_str);

  status("%s size distribution: "
         "min: %zu mean: %.1f median: %.1f mode: %zu max: %zu",
         insert_sizes ? "Insert" : "Seq error gap",
         min, mean, median, mode, max);

  status("  Gaps per read%s: %s / %s [%.3f]",
         insert_sizes ? " pair" : "", ngaps_str, ninputs_str,
         (double)ngaps / ninputs);

  StrBuf *csv_dump = strbuf_new();
  FILE *fout;

  if(!futil_generate_filename(base_fmt, csv_dump)) {
    warn("Cannot dump gapsize");
    strbuf_free(csv_dump);
    return;
  }

  if((fout = fopen(csv_dump->buff, "w")) == NULL) {
    warn("Cannot dump gapsize [cannot open: %s]", csv_dump->buff);
    strbuf_free(csv_dump);
    return;
  }

  fprintf(fout, "gap_in_kmers\tbp\tcount\n");

  if(arrlen > 0)
  {
    size_t start = 0, end = arrlen-1;

    while(start < arrlen && arr[start] == 0) start++;
    while(end > start && arr[end] == 0) end--;

    for(i = start; i <= end; i++) {
      fprintf(fout, "%4zu\t%4li\t%4zu\n",
              i, (long)i-(long)kmer_size, (size_t)arr[i]);
    }
  }

  status("Contig %s sizes dumped to %s\n",
         insert_sizes ? "insert" : "gap", csv_dump->buff);

  fclose(fout);
  strbuf_free(csv_dump);
}

// Get histogram array
const uint64_t* gen_paths_get_ins_gap(GenPathWorker *worker, size_t *len)
{
  *len = (size_t)worker->histgrm_len;
  return worker->gap_ins_histgrm;
}

const uint64_t* gen_paths_get_err_gap(GenPathWorker *worker, size_t *len)
{
  *len = (size_t)worker->histgrm_len;
  return worker->gap_err_histgrm;
}

void gen_paths_get_stats(GenPathWorker *worker, size_t num_workers,
                         SeqLoadingStats *stats)
{
  size_t i;
  for(i = 0; i < num_workers; i++)
    seq_loading_stats_sum(stats, worker[i].stats);
}
