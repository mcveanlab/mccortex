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
#include "path_store_thread_safe.h"
#include "async_read_io.h"

//
// Multithreaded code to add paths to the graph from sequence data
// Uses async_read_io to allow multiple readers, multiple workers
//


// #define CTXVERBOSE 1

// DEV: move gap_ins_histgrm, gap_err_histgrm to GeneratePathsTask
// DEV: support for read orientations other than FR

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
#define BRIDGE_GAP_TWO_WAY 1

// Should we print all paths?
boolean gen_paths_print_inserts = false;


static inline void merge_arrays(uint64_t *restrict a,
                                const uint64_t *restrict b,
                                size_t n1, size_t n2)
{
  size_t i, n = MIN2(n1, n2);
  for(i = 0; i < n; i++) a[i] += b[i];
  // In case n2 is greater than n2
  for(i = n1; i < n2; i++) a[n1-1] += b[i];
}

// Merge stats into workers[0]
void generate_paths_merge_stats(GenPathWorker *wrkrs, size_t nwrkrs)
{
  size_t i, n1, n2;

  if(nwrkrs <= 1) return;

  for(i = 1; i < nwrkrs; i++)
  {
    n1 = wrkrs[0].histgrm_len; n2 = wrkrs[i].histgrm_len;
    merge_arrays(wrkrs[0].gap_ins_histgrm, wrkrs[i].gap_ins_histgrm, n1, n2);
    merge_arrays(wrkrs[0].gap_err_histgrm, wrkrs[i].gap_err_histgrm, n1, n2);
  }
}


size_t gen_paths_worker_est_mem(const dBGraph *db_graph)
{
  size_t job_mem = 1024*4; // Assume 1024 bytes per read, 2 reads, seq+qual
  size_t graph_walker_mem = graph_walker_est_mem();
  size_t rpt_wlker_mem = rpt_walker_est_mem(db_graph->ht.capacity, 22);
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
  tmp.gap_ins_histgrm = calloc2(tmp.histgrm_len*2, sizeof(size_t));
  tmp.gap_err_histgrm = tmp.gap_ins_histgrm + tmp.histgrm_len;
  // Junction data
  tmp.junc_arrsize = INIT_BUFLEN;
  tmp.nuc_fw = malloc2(tmp.junc_arrsize * sizeof(Nucleotide));
  tmp.nuc_rv = malloc2(tmp.junc_arrsize * sizeof(Nucleotide));
  tmp.pos_fw = malloc2(tmp.junc_arrsize * sizeof(size_t));
  tmp.pos_rv = malloc2(tmp.junc_arrsize * sizeof(size_t));
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

// assume nbits > 0
#define bits_in_top_word(nbits) ((((nbits) - 1) & 3) + 1)

// Returns number of paths added
static inline size_t _juncs_to_paths(const size_t *pos_pl, const size_t *pos_mn,
                                     size_t num_pl, size_t num_mn,
                                     const Nucleotide *nuc_pl, boolean fw,
                                     GenPathWorker *wrkr)
{
  size_t num_added = 0;
  size_t ctpcol = wrkr->task.ctpcol;

  hkey_t node;
  Orientation orient;
  size_t start_mn, start_pl, prev_start_pl = 0, pos;
  PathLen plen;
  boolean added;

  const dBNode *nodes = wrkr->contig.data;
  dBGraph *db_graph = wrkr->db_graph;

  // Packed path sequence
  PathLen *plen_ptr = (PathLen*)wrkr->packed;
  uint8_t *packed_seq = wrkr->packed+sizeof(PathLen);
  pack_bases(packed_seq, nuc_pl, num_pl);

  // pl => plus in direction
  // mn => minus against direction

  for(start_pl = 0, start_mn = num_mn-1; start_mn != SIZE_MAX; start_mn--)
  {
    while(start_pl < num_pl && pos_pl[start_pl] > pos_mn[start_mn]) start_pl++;
    if(start_pl == num_pl) break;

    pos = (fw ? pos_mn[start_mn] - 1 : pos_mn[start_mn] + 1);
    start_pl -= (start_pl > 0 && pos_pl[start_pl-1] == pos);

    // Check path is not too long (MAX_PATHLEN is the limit)
    size_t njuncs = MIN2(num_pl - start_pl, MAX_PATHLEN);

    plen = (PathLen)njuncs;
    node = nodes[pos].key;
    orient = rev_orient(nodes[pos].orient);

    #ifdef CTXVERBOSE
      binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, str);
      printf(" %s:%i) start_pl: %zu start_mn: %zu {%zu}\n", str, orient,
             start_pl, start_mn, pos_mn[start_mn]);
    #endif

    *plen_ptr = plen | (PathLen)(orient << PATH_LEN_BITS);
    right_shift_packed_bases(packed_seq, prev_start_pl - start_pl, plen);
    prev_start_pl = start_pl;

    // mask top byte!
    size_t top_idx = (plen-1)/4;
    uint8_t top_byte = packed_seq[top_idx];
    packed_seq[top_idx] &= (uint8_t)(255 >> (8 - bits_in_top_word(plen)));

    added = path_store_mt_find_or_add(node, db_graph, ctpcol,
                                      wrkr->packed, plen);

    packed_seq[top_idx] = top_byte;

    // If the path already exists, all of its subpaths also already exist
    if(!added) break;
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

  status("num_fw: %zu num_rv: %zu", num_fw, num_rv);

  // Reverse pos_rv, nuc_rv
  size_t i, j, tmp_pos;
  Nucleotide tmp_nuc;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp_pos);
    SWAP(nuc_rv[i], nuc_rv[j], tmp_nuc);
  }

  // _junc_to_paths returns the number of paths added
  size_t n;

  n = _juncs_to_paths(pos_fw, pos_rv, num_fw, num_rv, nuc_fw, true, wrkr);
  if(n) _juncs_to_paths(pos_rv, pos_fw, num_rv, num_fw, nuc_rv, false, wrkr);
}

static void worker_contig_to_junctions(GenPathWorker *wrkr)
{
  worker_nuc_cap(wrkr, wrkr->contig.len);

  // Find forks in this colour
  Edges edges;
  size_t  i, num_fw = 0, num_rv = 0, addfw, addrv;
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

    addfw = (outdegree > 1 && i+1 < contig_len);
    addrv = (indegree > 1 && i > 0);

    if(addfw) {
      bkmer = db_node_bkmer(db_graph, nodes[i+1].key);
      nuc = db_node_last_nuc(bkmer, nodes[i+1].orient, kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
    }
    if(addrv) {
      bkmer = db_node_bkmer(db_graph, nodes[i-1].key);
      nuc = nodes[i-1].orient == FORWARD
              ? dna_nuc_complement(binary_kmer_first_nuc(bkmer, kmer_size))
              : binary_kmer_last_nuc(bkmer);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
    }
  }

  wrkr->num_fw = num_fw;
  wrkr->num_rv = num_rv;

  if(num_fw > 0 && num_rv > 0)
    worker_junctions_to_paths(wrkr);
}


static void prime_for_traversal(GraphWalker *wlk,
                                const dBNode *block, size_t n, boolean forward,
                                size_t ctxcol, size_t ctpcol, dBGraph *db_graph)
{
  dBNode node0;

  if(n > MAX_CONTEXT) {
    if(forward) block = block + n - MAX_CONTEXT;
    n = MAX_CONTEXT;
  }

  if(forward) { node0 = block[0]; block++; }
  else { node0 = db_node_reverse(block[n-1]); }

  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node0.key, node0.orient);
  graph_walker_fast_traverse(wlk, block, n-1, forward, false);
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

  boolean use[2] = {true};
  GraphWalker *wlk[2] = {wlk0, wlk1};
  RepeatWalker *rptwlk[2] = {rptwlk0, rptwlk1};
  dBNodeBuffer *contig[2] = {contig0, contig1};
  size_t i, gap_len = 0;

  while(gap_len < gap_max && (use[0] || use[1]) && !traversed) {
    for(i = 0; i < 2; i++) {
      if(use[i]) {
        use[i] = (graph_traverse(wlk[i]) &&
                  rpt_walker_attempt_traverse(rptwlk[i], wlk[i],
                                              wlk[i]->node, wlk[i]->orient,
                                              wlk[i]->bkmer));

        if(use[i]) {
          dBNode node = {.key = wlk[i]->node, .orient = wlk[i]->orient};
          contig[i]->data[contig[i]->len++] = node;
          if(wlk[0]->node == wlk[1]->node && wlk[0]->orient != wlk[1]->orient) {
            traversed = true;
            break;
          }
          gap_len++;
        }
      }
    }
  }

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

  db_node_buf_reset(&wrkr->revcontig);

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

  db_node_buf_reset(&wrkr->revcontig);

  prime_for_traversal(&wrkr->wlk, nodes->data+start_idx, block0len, true,
                      ctxcol, ctpcol, db_graph);
  prime_for_traversal(&wrkr->wlk2, nodes->data+gap_idx, block1len, false,
                      ctxcol, ctpcol, db_graph);

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

  // worker_generate_contigs ensures contig is at least nodes->len long
  boolean both_reads = (algnmnt->used_r1 && algnmnt->used_r2);
  size_t i, j, start_idx = 0, gap_idx, end_idx, block0len, block1len;
    size_t gap_est, gap_len, gap_min, gap_max;
  boolean traversed = false;

  // Assemble a contig in... contig
  db_node_buf_reset(contig);
  end_idx = db_alignment_next_gap(algnmnt, 0);
  block1len = end_idx;

  while(end_idx < num_align_nodes)
  {
    // block0 [start_idx..gap_idx-1], block1 [gap_idx..end_idx]
    gap_idx = end_idx;
    end_idx = db_alignment_next_gap(algnmnt, start_idx);
    block0len = block1len;
    block1len = end_idx - gap_idx;

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

    // Alternative traversing from both sides
    if(!BRIDGE_GAP_TWO_WAY)
      traversed = traverse_one_way(wrkr, start_idx, gap_idx, end_idx,
                                   gap_min, gap_max, &gap_len);
    else
      traversed = traverse_two_way(wrkr, start_idx, gap_idx, end_idx,
                                   gap_min, gap_max, &gap_len);

    if(traversed)
    {
      // reverse and copy from revcontig -> contig
      db_node_buf_ensure_capacity(contig, contig->len + revcontig->len);
      for(i = 0, j = contig->len+revcontig->len-1; i < revcontig->len; i++, j--)
        contig->data[j] = db_node_reverse(revcontig->data[i]);

      // Copy block1
      db_node_buf_append(contig, nodes->data+gap_idx, block1len);

      // DEV: update gap stats
    }
    else
    {
      // Process the contig we have constructed
      worker_contig_to_junctions(wrkr);

      // Move start up the the gap we are stuck on
      start_idx = gap_idx;
    }
  }

  if(!traversed) {
    db_node_buf_reset(contig);
    db_node_buf_append(contig, nodes->data+start_idx, block0len);
  }

  // Use contig
  worker_contig_to_junctions(wrkr);
}

static void worker_read_to_nodebuf(GenPathWorker *worker)
{
  AsyncIOData *data = &worker->data;
  read_t *r1 = &data->r1, *r2 = data->r2.seq.end == 0 ? NULL : &data->r2;

  uint8_t fq_cutoff1, fq_cutoff2;
  fq_cutoff1 = fq_cutoff2 = worker->task.fq_cutoff;

  if(fq_cutoff1 > 0)
  {
    fq_cutoff1 += worker->data.fq_offset1;
    fq_cutoff2 += worker->data.fq_offset2;
  }

  uint8_t hp_cutoff = worker->task.hp_cutoff;

  // DEV: might not be FR
  if(r2 != NULL)
    seq_read_reverse_complement(r2);

  db_alignment_from_reads(&worker->alignment, r1, r2,
                          fq_cutoff1, fq_cutoff2, hp_cutoff,
                          worker->db_graph);
}

static inline void worker_do_job(GenPathWorker *worker)
{
  worker_read_to_nodebuf(worker);
  worker_generate_contigs(worker);
}

// pthread method, loop: grabs job, does processing
static void* generate_paths_worker(void *ptr)
{
  GenPathWorker *wrkr = (GenPathWorker*)ptr;

  AsyncIOData data;
  while(msgpool_read(wrkr->pool, &data, &wrkr->data)) {
    wrkr->data = data;
    wrkr->task = *(GeneratePathsTask*)data.ptr;
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
  if(!num_tasks) return;
  assert(num_workers > 0);
  assert(workers[0].db_graph->path_kmer_locks != NULL);

  size_t i;
  MsgPool pool;
  msgpool_alloc_yield(&pool, MSGPOOLRSIZE, sizeof(AsyncIOData));

  for(i = 0; i < num_workers; i++) workers[i].pool = &pool;

  // Start async io reading
  AsyncIOWorker *asyncio_workers;
  AsyncIOReadTask *asyncio_tasks = malloc(num_tasks * sizeof(AsyncIOReadTask));

  for(i = 0; i < num_tasks; i++) {
    AsyncIOReadTask aio_task = {.file1 = tasks[i].file1,
                                .file2 = tasks[i].file2,
                                .fq_offset = tasks[i].fq_offset,
                                .ptr = &tasks[i]};
    asyncio_tasks[i] = aio_task;
  }

  asyncio_workers = asyncio_read_start(&pool, asyncio_tasks, num_tasks);

  // Run the workers until the pool is closed
  start_gen_path_workers(workers, num_workers);

  // Finish with the async io (waits until queue is empty)
  asyncio_read_finish(asyncio_workers, num_tasks);
  free(asyncio_tasks);
  msgpool_dealloc(&pool);

  // Merge stats into workers[0]
  generate_paths_merge_stats(workers, num_workers);
}



// Save gap size distribution
// base_fmt is the beginning of the file name - the reset is <num>.csv or something
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void gen_paths_dump_gap_lens(const char *base_fmt,
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

  char ngaps_str[100], nreads_str[100];
  ulong_to_str(ngaps, ngaps_str);
  ulong_to_str(nreads, nreads_str);

  status("%s size distribution: "
         "min: %zu mean: %.1f median: %.1f mode: %zu max: %zu; n=%s / %s [%zu%%]",
         insert_sizes ? "Insert" : "Seq error gap",
         min, mean, median, mode, max, ngaps_str, nreads_str,
         (size_t)((100.0 * ngaps) / nreads + 0.5));

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
