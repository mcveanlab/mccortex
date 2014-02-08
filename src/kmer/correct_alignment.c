#include "global.h"
#include "correct_alignment.h"

#define INIT_BUFLEN 1024

size_t correct_aln_worker_est_mem(const dBGraph *graph) {
  return 2*graph_walker_est_mem() + 2*rpt_walker_est_mem(graph->ht.capacity, 22) +
         db_alignment_est_mem() + 2*INIT_BUFLEN*sizeof(dBNode) +
         2*INIT_BUFLEN*sizeof(size_t) + sizeof(CorrectAlnWorker);
}

void correct_aln_worker_alloc(CorrectAlnWorker *wrkr, const dBGraph *db_graph)
{
  CorrectAlnWorker tmp = {.db_graph = db_graph, .aln = NULL,
                          .start_idx = 0, .gap_idx = 0, .end_idx = 0};

  // Graph traversal
  graph_walker_alloc(&tmp.wlk);
  graph_walker_alloc(&tmp.wlk2);
  rpt_walker_alloc(&tmp.rptwlk, db_graph->ht.capacity, 22); // 4MB
  rpt_walker_alloc(&tmp.rptwlk2, db_graph->ht.capacity, 22); // 4MB

  // Node buffers
  db_node_buf_alloc(&tmp.contig, INIT_BUFLEN);
  db_node_buf_alloc(&tmp.revcontig, INIT_BUFLEN);

  // Gap Histogram
  tmp.histgrm_len = INIT_BUFLEN;
  tmp.gap_ins_histgrm = calloc2(tmp.histgrm_len, sizeof(*tmp.gap_ins_histgrm));
  tmp.gap_err_histgrm = calloc2(tmp.histgrm_len, sizeof(*tmp.gap_err_histgrm));

  memcpy(wrkr, &tmp, sizeof(CorrectAlnWorker));
}

void correct_aln_worker_dealloc(CorrectAlnWorker *wrkr)
{
  graph_walker_dealloc(&wrkr->wlk);
  graph_walker_dealloc(&wrkr->wlk2);
  rpt_walker_dealloc(&wrkr->rptwlk);
  rpt_walker_dealloc(&wrkr->rptwlk2);
  db_node_buf_dealloc(&wrkr->contig);
  db_node_buf_dealloc(&wrkr->revcontig);
  free(wrkr->gap_ins_histgrm);
  free(wrkr->gap_err_histgrm);
}

static inline void worker_gap_cap(CorrectAlnWorker *wrkr, size_t max_gap)
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

void correct_alignment_init(CorrectAlnWorker *wrkr, const dBAlignment *aln,
                            CorrectAlnParam params)
{
  // Copy input
  wrkr->aln = aln;
  wrkr->params = params;

  // reset state
  wrkr->start_idx = wrkr->prev_start_idx = 0;
  wrkr->gap_idx = wrkr->end_idx = db_alignment_next_gap(aln, 0);
}

static void prime_for_traversal(GraphWalker *wlk,
                                const dBNode *block, size_t n,
                                size_t max_context, boolean forward,
                                size_t ctxcol, size_t ctpcol,
                                const dBGraph *db_graph)
{
  ctx_assert(n > 0);
  dBNode node0;

  // printf("Prime [%zu]: ", n);
  // db_nodes_print(block, n, db_graph, stdout);
  // printf("\n");

  db_node_check_nodes(block, n, db_graph);

  if(n > max_context) {
    if(forward) block = block + n - max_context;
    n = max_context;
  }

  if(forward) { node0 = block[0]; block++; }
  else { node0 = db_node_reverse(block[n-1]); }

  // status("prime_for_traversal() %zu:%i n=%zu %s", (size_t)node0.key, node0.orient,
  //        n, forward ? "forward" : "reverse");

  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node0);
  // graph_walker_fast_traverse(wlk, block, n-1, forward);
  graph_walker_slow_traverse(wlk, block, n-1, forward);


  // For debugging
  // graph_walker_print_state(wlk);
}

// Resets node buffer, GraphWalker and RepeatWalker after use
static boolean traverse_one_way2(dBNode end_node, dBNodeBuffer *contig,
                                 size_t gap_min, size_t gap_max,
                                 size_t *gap_len_ptr,
                                 GraphWalker *wlk, RepeatWalker *rptwlk)
{
  boolean traversed = false;

  BinaryKmer tmpbkmer = db_node_get_bkmer(wlk->db_graph, end_node.key);
  char tmpstr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(tmpbkmer, wlk->db_graph->kmer_size, tmpstr);
  // status("Endnode: %s contig->len: %zu", tmpstr, contig->len);

  size_t init_len = contig->len, max_len = contig->len + gap_max;
  db_node_buf_ensure_capacity(contig, contig->len + gap_max + 1);

  while(contig->len < max_len && graph_traverse(wlk) &&
        rpt_walker_attempt_traverse(rptwlk, wlk))
  {
    if(db_nodes_match(wlk->node, end_node)) {
      traversed = true;
      break;
    }
    contig->data[contig->len++] = wlk->node;
  }

  size_t gap_len = contig->len - init_len;

  // Clean up GraphWalker, RepeatWalker
  graph_walker_finish(wlk);
  rpt_walker_fast_clear(rptwlk, contig->data+init_len, gap_len);

  // Fail if we bridged but too short
  if(gap_len < gap_min) traversed = false;
  if(!traversed) contig->len = init_len;

  // status("gap_len: %zu contig->len: %zu", gap_len, contig->len);

  *gap_len_ptr = gap_len;
  return traversed;
}

// Traversed from both sides of the gap
// nodes in contig1 will be the *reverse*
// Resets node buffers, GraphWalkers and RepeatWalkers after use
static boolean traverse_two_way2(dBNodeBuffer *contig0, dBNodeBuffer *contig1,
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
  //   bkmer = db_node_get_bkmer(wlk0->db_graph, wlk[i]->node);
  //   binary_kmer_to_str(bkmer, wlk0->db_graph->kmer_size, str);
  //   status("primed %zu: %zu:%i %s\n", i, (size_t)wlk[i]->node, (int)wlk[i]->orient, str);
  // }
  //

  while(gap_len < gap_max && (use[0] || use[1])) {
    for(i = 0; i < 2; i++) {
      if(use[i]) {
        use[i] = (graph_traverse(wlk[i]) &&
                  rpt_walker_attempt_traverse(rptwlk[i], wlk[i]));

        if(use[i]) {
          //
          // bkmer = db_node_get_bkmer(wlk0->db_graph, wlk[i]->node);
          // binary_kmer_to_str(bkmer, wlk0->db_graph->kmer_size, str);
          // status("%zu: %zu:%i %s\n", i, (size_t)wlk[i]->node, (int)wlk[i]->orient, str);
          //

          if(wlk[0]->node.key == wlk[1]->node.key &&
             wlk[0]->node.orient != wlk[1]->node.orient) {
            traversed = true;
            use[0] = use[1] = false; // set both to false to exit loop
            break;
          }

          contig[i]->data[contig[i]->len++] = wlk[i]->node;
          gap_len++;
        }
      }
    }
  }

  // status("traverse_two_way2() gap_len:%zu gap_min:%zu gap_max:%zu %s",
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

static boolean traverse_one_way(CorrectAlnWorker *wrkr,
                                size_t gap_idx, size_t end_idx,
                                size_t gap_min, size_t gap_max,
                                size_t *gap_len_ptr)
{
  const CorrectAlnParam *params = &wrkr->params;
  const dBNode *aln_nodes = wrkr->aln->nodes.data;
  const dBNodeBuffer *nbuf = &wrkr->contig;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t block1len = end_idx-gap_idx;
  const size_t ctxcol = params->ctxcol, ctpcol = params->ctpcol;
  boolean traversed;


  // Start traversing forward
  prime_for_traversal(&wrkr->wlk, nbuf->data, nbuf->len,
                      params->max_context, true, ctxcol, ctpcol, db_graph);

  traversed
    = traverse_one_way2(aln_nodes[gap_idx],
                        &wrkr->contig, gap_min, gap_max, gap_len_ptr,
                        &wrkr->wlk, &wrkr->rptwlk);

  // status("We had %s %zu", traversed ? "Success" : "Failure", wrkr->contig.len);

  if(traversed) return true;

  // Start traversing backwards
  prime_for_traversal(&wrkr->wlk, aln_nodes+gap_idx, block1len,
                      params->max_context, false, ctxcol, ctpcol, db_graph);

  traversed
    = traverse_one_way2(db_node_reverse(aln_nodes[gap_idx-1]),
                        &wrkr->revcontig, gap_min, gap_max, gap_len_ptr,
                        &wrkr->wlk, &wrkr->rptwlk);

  return traversed;
}

static boolean traverse_two_way(CorrectAlnWorker *wrkr,
                                size_t gap_idx, size_t end_idx,
                                size_t gap_min, size_t gap_max,
                                size_t *gap_len_ptr)
{
  const CorrectAlnParam *params = &wrkr->params;
  const dBNode *nodes = wrkr->aln->nodes.data;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t block1len = end_idx-gap_idx;
  const size_t ctxcol = params->ctxcol, ctpcol = params->ctpcol;

  prime_for_traversal(&wrkr->wlk, wrkr->contig.data, wrkr->contig.len,
                      params->max_context, true, ctxcol, ctpcol, db_graph);
  prime_for_traversal(&wrkr->wlk2, nodes+gap_idx, block1len,
                      params->max_context, false, ctxcol, ctpcol, db_graph);

  // status("Traversing two-way... %zu[len:%zu]:%zu[len:%zu]",
  //        start_idx, block0len, gap_idx, block1len);

  return traverse_two_way2(&wrkr->contig, &wrkr->revcontig,
                           gap_min, gap_max, gap_len_ptr,
                           &wrkr->wlk, &wrkr->rptwlk,
                           &wrkr->wlk2, &wrkr->rptwlk2);
}

// Returns NULL if end of alignment
dBNodeBuffer* correct_alignment_nxt(CorrectAlnWorker *wrkr)
{
  if(wrkr->start_idx == wrkr->end_idx) return NULL;

  const dBAlignment *aln = wrkr->aln;
  const dBNodeBuffer *nodes = &aln->nodes;
  const uint32Buffer *gaps = &aln->gaps;
  const size_t num_align_nodes = nodes->len;

  const CorrectAlnParam *const params = &wrkr->params;
  const float gap_variance = params->gap_variance;
  const size_t gap_wiggle = params->gap_wiggle;
  const size_t ins_gap_min = params->ins_gap_min;
  const size_t ins_gap_max = params->ins_gap_max;

  // worker_generate_contigs ensures contig is at least nodes->len long
  boolean both_reads = (aln->used_r1 && aln->used_r2);
  size_t i, j, block0len, block1len;
  size_t gap_est, gap_len, gap_min, gap_max;
  boolean traversed;

  dBNodeBuffer *contig = &wrkr->contig, *revcontig = &wrkr->revcontig;

  block0len = wrkr->gap_idx - wrkr->start_idx;

  db_node_buf_reset(contig);
  db_node_buf_append(contig, nodes->data+wrkr->start_idx, block0len);

  while(wrkr->end_idx < num_align_nodes)
  {
    // block0 [start_idx..gap_idx-1], block1 [gap_idx..end_idx]
    wrkr->end_idx = db_alignment_next_gap(aln, wrkr->end_idx);
    block0len = wrkr->gap_idx - wrkr->start_idx;
    block1len = wrkr->end_idx - wrkr->gap_idx;

    // We've got a gap to traverse
    long gap_min_long, gap_max_long;
    boolean is_mp;

    // Get bound for acceptable bridge length (min,max length values)
    is_mp = (both_reads && wrkr->gap_idx == aln->r2strtidx);
    gap_est = gaps->data[wrkr->gap_idx];
    if(is_mp) gap_est += aln->r1enderr;

    gap_min_long = (long)gap_est - (long)(gap_est * gap_variance + gap_wiggle);
    gap_max_long = (long)gap_est + (long)(gap_est * gap_variance + gap_wiggle);

    if(is_mp) {
      gap_min_long += ins_gap_min;
      gap_max_long += ins_gap_max;
    }

    gap_min = (size_t)MAX2(0, gap_min_long);
    gap_max = (size_t)MAX2(0, gap_max_long);

    db_node_buf_reset(revcontig);

    // Alternative traversing from both sides
    if(params->one_way_gap_traverse)
      traversed = traverse_one_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                   gap_min, gap_max, &gap_len);
    else
      traversed = traverse_two_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                   gap_min, gap_max, &gap_len);

    // status("traversal: %s!\n", traversed ? "worked" : "failed");

    if(!traversed) break;

    // reverse and copy from revcontig -> contig
    db_node_buf_ensure_capacity(contig, contig->len + revcontig->len);
    for(i = 0, j = contig->len+revcontig->len-1; i < revcontig->len; i++, j--)
      contig->data[j] = db_node_reverse(revcontig->data[i]);

    contig->len += revcontig->len;

    // Copy block1
    db_node_buf_append(contig, nodes->data+wrkr->gap_idx, block1len);

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

    wrkr->gap_idx = wrkr->end_idx;
  }

  // status("start: %zu gap: %zu end: %zu",
  //        wrkr->start_idx, wrkr->gap_idx, wrkr->end_idx);
  // db_nodes_print(wrkr->contig.data, wrkr->contig.len, wrkr->db_graph, stdout);
  // printf("\n");

  wrkr->prev_start_idx = wrkr->start_idx;
  wrkr->start_idx = wrkr->gap_idx;
  wrkr->gap_idx = wrkr->end_idx;

  return &wrkr->contig;
}


size_t correct_alignment_get_strtidx(CorrectAlnWorker *wrkr) {
  return wrkr->prev_start_idx;
}

size_t correct_alignment_get_endidx(CorrectAlnWorker *wrkr) {
  return wrkr->gap_idx;
}

uint64_t* correct_alignment_get_errhist(CorrectAlnWorker *wrkr, size_t *n) {
  *n = wrkr->histgrm_len;
  return wrkr->gap_err_histgrm;
}

uint64_t* correct_alignment_get_inshist(CorrectAlnWorker *wrkr, size_t *n) {
  *n = wrkr->histgrm_len;
  return wrkr->gap_ins_histgrm;
}


static inline void merge_arrays(uint64_t *restrict a,
                                const uint64_t *restrict b,
                                size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) a[i] += b[i];
}

// copy to dst histrograms, zero src histograms
void correct_alignment_merge_hists(CorrectAlnWorker *dst, CorrectAlnWorker *src)
{
  ctx_assert(dst != src);
  worker_gap_cap(dst, src->histgrm_len);
  merge_arrays(dst->gap_err_histgrm, src->gap_err_histgrm, src->histgrm_len);
  merge_arrays(dst->gap_ins_histgrm, src->gap_ins_histgrm, src->histgrm_len);

  // Zero arrays
  memset(src->gap_err_histgrm, 0, src->histgrm_len * sizeof(uint64_t));
  memset(src->gap_ins_histgrm, 0, src->histgrm_len * sizeof(uint64_t));
}
