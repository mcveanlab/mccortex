#include "global.h"
#include "correct_alignment.h"

// Global variables to specify if we should print output - used for debugging only
// These are used in generate_paths.c
bool gen_paths_print_contigs = false, gen_paths_print_paths = false;
bool gen_paths_print_reads = false;


typedef struct {
  uint32_t gap_len;
  uint8_t traversed, paths_disagreed, gap_too_short;
} TraversalResult;

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

  correct_aln_stats_alloc(&tmp.gapstats);

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
  correct_aln_stats_dealloc(&wrkr->gapstats);
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

static inline void correct_aln_stats_update(CorrectAlnStats *stats,
                                            TraversalResult result)
{
  stats->num_gap_attempts++;
  stats->num_gap_successes += result.traversed;
  stats->num_paths_disagreed += result.paths_disagreed;
  stats->num_gaps_too_short += result.gap_too_short;
}

// block is nodes that we are walking towards
// if forward is true, we are walking left to right, otherwise right to left
// Resets node buffer, GraphWalker and RepeatWalker after use
static TraversalResult traverse_one_way2(const dBNode *block, size_t n,
                                         bool forward, dBNodeBuffer *contig,
                                         size_t gap_min, size_t gap_max,
                                         GraphWalker *wlk, RepeatWalker *rptwlk,
                                         bool only_in_col, bool do_paths_check)
{
  dBNode end_node = forward ? block[0] : db_node_reverse(block[n-1]);

  // BinaryKmer tmpbkmer = db_node_get_bkmer(wlk->db_graph, end_node.key);
  // char tmpstr[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(tmpbkmer, wlk->db_graph->kmer_size, tmpstr);
  // status("Endnode: %s contig->len: %zu", tmpstr, contig->len);

  // max_len allows for node on other side of gap
  size_t init_len = contig->len, max_len = contig->len + gap_max + 1;
  db_node_buf_ensure_capacity(contig, max_len);

  TraversalResult result = {.traversed = false, .paths_disagreed = false,
                            .gap_too_short = false, .gap_len = 0};

  // DEBUG
  // fprintf(stderr, "\ntraverse_one_way2:\n\n");
  // db_nodes_print_verbose(contig->data, contig->len, wlk->db_graph, stderr);

  while(contig->len < max_len && graph_walker_next(wlk) &&
        rpt_walker_attempt_traverse(rptwlk, wlk) &&
        (!only_in_col || wlk->last_step.node_has_col))
  {
    // DEBUG
    // db_nodes_print_verbose(&wlk->node, 1, wlk->db_graph, stderr);

    if(db_nodes_are_equal(wlk->node, end_node)) {
      result.traversed = true;
      break;
    }
    contig->data[contig->len++] = wlk->node;
  }

  result.gap_len = contig->len - init_len;
  rpt_walker_fast_clear(rptwlk, contig->data+init_len, result.gap_len);

  // Check paths match remaining nodes
  if(result.traversed && do_paths_check) {
    if(( forward && !graph_walker_agrees_contig(wlk, block+1, n-1, true)) ||
       (!forward && !graph_walker_agrees_contig(wlk, block, n-1, false))) {
      // printf("1) Paths don't agree!\n");
      result.traversed = false;
      result.paths_disagreed = true;
    }
  }

  // printf("gap_len: %zu contig->len: %zu; success: %i gap_min: %zu\n",
  //        result.gap_len, contig->len, (int)traversed, gap_min);

  // Clean up GraphWalker, RepeatWalker
  graph_walker_finish(wlk);

  // Fail if we bridged but too short
  if(result.gap_len < gap_min) {
    result.traversed = false;
    result.gap_too_short = true;
  }

  if(!result.traversed) contig->len = init_len;

  return result;
}

// Traversed from both sides of the gap
// nodes in contig1 will be the *reverse*
// Resets node buffers, GraphWalkers and RepeatWalkers after use
static TraversalResult traverse_two_way2(dBNodeBuffer *contig0,
                                         dBNodeBuffer *contig1,
                                         const dBNode *rhs_block, size_t rhs_n,
                                         size_t gap_min, size_t gap_max,
                                         GraphWalker *wlk0, RepeatWalker *rptwlk0,
                                         GraphWalker *wlk1, RepeatWalker *rptwlk1,
                                         bool only_in_col, bool do_paths_check)
{
  const size_t init_len0 = contig0->len, init_len1 = contig1->len;

  // graph_walker_print_state(wlk0, stdout);
  // graph_walker_print_state(wlk1, stdout);

  // +1 to allow for node on other side of gap
  db_node_buf_ensure_capacity(contig0, contig0->len + gap_max + 1);
  db_node_buf_ensure_capacity(contig1, contig1->len + gap_max + 1);

  bool use[2] = {true,true};
  GraphWalker *wlk[2] = {wlk0, wlk1};
  RepeatWalker *rptwlk[2] = {rptwlk0, rptwlk1};
  dBNodeBuffer *contig[2] = {contig0, contig1};
  dBNode nodes[2] = {wlk0->node, wlk1->node};
  size_t i, gap_len = 0;

  TraversalResult result = {.traversed = false, .paths_disagreed = false,
                            .gap_too_short = false, .gap_len = 0};

  while(gap_len <= gap_max && (use[0] || use[1])) {
    for(i = 0; i < 2; i++) {
      use[i] = (use[i] && graph_walker_next(wlk[i]) &&
                (!only_in_col || wlk[i]->last_step.node_has_col));

      if(use[i])
      {
        // If one of the sides gets stuck in a repeat, we should abandon
        if(!rpt_walker_attempt_traverse(rptwlk[i], wlk[i])) {
          use[0] = use[1] = false; // set both to false to exit loop
          break;
        }

        nodes[i] = wlk[i]->node;

        if(db_nodes_are_equal(nodes[0], db_node_reverse(nodes[1]))) {
          result.traversed = (gap_len <= gap_max);
          use[0] = use[1] = false; // set both to false to exit loop
          break;
        }

        contig[i]->data[contig[i]->len++] = nodes[i];
        gap_len++;
      }
    }
  }

  // printf("2-way Traversal %s\n", result.traversed ? "Worked" : "Failed");

  if(result.traversed && do_paths_check)
  {
    // Check paths match remaining nodes
    dBNode *left_contig = contig0->data; size_t left_n = contig0->len;
    dBNode *right_contig = contig1->data; size_t right_n = contig1->len;

    if(db_nodes_are_equal(wlk[0]->node, db_node_reverse(right_contig[right_n-1]))) {
      right_n--;
    } else {
      ctx_assert(db_nodes_are_equal(db_node_reverse(wlk[1]->node), left_contig[left_n-1]));
      left_n--;
    }

    if(!graph_walker_agrees_contig(wlk[0], right_contig, right_n, false) ||
       !graph_walker_agrees_contig(wlk[0], rhs_block, rhs_n, true) ||
       !graph_walker_agrees_contig(wlk[1], left_contig, left_n, false))
    {
      // printf("2) Paths don't agree!\n");
      result.traversed = false;
      result.paths_disagreed = true;
    }
  }

  // Clear RepeatWalker
  rpt_walker_fast_clear(rptwlk0, contig0->data+init_len0, contig0->len-init_len0);
  rpt_walker_fast_clear(rptwlk1, contig1->data+init_len1, contig1->len-init_len1);

  // Clean up GraphWalker
  graph_walker_finish(wlk0);
  graph_walker_finish(wlk1);

  // Fail if we bridged but too short
  if(gap_len < gap_min) {
    result.traversed = false;
    result.gap_too_short = true;
  }

  if(!result.traversed) {
    contig0->len = init_len0;
    contig1->len = init_len1;
  }

  result.gap_len = gap_len;
  return result;
}

static TraversalResult traverse_one_way(CorrectAlnWorker *wrkr,
                                        size_t gap_idx, size_t end_idx,
                                        size_t gap_min, size_t gap_max)
{
  const CorrectAlnParam *params = &wrkr->params;
  const int aln_colour = wrkr->aln->colour; // -1 for all
  const dBNode *aln_nodes = wrkr->aln->nodes.data;
  const dBNodeBuffer *nbuf = &wrkr->contig;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t block1len = end_idx-gap_idx;
  const size_t ctxcol = params->ctxcol, ctpcol = params->ctpcol;
  TraversalResult result;

  bool only_in_one_col = aln_colour != -1;

  ctx_assert(!only_in_one_col || (size_t)aln_colour == ctxcol);
  ctx_assert(db_nodes_are_equal(aln_nodes[gap_idx-1], nbuf->data[nbuf->len-1]));

  // Start traversing forward
  graph_walker_prime(&wrkr->wlk, nbuf->data, nbuf->len,
                     params->max_context, true, ctxcol, ctpcol, db_graph);

  // left to right
  result = traverse_one_way2(aln_nodes+gap_idx, block1len, true,
                             &wrkr->contig, gap_min, gap_max,
                             &wrkr->wlk, &wrkr->rptwlk,
                             only_in_one_col, params->use_end_check);

  correct_aln_stats_update(&wrkr->gapstats, result);

  if(result.traversed) return result;

  // right to left
  graph_walker_prime(&wrkr->wlk, aln_nodes+gap_idx, block1len,
                     params->max_context, false, ctxcol, ctpcol, db_graph);

  ctx_assert(wrkr->revcontig.len == 0);

  result = traverse_one_way2(nbuf->data, nbuf->len, false,
                             &wrkr->revcontig, gap_min, gap_max,
                             &wrkr->wlk, &wrkr->rptwlk,
                             only_in_one_col, params->use_end_check);

  correct_aln_stats_update(&wrkr->gapstats, result);

  return result;
}

static TraversalResult traverse_two_way(CorrectAlnWorker *wrkr,
                                        size_t gap_idx, size_t end_idx,
                                        size_t gap_min, size_t gap_max)
{
  const CorrectAlnParam *params = &wrkr->params;
  const int aln_colour = wrkr->aln->colour; // -1 for all
  const dBNode *nodes = wrkr->aln->nodes.data;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t block1len = end_idx-gap_idx;
  const size_t ctxcol = params->ctxcol, ctpcol = params->ctpcol;
  TraversalResult result;

  ctx_assert(aln_colour == -1 || (size_t)aln_colour == ctxcol);

  graph_walker_prime(&wrkr->wlk, wrkr->contig.data, wrkr->contig.len,
                     params->max_context, true, ctxcol, ctpcol, db_graph);
  graph_walker_prime(&wrkr->wlk2, nodes+gap_idx, block1len,
                     params->max_context, false, ctxcol, ctpcol, db_graph);

  result = traverse_two_way2(&wrkr->contig, &wrkr->revcontig,
                             nodes+gap_idx, block1len,
                             gap_min, gap_max,
                             &wrkr->wlk, &wrkr->rptwlk,
                             &wrkr->wlk2, &wrkr->rptwlk2,
                             aln_colour != -1, params->use_end_check);

  correct_aln_stats_update(&wrkr->gapstats, result);

  return result;
}

// Returns NULL if end of alignment
dBNodeBuffer* correct_alignment_nxt(CorrectAlnWorker *wrkr)
{
  if(wrkr->start_idx == wrkr->end_idx) return NULL;

  const dBAlignment *aln = wrkr->aln;
  const dBNodeBuffer *nodes = &aln->nodes;
  const Uint32Buffer *gaps = &aln->gaps;
  const size_t num_align_nodes = nodes->len;

  const CorrectAlnParam *const params = &wrkr->params;
  const float gap_variance = params->gap_variance;
  const size_t gap_wiggle = params->gap_wiggle;
  const size_t ins_gap_min = params->ins_gap_min;
  const size_t ins_gap_max = params->ins_gap_max;

  // worker_generate_contigs ensures contig is at least nodes->len long
  bool both_reads = (aln->used_r1 && aln->used_r2);
  size_t i, j, block0len, block1len;
  size_t gap_est, gap_min, gap_max;

  dBNodeBuffer *contig = &wrkr->contig, *revcontig = &wrkr->revcontig;

  block0len = wrkr->gap_idx - wrkr->start_idx;

  db_node_buf_reset(contig);
  db_node_buf_append(contig, nodes->data+wrkr->start_idx, block0len);

  while(wrkr->end_idx < num_align_nodes)
  {
    // block0 [start_idx..gap_idx-1], block1 [gap_idx..end_idx]
    wrkr->end_idx = db_alignment_next_gap(aln, wrkr->end_idx);
    // block0len = wrkr->gap_idx - wrkr->start_idx;
    block1len = wrkr->end_idx - wrkr->gap_idx;

    // We've got a gap to traverse
    long gap_min_long, gap_max_long;
    bool is_mp;

    // Get bound for acceptable bridge length (min,max length values)
    is_mp = (both_reads && wrkr->gap_idx == aln->r2strtidx);

    // gap_est is how many kmers we lost through low qual scores, hp runs etc.
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
    // gap len is the number of kmers filling the gap
    TraversalResult result;

    if(params->one_way_gap_traverse)
      result = traverse_one_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                gap_min, gap_max);
    else
      result = traverse_two_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                gap_min, gap_max);

    // status("traversal: %s!\n", result.traversed ? "worked" : "failed");

    if(!result.traversed) break;

    // reverse and copy from revcontig -> contig
    db_node_buf_ensure_capacity(contig, contig->len + revcontig->len);

    // reverse order and orientation of nodes
    for(i = 0, j = contig->len+revcontig->len-1; i < revcontig->len; i++, j--) {
      contig->data[j] = db_node_reverse(revcontig->data[i]);
    }

    contig->len += revcontig->len;

    // Copy block1
    db_node_buf_append(contig, nodes->data+wrkr->gap_idx, block1len);

    // Update gap stats
    correct_aln_stats_cap(&wrkr->gapstats, result.gap_len);

    // gap_est is the sequence gap (number of missing kmers)
    if(is_mp) {
      size_t ins_gap = (size_t)MAX2((long)result.gap_len - (long)gap_est, 0);
      wrkr->gapstats.gap_ins_histgrm[ins_gap]++;
      // size_t err_gap = result.gap_len - ins_gap;
      // if(err_gap) wrkr->gapstats.gap_err_histgrm[err_gap]++;
    } else {
      wrkr->gapstats.gap_err_histgrm[result.gap_len]++;
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

  // db_nodes_print_verbose(contig->data, contig->len, wrkr->db_graph, stdout);
  ctx_check(db_node_check_nodes(contig->data, contig->len, wrkr->db_graph));

  return contig;
}

// Called after correct_alignment_nxt()
size_t correct_alignment_get_strtidx(CorrectAlnWorker *wrkr) {
  return wrkr->prev_start_idx;
}

// Called after correct_alignment_nxt()
size_t correct_alignment_get_endidx(CorrectAlnWorker *wrkr) {
  return wrkr->start_idx; // start has moved on, so is now where gap_idx was
}
