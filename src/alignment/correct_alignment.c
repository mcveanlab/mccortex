#include "global.h"
#include "correct_alignment.h"

// Global variables to specify if we should print output - used for debugging only
// These are used in generate_paths.c
bool gen_paths_print_contigs = false, gen_paths_print_paths = false;
bool gen_paths_print_reads = false;

#define INIT_BUFLEN 1024

size_t correct_aln_worker_est_mem(const dBGraph *graph) {
  return 2*graph_walker_est_mem() + 2*rpt_walker_est_mem(graph->ht.capacity, 22) +
         db_alignment_est_mem() + 2*INIT_BUFLEN*sizeof(dBNode) +
         2*INIT_BUFLEN*sizeof(size_t) + sizeof(CorrectAlnWorker);
}

void correct_aln_worker_alloc(CorrectAlnWorker *wrkr, bool store_contig_lens,
                              const dBGraph *db_graph)
{
  CorrectAlnWorker tmp = {.db_graph = db_graph,
                          .start_idx = 0, .gap_idx = 0, .end_idx = 0,
                          .store_contig_lens = store_contig_lens};

  // Graph alignment of reads
  db_alignment_alloc(&tmp.aln);

  // Graph traversal
  graph_walker_alloc(&tmp.wlk);
  graph_walker_alloc(&tmp.wlk2);
  rpt_walker_alloc(&tmp.rptwlk, db_graph->ht.capacity, 22); // 4MB
  rpt_walker_alloc(&tmp.rptwlk2, db_graph->ht.capacity, 22); // 4MB

  // Node buffers
  db_node_buf_alloc(&tmp.contig, INIT_BUFLEN);
  db_node_buf_alloc(&tmp.revcontig, INIT_BUFLEN);
  int32_buf_alloc(&tmp.rpos, INIT_BUFLEN);

  correct_aln_stats_alloc(&tmp.aln_stats);
  loading_stats_init(&tmp.load_stats);

  memcpy(wrkr, &tmp, sizeof(CorrectAlnWorker));
}

void correct_aln_worker_dealloc(CorrectAlnWorker *wrkr)
{
  db_alignment_dealloc(&wrkr->aln);
  graph_walker_dealloc(&wrkr->wlk);
  graph_walker_dealloc(&wrkr->wlk2);
  rpt_walker_dealloc(&wrkr->rptwlk);
  rpt_walker_dealloc(&wrkr->rptwlk2);
  db_node_buf_dealloc(&wrkr->contig);
  db_node_buf_dealloc(&wrkr->revcontig);
  int32_buf_dealloc(&wrkr->rpos);
  correct_aln_stats_dealloc(&wrkr->aln_stats);
}

/*!
  @param params Settings for correction - needs to be passed since we don't know
                which source the reads came from and diff input sources have
                different requirements (e.g. expected insert size)
 */
void correct_alignment_init(CorrectAlnWorker *wrkr,
                            const CorrectAlnParam *params,
                            const read_t *r1, const read_t *r2,
                            uint8_t fq_cutoff1, uint8_t fq_cutoff2,
                            int8_t hp_cutoff)
{
  ctx_assert(params->ctxcol == params->ctpcol);

  db_alignment_from_reads(&wrkr->aln, r1, r2,
                          fq_cutoff1, fq_cutoff2, hp_cutoff,
                          wrkr->db_graph, params->ctxcol);

  const dBAlignment *aln = &wrkr->aln;

  // Copy parameters
  wrkr->params = *params;

  // reset state
  wrkr->start_idx = wrkr->prev_start_idx = 0;
  wrkr->gap_idx = wrkr->end_idx = db_alignment_next_gap(aln, 0);

  // Update stats
  wrkr->load_stats.total_bases_read += aln->r1bases + aln->r2bases;
  wrkr->load_stats.num_kmers_parsed += aln->nodes.len;

  if(aln->passed_r2) wrkr->load_stats.num_pe_reads += 2;
  else               wrkr->load_stats.num_se_reads++;

  wrkr->aln_stats.num_ins_gaps += aln->used_r2;
}

// Merge stats into dst and reset src
void correct_aln_merge_stats(CorrectAlnWorker *restrict dst,
                             CorrectAlnWorker *restrict src)
{
  correct_aln_stats_merge(&dst->aln_stats, &src->aln_stats);
  correct_aln_stats_dealloc(&src->aln_stats);
  loading_stats_merge(&dst->load_stats, &src->load_stats);
  loading_stats_init(&src->load_stats);
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
  db_node_buf_capacity(contig, max_len);

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
  db_node_buf_capacity(contig0, contig0->len + gap_max + 1);
  db_node_buf_capacity(contig1, contig1->len + gap_max + 1);

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

        if(db_nodes_are_equal(nodes[0], db_node_reverse(nodes[1])))
        {
          // gap_len may now be > gap_max since we loop twice before
          // checking outer loop condition
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
  const int aln_colour = wrkr->aln.colour; // -1 for all
  const dBNode *aln_nodes = wrkr->aln.nodes.data;
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

  correct_aln_stats_update(&wrkr->aln_stats, result);

  if(result.traversed) return result;

  // right to left
  graph_walker_prime(&wrkr->wlk, aln_nodes+gap_idx, block1len,
                     params->max_context, false, ctxcol, ctpcol, db_graph);

  ctx_assert(wrkr->revcontig.len == 0);

  result = traverse_one_way2(nbuf->data, nbuf->len, false,
                             &wrkr->revcontig, gap_min, gap_max,
                             &wrkr->wlk, &wrkr->rptwlk,
                             only_in_one_col, params->use_end_check);

  correct_aln_stats_update(&wrkr->aln_stats, result);

  return result;
}

static TraversalResult traverse_two_way(CorrectAlnWorker *wrkr,
                                        size_t gap_idx, size_t end_idx,
                                        size_t gap_min, size_t gap_max)
{
  const CorrectAlnParam *params = &wrkr->params;
  const int aln_colour = wrkr->aln.colour; // -1 for all
  const dBNode *nodes = wrkr->aln.nodes.data;
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

  correct_aln_stats_update(&wrkr->aln_stats, result);

  return result;
}

// @return NULL if end of alignment, otherwise returns pointer to wrkr->contig
dBNodeBuffer* correct_alignment_nxt(CorrectAlnWorker *wrkr)
{
  if(wrkr->start_idx == wrkr->end_idx) return NULL;

  const size_t kmer_size = wrkr->db_graph->kmer_size;

  const dBAlignment *aln = &wrkr->aln;
  const CorrectAlnParam params = wrkr->params;
  const size_t num_align_nodes = aln->nodes.len;

  // Get arrays for nodes and kmer positions in alignment
  const dBNode *aln_nodes = aln->nodes.data;
  const int32_t *aln_rpos = aln->rpos.data;

  // worker_generate_contigs ensures contig is at least aln_nodes->len long
  bool both_reads = (aln->used_r1 && aln->used_r2);
  size_t i, j, block0len, block1len;
  size_t gap_est, gap_min, gap_max;

  dBNodeBuffer *contig = &wrkr->contig, *revcontig = &wrkr->revcontig;
  Int32Buffer *contig_rpos = &wrkr->rpos;

  block0len = wrkr->gap_idx - wrkr->start_idx;

  db_node_buf_reset(contig);
  int32_buf_reset(contig_rpos);
  db_node_buf_append(contig,    aln_nodes+wrkr->start_idx, block0len);
  int32_buf_append(contig_rpos, aln_rpos +wrkr->start_idx, block0len);

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
    if(is_mp) {
      ctx_assert(aln->r1enderr == db_aln_r1enderr(aln,kmer_size));
      gap_est = aln->r1enderr + aln_rpos[wrkr->gap_idx];
      wrkr->aln_stats.num_ins_gaps++;
    } else {
      gap_est = aln_rpos[wrkr->gap_idx] - aln_rpos[wrkr->gap_idx-1];
      wrkr->aln_stats.num_mid_gaps++;
    }

    // Wiggle due to variation / error
    long wiggle = gap_est * params.gap_variance + params.gap_wiggle;

    gap_min_long = gap_est - wiggle;
    gap_max_long = gap_est + wiggle;

    // Note: do not apply wiggle to insert gap
    // Work out insert gap kmers using fragment length
    if(is_mp) {
      // Work out insert-gap kmers using fragment length
      size_t sum_read_bases = aln->r1bases + aln->r2bases;
      gap_min_long += params.frag_len_min - sum_read_bases + kmer_size-1;
      gap_max_long += params.frag_len_max - sum_read_bases + kmer_size-1;
    }

    if(gap_max_long < 0) break;

    // Note: we don't deal with negative gaps currently
    gap_min = (size_t)MAX2(0, gap_min_long);
    gap_max = (size_t)MAX2(0, gap_max_long);

    db_node_buf_reset(revcontig);

    // Alternative traversing from both sides
    // gap len is the number of kmers filling the gap
    TraversalResult result;

    if(params.one_way_gap_traverse)
      result = traverse_one_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                gap_min, gap_max);
    else
      result = traverse_two_way(wrkr, wrkr->gap_idx, wrkr->end_idx,
                                gap_min, gap_max);

    // status("traversal: %s!\n", result.traversed ? "worked" : "failed");

    if(!result.traversed) break;

    if(is_mp) wrkr->aln_stats.num_ins_traversed++;
    else      wrkr->aln_stats.num_mid_traversed++;

    // reverse and copy from revcontig -> contig
    db_node_buf_capacity(contig, contig->len + revcontig->len);

    // reverse order and orientation of nodes
    size_t new_contig_len = contig->len + revcontig->len;
    for(i = contig->len, j = revcontig->len-1; i < new_contig_len; i++, j--)
      contig->data[i] = db_node_reverse(revcontig->data[j]);

    contig->len += revcontig->len;

    // Append -1 values to rpos for gap
    size_t len_and_gap = contig_rpos->len + result.gap_len;
    int32_buf_capacity(contig_rpos, len_and_gap);
    while(contig_rpos->len < len_and_gap)
      contig_rpos->data[contig_rpos->len++] = -1;

    // Copy block1, after gap
    db_node_buf_append(contig,    aln_nodes+wrkr->gap_idx, block1len);
    int32_buf_append(contig_rpos, aln_rpos+wrkr->gap_idx,  block1len);

    // Update gap stats
    if(is_mp) {
      correct_aln_stats_add_mp(&wrkr->aln_stats, result.gap_len,
                               aln->r1bases, aln->r2bases, kmer_size);
    } else {
      correct_aln_stats_add(&wrkr->aln_stats, gap_est, result.gap_len);
    }

    wrkr->gap_idx = wrkr->end_idx;
  }

  wrkr->prev_start_idx = wrkr->start_idx;
  wrkr->start_idx = wrkr->gap_idx;
  wrkr->gap_idx = wrkr->end_idx;

  // db_nodes_print_verbose(contig->data, contig->len, wrkr->db_graph, stdout);
  ctx_check(db_node_check_nodes(contig->data, contig->len, wrkr->db_graph));

  size_t contig_length_bp = contig->len + kmer_size - 1;

  wrkr->load_stats.contigs_parsed++;
  wrkr->load_stats.num_kmers_loaded += contig->len;
  wrkr->load_stats.total_bases_loaded += contig_length_bp;

  if(wrkr->store_contig_lens)
    correct_aln_stats_add_contig(&wrkr->aln_stats, contig_length_bp);

  return contig;
}

/*!
  Correct a whole read, filling in gaps caused by sequencing error with the
  graph
  @param wrkr     Initialised temporary memory to use in doing alignment
  @param r        Read to align to the graph
  @param nodebuf  Store nodebuf from read and inferred in gaps in buffer
  @param posbuf   Positions in the read of kmers (-1 means inferred from graph)
 */
void correct_aln_read(CorrectAlnWorker *wrkr, const CorrectAlnParam *params,
                      const read_t *r, uint8_t fq_cutoff, uint8_t hp_cutoff,
                      dBNodeBuffer *nodebuf, Int32Buffer *posbuf)
{
  size_t i, j;

  db_node_buf_reset(nodebuf);
  int32_buf_reset(posbuf);

  db_node_buf_capacity(nodebuf, r->seq.end);
  int32_buf_capacity(posbuf, r->seq.end);

  // Correct sequence errors in the alignment
  correct_alignment_init(wrkr, params, r, NULL, fq_cutoff, 0, hp_cutoff);

  const dBAlignment *aln = &wrkr->aln;

  // All or no nodes aligned
  if(db_alignment_is_perfect(aln) || aln->nodes.len == 0) {
    memcpy(nodebuf->data, aln->nodes.data, aln->nodes.len * sizeof(dBNode));
    for(i = 0; i < aln->nodes.len; i++) posbuf->data[i] = i;
    nodebuf->len = posbuf->len = aln->nodes.len;
    return;
  }

  size_t left_gap = aln->rpos.data[0], right_gap = aln->r1enderr;

  // Temporary walkers used to traverse the graph
  GraphWalker *wlk = &wrkr->wlk;
  RepeatWalker *rptwlk = &wrkr->rptwlk;

  // Get first alignment, if there are gaps, backtrack
  // nbuf = correct_alignment_nxt(wrkr);
  // ctx_assert(nbuf != NULL); // We've already checked aln->nodes.len > 0

  // Append contigs
  while(correct_alignment_nxt(wrkr) != NULL)
  {
    // printf("got thing!\n");
    ctx_assert(wrkr->contig.len ==  wrkr->rpos.len);
    db_node_buf_append(nodebuf, wrkr->contig.data, wrkr->contig.len);
    int32_buf_append(posbuf, wrkr->rpos.data, wrkr->rpos.len);
  }

  // Couldn't get any kmers
  if(nodebuf->len == 0) {
    // printf("\nGot no kmers!\n\n");
    return;
  }

  if(left_gap > 0)
  {
    // Extend to the left
    wrkr->aln_stats.num_end_gaps++;

    dBNodeBuffer *revcontig = &wrkr->revcontig;
    db_node_buf_reset(revcontig);
    db_node_buf_capacity(revcontig, left_gap);

    // Try to fill in missing kmers
    size_t n = 1;
    while(n < nodebuf->len && posbuf->data[n] == posbuf->data[n-1]+1) { n++; }

    graph_walker_prime(wlk, nodebuf->data, n,
                       params->max_context, false,
                       params->ctxcol, params->ctpcol,
                       wrkr->db_graph);

    for(i = 0; i < left_gap && graph_walker_next(wlk) &&
                               rpt_walker_attempt_traverse(rptwlk, wlk);  i++)
    {
      revcontig->data[i] = wlk->node;
    }
    revcontig->len = i;

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, revcontig->data, revcontig->len);

    if(revcontig->len > 0)
      wrkr->aln_stats.num_end_traversed++;

    // Shift away from zero
    db_node_buf_shift_left(nodebuf, revcontig->len);
    int32_buf_shift_left(  posbuf,  revcontig->len);

    // copy in reverse order
    for(i = 0, j = revcontig->len-1; i < revcontig->len; i++, j--) {
      // printf("[add] %zu: %zu\n", j, revcontig->data[j].key);
      nodebuf->data[i] = db_node_reverse(revcontig->data[j]);
      posbuf->data[i] = -1;
    }
  }

  if(right_gap > 0)
  {
    // Extend to the right
    wrkr->aln_stats.num_end_gaps++;

    size_t n = nodebuf->len-1;
    while(n > 0 && posbuf->data[n] == posbuf->data[n-1]+1) { n--; }

    graph_walker_prime(wlk, nodebuf->data + n, nodebuf->len - n,
                       params->max_context, true,
                       params->ctxcol, params->ctpcol,
                       wrkr->db_graph);

    size_t orig_len = nodebuf->len;
    size_t end_len = nodebuf->len + right_gap;
    db_node_buf_capacity(nodebuf, end_len);
    int32_buf_capacity(  posbuf,  end_len);

    for(i = nodebuf->len; i < end_len && graph_walker_next(wlk) &&
                          rpt_walker_attempt_traverse(rptwlk, wlk); i++)
    {
      nodebuf->data[i] = wlk->node;
      posbuf->data[i]  = -1;
    }
    nodebuf->len = posbuf->len = i;

    if(nodebuf->len > orig_len)
      wrkr->aln_stats.num_end_traversed++;

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, nodebuf->data+orig_len, nodebuf->len-orig_len);
  }
}
