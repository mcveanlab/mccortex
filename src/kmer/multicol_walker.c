#include "global.h"
#include "multicol_walker.h"

// Returns number of colours in the contig
size_t multicol_walker_assemble_contig(MulticolWalker *walker,
                                       const size_t *colours, size_t num_cols,
                                       size_t *cols_used, size_t *col_lengths,
                                       dBNodeBuffer *nbuf)
{
  ctx_assert(nbuf->len > 0);

  // colours[0] is the primary colour
  const dBGraph *db_graph = walker->db_graph;
  size_t i, j;
  int last_idx;

  for(i = 0; i < num_cols; i++) {
    graph_walker_init(&walker->wlks[i], walker->db_graph, colours[i], colours[i],
                      nbuf->data[0]);
    graph_walker_slow_traverse(&walker->wlks[i], nbuf->data+1, nbuf->len-1, true);
  }

  // GO!
  GraphWalker *wlks = walker->wlks;
  RepeatWalker *rptwlk = &walker->rptwlk;
  size_t num_remaining_wlkrs, num_forking_cols = 0, num_finished_cols = 0;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  size_t num_next;
  BinaryKmer bkmer;
  Edges edges;
  dBNode node = nbuf->data[nbuf->len-1];
  uint32_t first_finished;

  num_remaining_wlkrs = num_cols;

  while(num_remaining_wlkrs > 0)
  {
    bkmer = db_node_get_bkmer(db_graph, node.key);
    edges = db_node_get_edges(db_graph, node.key, 0);
    num_next = db_graph_next_nodes(db_graph, bkmer, node.orient, edges,
                                   next_nodes, next_nucs);

    first_finished = UINT32_MAX;
    last_idx = -1;

    for(i = 0; i < num_remaining_wlkrs; i++)
    {
      if(!graph_traverse_nodes(&wlks[i], num_next, next_nodes, next_nucs))
      {
        // colour finished
        cols_used[num_finished_cols] = wlks[i].ctxcol;
        col_lengths[num_finished_cols] = nbuf->len;
        first_finished = MIN2(i, first_finished);
        num_finished_cols++;
      }
      else if(last_idx == -1) last_idx = wlks[i].last_step.idx;
      else if(wlks[i].last_step.idx != last_idx) {
        // colour diverged
        first_finished = MIN2(i, first_finished);
        num_forking_cols++;
      }
    }

    if(last_idx == -1) { break; }
    else if(first_finished < UINT32_MAX) {
      for(i = first_finished+1, j = first_finished; i < num_remaining_wlkrs; i++) {
        if(wlks[i].last_step.idx != last_idx)
          graph_walker_finish(&walker->wlks[i]);
        else
          memcpy(&wlks[j++], &wlks[i], sizeof(GraphWalker));
      }
      num_remaining_wlkrs = j;
    }

    if(!rpt_walker_attempt_traverse(rptwlk, &wlks[0])) {
      break;
    }

    db_node_buf_add(nbuf, wlks[0].node);
  }

  for(i = 0; i < num_remaining_wlkrs; i++)
    graph_walker_finish(&walker->wlks[i]);

  rpt_walker_fast_clear(rptwlk, nbuf->data, nbuf->len);

  return num_remaining_wlkrs + num_finished_cols;
}

// Returns number of remaining colours
size_t multicol_walker_rem_cols(size_t *colours, size_t num_cols,
                                const size_t *used_cols, size_t num_used)
{
  size_t i, j, k;
  for(i = j = k = 0; i < num_cols && j < num_used; i++) {
    if(colours[i] == used_cols[j]) j++;
    else colours[k++] = colours[i];
  }
  // Get remaining colours
  while(i < num_cols) colours[k++] = colours[i++];
  return k;
}
