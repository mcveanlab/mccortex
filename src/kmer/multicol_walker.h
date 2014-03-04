#ifndef MULTICOL_H_
#define MULTICOL_H_

#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "repeat_walker.h"

typedef struct
{
  const dBGraph *db_graph;
  GraphWalker *wlks;
  RepeatWalker rptwlk;
} MulticolWalker;

static inline void multicol_walker_alloc(MulticolWalker *walker,
                                         const dBGraph *db_graph)
{
  size_t i;
  *walker = (MulticolWalker){.db_graph = db_graph};
  walker->wlks = malloc2(db_graph->num_of_cols * sizeof(GraphWalker));
  for(i = 0; i < db_graph->num_of_cols; i++)
    graph_walker_alloc(&walker->wlks[i]);
  rpt_walker_alloc(&walker->rptwlk, db_graph->ht.capacity, 22); // use 4MB
}

static inline void multicol_walker_dealloc(MulticolWalker *walker)
{
  size_t i, num_of_cols = walker->db_graph->num_of_cols;
  for(i = 0; i < num_of_cols; i++)
    graph_walker_dealloc(&walker->wlks[i]);
  free(walker->wlks);
  rpt_walker_dealloc(&walker->rptwlk);
}

size_t multicol_walker_assemble_contig(MulticolWalker *walker,
                                       const size_t *colours, size_t num_cols,
                                       size_t *cols_used, size_t *col_lengths,
                                       dBNodeBuffer *nbuf);

// Returns number of remaining colours
size_t multicol_walker_rem_cols(size_t *colours, size_t num_cols,
                                const size_t *used_cols, size_t num_used);

#endif /* MULTICOL_H_ */
