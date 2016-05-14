#include "global.h"
#include "unitig_graph.h"
#include "db_node.h"
#include "db_unitig.h"

// Store ends of unitig currently stored in `nodes` and `orients` arrays
size_t unitig_graph_store_end_mt(const dBNode *nodes, size_t n,
                                 UnitigKmerGraph *ugraph)
{
  size_t unitig_id = __sync_fetch_and_add((volatile size_t*)&ugraph->num_unitigs,
                                          1);
  volatile UnitigEnd *unitig_ends = ugraph->unitig_ends;

  ctx_assert(n > 0);
  ctx_assert(unitig_ends[nodes[0].key].assigned == 0);
  ctx_assert(unitig_ends[nodes[n-1].key].assigned == 0);

  /*
  // This is not required by any code atm, commented out for speedup
  // Label intermediate kmers
  for(i = 1; i+1 < n; i++) {
    ctx_assert(unitig_ends[nodes[i].key].assigned == 0);
    unitig_ends[nodes[i].key] = (UnitigEnd){.unitigid = unitig_id,
                                            .assigned = 0,
                                            .left = 0, .right = 0,
                                            .lorient = 0, .rorient = 0};
  }
  */

  UnitigEnd end0 = {.unitigid = unitig_id, .assigned = 1,
                    .left = 1, .right = (n == 1),
                    .lorient = nodes[0].orient,
                    .rorient = nodes[n-1].orient};

  UnitigEnd end1 = {.unitigid = unitig_id, .assigned = 1,
                    .left = (n == 1), .right = 1,
                    .lorient = nodes[0].orient,
                    .rorient = nodes[n-1].orient};

  unitig_ends[nodes[0].key] = end0;
  unitig_ends[nodes[n-1].key] = end1;

  return unitig_id;
}

static void _create_unitig(dBNodeBuffer nbuf, size_t threadid, void *arg)
{
  (void)threadid; // unused
  UnitigKmerGraph *ugraph = (UnitigKmerGraph*)arg;
  db_unitig_normalise(nbuf.b, nbuf.len, ugraph->db_graph);
  size_t uidx = unitig_graph_store_end_mt(nbuf.b, nbuf.len, ugraph);
  if(ugraph->per_untig) {
    ugraph->per_untig(nbuf.b, nbuf.len, uidx, ugraph->per_untig_arg);
  }
}

/**
 * @param visited must be initialised to zero, will be dirty upon return
 **/
void unitig_graph_create(UnitigKmerGraph *ugraph,
                         size_t nthreads,
                         uint8_t *visited,
                         void (*per_untig)(const dBNode *nodes, size_t n,
                                           size_t uidx, void *arg),
                        void *per_untig_arg)
{
  ugraph->per_untig = per_untig;
  ugraph->per_untig_arg = per_untig_arg;
  db_unitigs_iterate(nthreads, visited, ugraph->db_graph, _create_unitig, ugraph);
}

void unitig_graph_alloc(UnitigKmerGraph *ugraph, const dBGraph *db_graph)
{
  UnitigEnd *unitig_ends = ctx_calloc(db_graph->ht.capacity, sizeof(UnitigEnd));
  UnitigKmerGraph tmp = {.unitig_ends = unitig_ends,
                         .num_unitigs = 0,
                         .db_graph = db_graph,
                         .per_untig = NULL,
                         .per_untig_arg = NULL};
  memcpy(ugraph, &tmp, sizeof(tmp));
}

void unitig_graph_dealloc(UnitigKmerGraph *ugraph)
{
  ctx_free(ugraph->unitig_ends);
}
