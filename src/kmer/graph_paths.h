#ifndef GRAPH_PATH_H_
#define GRAPH_PATH_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"

typedef struct {
  Colour *ctxcols, *ctpcols;
  size_t n;
} GraphPathPairing;

static inline void gp_alloc(GraphPathPairing *gp, size_t n) {
  gp->ctxcols = malloc2(2*n*sizeof(*gp->ctxcols));
  gp->ctpcols = gp->ctxcols+n;
  gp->n = n;
}

static inline void gp_dealloc(GraphPathPairing *gp) {
  free(gp->ctxcols);
  gp->n = 0; gp->ctxcols = gp->ctpcols = NULL;
}

// Returns true if new to colour, false otherwise
// packed points to <PathLen><PackedSeq>
// Returns address of path in PathStore by setting newidx
boolean graph_paths_find_or_add_mt(hkey_t hkey, dBGraph *db_graph, Colour ctpcol,
                                   const uint8_t *packed, size_t plen,
                                   PathIndex *newidx);

//
// Functions on graph+paths
//

// col is graph colour
// packed is just <PackedBases>
void graph_path_check_valid(dBNode node, size_t ctxcol, const uint8_t *packed,
                            size_t nbases, const dBGraph *db_graph);

void graph_paths_check_all_paths(const GraphPathPairing *gp,
                                 const dBGraph *db_graph);

void graph_path_check_path(hkey_t node, PathIndex pindex,
                           const GraphPathPairing *gp,
                           const dBGraph *db_graph);

//
// Check
//
// For debugging
void graph_paths_check_counts(const dBGraph *db_graph);

#endif /* GRAPH_PATH_H_ */
