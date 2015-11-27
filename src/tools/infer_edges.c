#include "global.h"
#include "infer_edges.h"
#include "db_node.h"
#include "db_graph.h"

static inline void _add_edge_to_colours(hkey_t next_hkey,
                                        const Covg *covgs, Edges *edges,
                                        Edges new_edge,
                                        const dBGraph *db_graph)
{
  size_t col, ncols = db_graph->num_of_cols;

  if(db_graph->col_covgs != NULL) {
    for(col = 0; col < ncols; col++) {
      if(covgs[col] > 0 && db_node_covg(db_graph, next_hkey, col)) {
        edges[col] |= new_edge;
      }
    }
  }
  else {
    for(col = 0; col < ncols; col++) {
      if(covgs[col] > 0 && db_node_has_col(db_graph, next_hkey, col)) {
        edges[col] |= new_edge;
      }
    }
  }
}

// If two kmers are in a sample and the population has an edges between them,
// Add edge to sample

// Return 1 if changed; 0 otherwise
bool infer_pop_edges(const BinaryKmer node_bkey, Edges *edges,
                     const Covg *covgs, const dBGraph *db_graph)
{
  Edges uedges = 0, iedges = 0xf, add_edges, edge;
  size_t orient, nuc, col, kmer_size = db_graph->kmer_size;
  const size_t ncols = db_graph->num_of_cols;
  BinaryKmer bkey, bkmer;
  hkey_t next;
  Edges newedges[ncols];
  memcpy(newedges, edges, ncols * sizeof(Edges));

  // char tmp[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(node_bkey, db_graph->kmer_size, tmp);
  // status("Inferring %s", tmp);

  for(col = 0; col < ncols; col++) {
    uedges |= edges[col]; // union of edges
    iedges &= edges[col]; // intersection of edges
  }

  add_edges = uedges & ~iedges;

  if(!add_edges) return 0;

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                               : binary_kmer_right_shift_one_base(node_bkey));

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);
      if(add_edges & edge)
      {
        // get next bkmer, look up in graph
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);

        bkey = binary_kmer_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);
        ctx_assert(next != HASH_NOT_FOUND);

        _add_edge_to_colours(next, covgs, newedges, edge, db_graph);
      }
    }
  }

  int cmp = memcmp(edges, newedges, sizeof(Edges)*ncols);
  memcpy(edges, newedges, sizeof(Edges)*ncols);
  return (cmp != 0);
}

// Return 1 if changed; 0 otherwise
bool infer_all_edges(const BinaryKmer node_bkey, Edges *edges,
                     const Covg *covgs, const dBGraph *db_graph)
{
  Edges iedges = 0xff, edge;
  size_t orient, nuc, col, kmer_size = db_graph->kmer_size;
  const size_t ncols = db_graph->num_of_cols;
  BinaryKmer bkey, bkmer;
  hkey_t next;

  Edges newedges[ncols];
  memcpy(newedges, edges, ncols * sizeof(Edges));

  // intersection of edges
  for(col = 0; col < ncols; col++) iedges &= edges[col];

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                               : binary_kmer_right_shift_one_base(node_bkey));

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);
      if(!(iedges & edge))
      {
        // edges are missing from some samples
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);

        bkey = binary_kmer_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);

        if(next != HASH_NOT_FOUND) {
          _add_edge_to_colours(next, covgs, newedges, edge, db_graph);
        }
      }
    }
  }

  // Check if we changed the edges
  int cmp = memcmp(edges, newedges, ncols*sizeof(Edges));
  memcpy(edges, newedges, ncols*sizeof(Edges));
  return (cmp != 0);
}

static inline int infer_edges_node(hkey_t hkey,
                                   bool add_all_edges,
                                   Covg *tmp_covgs,
                                   const dBGraph *db_graph,
                                   size_t *num_nodes_modified)
{
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  Edges *edges = &db_node_edges(db_graph, hkey, 0);
  size_t col;

  // Create coverages that are zero or one depending on if node has colour
  if(db_graph->col_covgs == NULL) {
    for(col = 0; col < db_graph->num_of_cols; col++)
      tmp_covgs[col] = db_node_has_col(db_graph, hkey, col);
  } else {
    tmp_covgs = &db_node_covg(db_graph, hkey, 0);
  }

  (*num_nodes_modified)
    += (add_all_edges ? infer_all_edges(bkmer, edges, tmp_covgs, db_graph)
                      : infer_pop_edges(bkmer, edges, tmp_covgs, db_graph));

  return 0; // => keep iterating
}

typedef struct {
  const size_t nthreads;
  const bool add_all_edges;
  const dBGraph *db_graph;
  size_t num_nodes_modified;
} InferringEdges;

static void infer_edges_worker(void *arg, size_t threadid)
{
  InferringEdges *wrkr = (InferringEdges*)arg;
  size_t num_modified = 0;
  Covg covgs[wrkr->db_graph->num_of_cols];

  HASH_ITERATE_PART(&wrkr->db_graph->ht, threadid, wrkr->nthreads,
                    infer_edges_node,
                    wrkr->add_all_edges, covgs, wrkr->db_graph,
                    &num_modified);

  __sync_fetch_and_add((volatile size_t *)&wrkr->num_nodes_modified, num_modified);
}

size_t infer_edges(size_t nthreads, bool add_all_edges, const dBGraph *db_graph)
{
  ctx_assert(db_graph->node_in_cols != NULL || db_graph->col_covgs != NULL);
  ctx_assert(db_graph->col_edges != NULL);

  status("[inferedges] Processing stream");

  InferringEdges infedges = {.nthreads = nthreads,
                             .add_all_edges = add_all_edges,
                             .db_graph = db_graph,
                             .num_nodes_modified = 0};

  util_multi_thread(&infedges, nthreads, infer_edges_worker);

  return infedges.num_nodes_modified;
}
