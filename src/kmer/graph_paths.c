#include "global.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_paths.h"

// 1) check node after node has indegree >1 in sample ctxcol
// 2) follow path, check each junction matches up with a node with outdegree >1
void graph_path_check_valid(const dBGraph *db_graph, dBNode node, size_t col,
                            const Nucleotide *bases, size_t nbases)
{
  assert(db_graph->num_edge_cols == db_graph->num_of_cols ||
         db_graph->node_in_cols != NULL);

  BinaryKmer bkmer;
  Edges edges;
  hkey_t nodes[4];
  Orientation orients[4];
  Nucleotide nucs[4];
  size_t i, j, n, edgecol = db_graph->num_edge_cols > 1 ? col : 0;
  // length is kmers and juctions
  size_t klen, plen;

  // Check node is in this colour
  if(db_graph->node_in_cols != NULL) {
    assert(db_node_has_col(db_graph, col, node.key));
  } else if(db_graph->col_covgs != NULL) {
    assert(db_node_covg(db_graph, col, node.key) > 0);
  }

  for(klen = 0, plen = 0; plen < nbases; klen++)
  {
    bkmer = db_node_bkmer(db_graph, node.key);
    edges = db_node_edges(db_graph, edgecol, node.key);

    // char bkmerstr[MAX_KMER_SIZE+1];
    // binary_kmer_to_str(bkmer, db_graph->kmer_size, bkmerstr);
    // status("klen: %zu plen: %zu %zu:%i %s",
    //        klen, plen, (size_t)node.key, node.orient, bkmerstr);

    if(klen == 1) {
      dBNode rnode = db_node_reverse(node);
      Edges backedges = db_node_oriented_edges_in_col(rnode, col, db_graph);
      assert(edges_get_outdegree(backedges, FORWARD) > 1);
    }

    n = db_graph_next_nodes(db_graph, bkmer, node.orient,
                            edges, nodes, orients, nucs);

    // Reduce to nodes in our colour if edges limited
    if(db_graph->num_edge_cols == 1 && db_graph->node_in_cols != NULL) {
      for(i = 0, j = 0; i < n; i++) {
        if(db_node_has_col(db_graph, col, nodes[i])) {
          nodes[j] = nodes[i];
          orients[j] = orients[i];
          nucs[j] = nucs[i];
          j++;
        }
      }
      n = j; // update number of next nodes
    }

    assert(n > 0);

    // If fork check nucleotide
    if(n > 1) {
      assert(bases[plen] < 4);
      for(i = 0; i < n && nucs[i] != bases[plen]; i++);
      assert(i < n && nucs[i] == bases[plen]);
      node.key = nodes[i];
      node.orient = orients[i];
      plen++;
    }
    else {
      node.key = nodes[0];
      node.orient = orients[0];
    }
  }
}
