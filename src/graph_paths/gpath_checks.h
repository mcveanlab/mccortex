#ifndef GPATH_CHECKS_H_
#define GPATH_CHECKS_H_

#include "gpath.h"
#include "db_graph.h"
#include "graph_file_reader.h"
#include "gpath_reader.h"

// Similar to path_file_reader.c:path_file_load_check()
// Check kmer size matches and sample names match
void graphs_gpaths_compatible(const GraphFileReader *graphs, size_t num_graphs,
                              const GPathReader *gpaths, size_t num_gpaths);

// 1) check dBNode following `node` has indegree >1 in sample ctxcol
// 2) follow path, check each junction matches up with a node with outdegree >1
// col is graph colour
bool gpath_checks_path_col(dBNode node, const GPath *gpath, int exp_klen,
                           size_t ctxcol, const dBGraph *db_graph);

// Returns false on first error
bool gpath_checks_path(hkey_t hkey, const GPath *gpath, int exp_klen,
                       const dBGraph *db_graph);

// Returns false on error
bool gpath_checks_all_paths(const dBGraph *db_graph, size_t nthreads);

// Dies on first error
void gpath_checks_counts(const dBGraph *db_graph);

#endif /* GPATH_CHECKS_H_ */