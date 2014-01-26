#ifndef GRAPH_PATH_H_
#define GRAPH_PATH_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"

//
// Functions on graph+paths
//

void graph_path_check_valid(const dBGraph *db_graph, dBNode node, size_t col,
                            const Nucleotide *bases, size_t nbases);

void graph_paths_check_all_paths(const dBGraph *db_graph);

#endif /* GRAPH_PATH_H_ */
