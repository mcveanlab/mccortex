#ifndef GRAPH_SEARCH_H_
#define GRAPH_SEARCH_H_

#include "cortex_types.h"
#include "binary_kmer.h"
#include "graph_file_reader.h"

//
// Search a sorted graph file on disk
//

typedef struct GraphFileSearch GraphFileSearch;

GraphFileSearch *graph_search_new(GraphFileReader *file);
void graph_search_destroy(GraphFileSearch *gs);

bool graph_search_find(GraphFileSearch *gs, BinaryKmer bkey,
                       Covg *covgs, Edges *edges);

void graph_search_fetch(GraphFileSearch *gs, size_t idx,
                        BinaryKmer *bkey, Covg *covgs, Edges *edges);

void graph_search_rand(GraphFileSearch *gs,
                       BinaryKmer *bkey, Covg *covgs, Edges *edges);

#endif /* GRAPH_SEARCH_H_ */
