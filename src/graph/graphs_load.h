#ifndef GRAPHS_LOAD_H_
#define GRAPHS_LOAD_H_

#include "db_graph.h"
#include "graph_file_reader.h"

//
// Load graphs into the de Bruijn graph in memory from disk
//

// Stucture for specifying how to load data
typedef struct
{
  dBGraph *db_graph;
  bool boolean_covgs; // Update covg by at most 1
  bool must_exist_in_graph;
  const Edges *must_exist_in_edges;
  // if empty_colours is true an error is thrown if a kmer from a graph file
  // is already in the graph
  bool empty_colours;
} GraphLoadingPrefs;

typedef struct
{
  uint64_t nkmers_read, nkmers_loaded, nkmers_novel;
  uint64_t *nkmers, *sumcov; // one per colour
  size_t ncols; // max colour + 1 that we loaded into
} GraphLoadingStats;

void graph_loading_stats_capacity(GraphLoadingStats *s, size_t n);
void graph_loading_stats_destroy(GraphLoadingStats *stats);

static inline GraphLoadingPrefs graph_loading_prefs(dBGraph *graph)
{
  GraphLoadingPrefs prefs =
  {
    .db_graph = graph,
    .boolean_covgs = false,
    .must_exist_in_graph = false,
    .must_exist_in_edges = NULL,
    .empty_colours = false
  };
  return prefs;
}

// Print loading message
void graph_loading_print_status(const GraphFileReader *file);

// if only_load_if_in_colour is >= 0, only kmers with coverage in existing
// colour only_load_if_in_colour will be loaded.
// if clean_colours != 0 an error is thrown if a node already exists
// returns the number of colours in the binary
// If stats != NULL, updates:
//   stats->num_kmers_loaded
//   stats->total_bases_read
// If header is != NULL, header will be stored there.  Be sure to free.
size_t graph_load(GraphFileReader *file, const GraphLoadingPrefs prefs,
                  GraphLoadingStats *stats);

// Load all files into colour 0
void graphs_load_files_flat(GraphFileReader *gfiles, size_t num_files,
                           GraphLoadingPrefs prefs, GraphLoadingStats *stats);

#endif /* GRAPHS_LOAD_H_ */
