#include "global.h"
#include "graphs_load.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "range.h"

// Memory mapped files used in graph_writer_merge()
#include <sys/mman.h>

// Print loading message
void graph_loading_print_status(const GraphFileReader *file)
{
  const FileFilter *fltr = &file->fltr;
  char nkmers_str[100], filesize_str[100];

  file_filter_status(fltr);

  if(isatty(fileno(file->fh))) status("  reading from a stream.");
  else {
    ulong_to_str(file->num_of_kmers, nkmers_str);
    bytes_to_str(file->file_size, 1, filesize_str);
    status("[GReader] %s kmers, %s filesize", nkmers_str, filesize_str);
  }
}

// if only_load_if_in_colour is >= 0, only kmers with coverage in existing
// colour only_load_if_in_colour will be loaded.
// We assume only_load_if_in_colour < load_first_colour_into
// if all_kmers_are_unique != 0 an error is thrown if a node already exists
// If stats != NULL, updates:
//   stats->num_kmers_loaded
//   stats->total_bases_read

/*!
  @return number of kmers loaded
 */
size_t graph_load(GraphFileReader *file, const GraphLoadingPrefs prefs,
                  LoadingStats *stats)
{
  ctx_assert(!prefs.must_exist_in_graph || prefs.db_graph->col_edges == NULL ||
             prefs.must_exist_in_edges != NULL);
  ctx_assert(!prefs.must_exist_in_graph || !prefs.empty_colours);

  dBGraph *graph = prefs.db_graph;
  GraphInfo *ginfo = graph->ginfo;
  size_t i, ncols = file_filter_into_ncols(&file->fltr), fromcol, intocol;
  FileFilter *fltr = &file->fltr;
  GraphFileHeader *hdr = &file->hdr;

  ctx_assert(file_filter_num(fltr) > 0);

  // Print status
  graph_loading_print_status(file);

  // Functions such as merging multiple coloured graphs required us to load
  // each graph twice. It is convenient to do the fseek here.
  if(!file_filter_isstdin(fltr)) {
    if(graph_file_fseek(file, file->hdr_size, SEEK_SET) != 0)
      die("fseek failed: %s", strerror(errno));
  }

  // Check we can load this graph file into db_graph (kmer size + num colours)
  if(hdr->kmer_size != graph->kmer_size)
  {
    die("Graph has different kmer size [kmer_size: %u vs %zu; path: %s]",
        hdr->kmer_size, graph->kmer_size, fltr->path.b);
  }

  if(ncols > graph->num_of_cols)
  {
    die("Program has not assigned enough colours! "
        "[colours in graph: %zu vs file: %zu; path: %s]",
        graph->num_of_cols, ncols, fltr->path.b);
  }

  for(i = 0; i < file_filter_num(fltr); i++) {
    fromcol = file_filter_fromcol(fltr, i);
    intocol = file_filter_intocol(fltr, i);
    graph_info_merge(ginfo+intocol, hdr->ginfo+fromcol);
  }

  // Update number of colours loaded
  graph->num_of_cols_used = MAX2(graph->num_of_cols_used, ncols);

  // Read kmers, align colours to those they are updating
  //  e.g. covgs[i] -> colour i in the graph
  BinaryKmer bkmer;
  Covg covgs[ncols];
  Edges edges[ncols];

  size_t nkmers_parsed, num_of_kmers_loaded = 0;
  uint64_t num_of_kmers_already_loaded = graph->ht.num_kmers;

  for(nkmers_parsed = 0;
      graph_file_read_reset(file, &bkmer, covgs, edges);
      nkmers_parsed++)
  {
    // If kmer has no covg or edges -> don't load
    Covg keep_kmer = 0;
    for(i = 0; i < ncols; i++) keep_kmer |= covgs[i] | edges[i];
    if(keep_kmer == 0) continue;

    if(prefs.boolean_covgs)
      for(i = 0; i < ncols; i++)
        covgs[i] = covgs[i] > 0;

    // Fetch node in the de bruijn graph
    hkey_t hkey;

    if(prefs.must_exist_in_graph)
    {
      hkey = hash_table_find(&graph->ht, bkmer);
      if(hkey == HASH_NOT_FOUND) continue;

      // Edges union_edges = db_node_get_edges_union(graph, hkey);
      Edges union_edges = prefs.must_exist_in_edges[hkey];

      for(i = 0; i < ncols; i++) edges[i] &= union_edges;
    }
    else
    {
      bool found;
      hkey = hash_table_find_or_insert(&graph->ht, bkmer, &found);

      if(prefs.empty_colours && found)
        die("Duplicate kmer loaded");
    }

    // Set presence in colours
    if(graph->node_in_cols != NULL) {
      for(i = 0; i < ncols; i++) {
        db_node_or_col(graph, hkey, i, (covgs[i] || edges[i]));
      }
    }

    if(graph->col_covgs != NULL) {
      for(i = 0; i < ncols; i++)
        db_node_add_col_covg(graph, hkey, i, covgs[i]);
    }

    // Merge all edges into one colour
    if(graph->col_edges != NULL)
    {
      Edges *col_edges = &db_node_edges(graph, hkey, 0);

      if(graph->num_edge_cols == 1) {
        for(i = 0; i < ncols; i++)
          col_edges[0] |= edges[i];
      }
      else {
        for(i = 0; i < ncols; i++)
          col_edges[i] |= edges[i];
      }
    }

    num_of_kmers_loaded++;
  }

  if(file->num_of_kmers >= 0 && nkmers_parsed != (uint64_t)file->num_of_kmers)
  {
    warn("%s kmers in the graph file than expected "
         "[exp: %zu; act: %zu; path: %s]",
         nkmers_parsed > (uint64_t)file->num_of_kmers ? "More" : "Fewer",
         (size_t)file->num_of_kmers, nkmers_parsed, fltr->path.b);
  }

  if(stats != NULL)
  {
    stats->num_kmers_loaded += num_of_kmers_loaded;
    stats->num_kmers_novel += graph->ht.num_kmers - num_of_kmers_already_loaded;
    for(i = 0; i < file_filter_num(fltr); i++) {
      fromcol = file_filter_fromcol(fltr,i);
      stats->total_bases_read += hdr->ginfo[fromcol].total_sequence;
    }
  }

  char parsed_nkmers_str[100], loaded_nkmers_str[100];
  double loaded_nkmers_pct = 0;
  if(nkmers_parsed)
    loaded_nkmers_pct = (100.0 * num_of_kmers_loaded) / nkmers_parsed;

  ulong_to_str(num_of_kmers_loaded, loaded_nkmers_str);
  ulong_to_str(nkmers_parsed, parsed_nkmers_str);
  status("[GReader] Loaded %s / %s (%.2f%%) of kmers parsed",
         loaded_nkmers_str, parsed_nkmers_str, loaded_nkmers_pct);

  return num_of_kmers_loaded;
}

// Load all files into colour 0
void graphs_load_files_flat(GraphFileReader *gfiles, size_t num_files,
                            GraphLoadingPrefs prefs, LoadingStats *stats)
{
  size_t i;
  FileFilter origfltr;
  memset(&origfltr, 0, sizeof(origfltr));

  for(i = 0; i < num_files; i++)
  {
    file_filter_copy(&origfltr, &gfiles[i].fltr);
    file_filter_flatten(&gfiles[i].fltr, 0);
    graph_load(&gfiles[i], prefs, stats);
    file_filter_copy(&gfiles[i].fltr, &origfltr);
  }

  file_filter_close(&origfltr);
}

