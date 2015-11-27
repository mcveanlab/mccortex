#include "global.h"
#include "graphs_load.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"

//
// Graph loading stats
//

void graph_loading_stats_destroy(GraphLoadingStats *s)
{
  ctx_free(s->nkmers);
  ctx_free(s->sumcov);
  memset(s, 0, sizeof(*s));
}

void graph_loading_stats_capacity(GraphLoadingStats *s, size_t n)
{
  if(n > s->ncols) {
    s->nkmers = ctx_recallocarray(s->nkmers, s->ncols, n, sizeof(s->nkmers[0]));
    s->sumcov = ctx_recallocarray(s->sumcov, s->ncols, n, sizeof(s->sumcov[0]));
    s->ncols = n;
  }
}


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
                  GraphLoadingStats *stats)
{
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
  hkey_t hkey;
  size_t nkmers_read = 0, nkmers_loaded = 0, nkmers_novel = 0;

  if(stats) graph_loading_stats_capacity(stats, ncols);

  for(; graph_file_read_reset(file, &bkmer, covgs, edges); nkmers_read++)
  {
    // If kmer has no covg -> don't load
    Covg keep_kmer = 0;
    for(i = 0; i < ncols; i++) keep_kmer |= covgs[i];
    if(keep_kmer == 0) continue;

    if(stats) {
      for(i = 0; i < ncols; i++) {
        stats->nkmers[i] += covgs[i] > 0;
        stats->sumcov[i] += covgs[i];
      }
    }

    if(prefs.boolean_covgs)
      for(i = 0; i < ncols; i++)
        covgs[i] = covgs[i] > 0;

    // Fetch node in the de bruijn graph
    if(prefs.must_exist_in_graph)
    {
      if((hkey = hash_table_find(&graph->ht, bkmer)) == HASH_NOT_FOUND) continue;
    }
    else
    {
      bool found;
      hkey = hash_table_find_or_insert(&graph->ht, bkmer, &found);
      if(prefs.empty_colours && found) die("Duplicate kmer loaded");
      nkmers_novel += !found;
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
      // Edges edge_mask = db_node_get_edges_union(graph, hkey);
      Edges edge_mask = 0xff;

      if(prefs.must_exist_in_edges)
        edge_mask = prefs.must_exist_in_edges[hkey];
      else if(prefs.must_exist_in_graph)
        edge_mask = db_node_get_edges_union(graph, hkey);

      Edges *col_edges = &db_node_edges(graph, hkey, 0);

      if(graph->num_edge_cols == 1) {
        for(i = 0; i < ncols; i++)
          col_edges[0] |= edges[i] & edge_mask;
      }
      else {
        for(i = 0; i < ncols; i++)
          col_edges[i] |= edges[i] & edge_mask;
      }
    }

    nkmers_loaded++;
  }

  if(file->num_of_kmers >= 0 && nkmers_read != (uint64_t)file->num_of_kmers)
  {
    warn("%s kmers in the graph file than expected "
         "[exp: %zu; act: %zu; path: %s]",
         nkmers_read > (uint64_t)file->num_of_kmers ? "More" : "Fewer",
         (size_t)file->num_of_kmers, nkmers_read, fltr->path.b);
  }

  if(stats != NULL)
  {
    stats->nkmers_read += nkmers_read;
    stats->nkmers_loaded += nkmers_loaded;
    stats->nkmers_novel += nkmers_novel;
    // for(i = 0; i < file_filter_num(fltr); i++) {
    //   fromcol = file_filter_fromcol(fltr,i);
    //   stats->total_bases_read += hdr->ginfo[fromcol].total_sequence;
    // }
  }

  char n0[50], n1[50];
  status("[GReader] Loaded %s / %s (%.2f%%) of kmers parsed",
         ulong_to_str(nkmers_loaded, n0),
         ulong_to_str(nkmers_read, n1),
         safe_percent(nkmers_loaded, nkmers_read));

  return nkmers_loaded;
}

// Load all files into colour 0
void graphs_load_files_flat(GraphFileReader *gfiles, size_t num_files,
                            GraphLoadingPrefs prefs, GraphLoadingStats *stats)
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

