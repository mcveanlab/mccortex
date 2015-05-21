#include "global.h"
#include "graph_file_reader.h"
#include "graph_format.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "range.h"

// Memory mapped files used in graph_files_merge()
#include <sys/mman.h>

void graph_header_alloc(GraphFileHeader *h, size_t num_of_cols)
{
  size_t i;

  if(num_of_cols > h->capacity) {
    h->ginfo = ctx_recallocarray(h->ginfo, h->capacity, num_of_cols, sizeof(GraphInfo));
    for(i = h->capacity; i < num_of_cols; i++)
      graph_info_alloc(&h->ginfo[i]);
    h->capacity = num_of_cols;
  }
}

void graph_header_dealloc(GraphFileHeader *h)
{
  size_t i;
  for(i = 0; i < h->capacity; i++)
    graph_info_dealloc(&h->ginfo[i]);
  ctx_free(h->ginfo);
  memset(h, 0, sizeof(*h));
}

void graph_header_print(const GraphFileHeader *header)
{
  printf("HEADER\n");
  printf("  version: %u\n", header->version);
  printf("  kmer_size: %u\n", header->kmer_size);
  printf("  num_of_bitfields: %u\n", header->num_of_bitfields);
  printf("  num_of_cols: %u\n", header->num_of_cols);
  printf("  [capacity: %zu]\n", header->capacity);
}

// Merge headers and set intersect name (if intersect_gname != NULL)
void graph_reader_merge_headers(GraphFileHeader *hdr,
                                const GraphFileReader *files, size_t num_files,
                                const char *intersect_gname)
{
  size_t i, j, ncols = 0, intocol, fromcol;

  for(i = 0; i < num_files; i++)
    ncols = MAX2(ncols, file_filter_into_ncols(&files[i].fltr));

  hdr->version = files[0].hdr.version;
  hdr->kmer_size = files[0].hdr.kmer_size;
  hdr->num_of_bitfields = files[0].hdr.num_of_bitfields;
  hdr->num_of_cols = ncols;
  graph_header_alloc(hdr, ncols);

  for(i = 0; i < num_files; i++) {
    for(j = 0; j < file_filter_num(&files[i].fltr); j++) {
      intocol = file_filter_intocol(&files[i].fltr, j);
      fromcol = file_filter_fromcol(&files[i].fltr, j);
      graph_info_merge(&hdr->ginfo[intocol], &files[i].hdr.ginfo[fromcol]);
    }
  }

  if(intersect_gname != NULL) {
    for(i = 0; i < ncols; i++) {
      if(graph_file_is_colour_loaded(i, files, num_files))
        graph_info_append_intersect(&hdr->ginfo[i].cleaning, intersect_gname);
    }
  }
}

// Return number of bytes read or die() with error
size_t graph_file_read_header(FILE *fh, GraphFileHeader *h, const char *path)
{
  size_t i;
  int bytes_read = 0;

  char magic_word[7];
  magic_word[6] = '\0';

  safe_fread(fh, magic_word, strlen("CORTEX"), "Magic word", path);
  if(strcmp(magic_word, "CORTEX") != 0) {
    die("Magic word doesn't match '%s' (start): %s", "CORTEX", path);
  }
  bytes_read += strlen("CORTEX");

  // Read version number, kmer_size, num bitfields, num colours
  safe_fread(fh, &h->version, sizeof(uint32_t), "graph version", path);
  safe_fread(fh, &h->kmer_size, sizeof(uint32_t), "kmer size", path);
  safe_fread(fh, &h->num_of_bitfields, sizeof(uint32_t), "num of bitfields", path);
  safe_fread(fh, &h->num_of_cols, sizeof(uint32_t), "number of colours", path);
  bytes_read += 4*sizeof(uint32_t);

  // Checks
  if(h->version > 7 || h->version < 4)
  {
    die("Sorry, we only support graph file versions 4, 5, 6 & 7 "
        "[version: %u; path: %s]\n", h->version, path);
  }

  if(h->kmer_size % 2 == 0)
  {
    die("kmer size is not an odd number [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3)
  {
    die("kmer size is less than three [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_bitfields * 32 < h->kmer_size)
  {
    die("Not enough bitfields for kmer size "
        "[kmer_size: %u; bitfields: %u; path: %s]\n",
        h->kmer_size, h->num_of_bitfields, path);
  }

  if((h->num_of_bitfields-1)*32 >= h->kmer_size) {
    die("using more than the minimum number of bitfields [path: %s]\n", path);
  }

  if(h->num_of_cols == 0)
    die("number of colours is zero [path: %s]\n", path);
  if(h->num_of_cols > 10000)
    die("Very high number of colours: %zu [path: %s]", (size_t)h->num_of_cols, path);

  // graph_header_alloc will only alloc or realloc if it needs to
  graph_header_alloc(h, h->num_of_cols);

  // Assume to be graph file now, any error is therefore fatal

  for(i = 0; i < h->num_of_cols; i++) {
    safe_fread(fh, &h->ginfo[i].mean_read_length, sizeof(uint32_t),
               "mean read length for each colour", path);
  }

  for(i = 0; i < h->num_of_cols; i++) {
    safe_fread(fh, &h->ginfo[i].total_sequence, sizeof(uint64_t),
               "total sequance loaded for each colour", path);
  }

  bytes_read += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint64_t));

  if(h->version >= 6)
  {
    for(i = 0; i < h->num_of_cols; i++)
    {
      uint32_t len;
      safe_fread(fh, &len, sizeof(uint32_t), "sample name length", path);
      if(len > 10000) die("Very big sample name. Length: %u", len);

      StrBuf *sbuf = &h->ginfo[i].sample_name;
      strbuf_ensure_capacity(sbuf, len);

      safe_fread(fh, sbuf->b, len, "sample name", path);
      bytes_read += sizeof(uint32_t) + len;

      sbuf->b[len] = '\0';

      size_t len2 = strlen(sbuf->b);
      sbuf->end = len2;
      sbuf->b[len2] = '\0';

      // Check sample length is as long as we were told
      if(len != len2)
      {
        warn("Sample %zu name has length %u but is only %zu chars long "
             "(premature '\\0') [path: %s]\n", i, len, len2, path);
      }
    }

    for(i = 0; i < h->num_of_cols; i++) {
      safe_fread(fh, &h->ginfo[i].seq_err, sizeof(long double),
                 "seq error rates", path);
    }

    bytes_read += sizeof(long double) * h->num_of_cols;

    for(i = 0; i < h->num_of_cols; i++)
    {
      ErrorCleaning *err_cleaning = &h->ginfo[i].cleaning;

      safe_fread(fh, &(err_cleaning->cleaned_tips),
                 sizeof(uint8_t), "tip cleaning", path);
      safe_fread(fh, &(err_cleaning->cleaned_snodes),
                 sizeof(uint8_t), "remove low covg supernodes", path);
      safe_fread(fh, &(err_cleaning->cleaned_kmers),
                 sizeof(uint8_t), "remove low covg kmers", path);
      safe_fread(fh, &(err_cleaning->is_graph_intersection),
                 sizeof(uint8_t), "cleaned against graph", path);

      uint32_t clean_snodes_thresh = 0, clean_kmers_thresh = 0;
      safe_fread(fh, &clean_snodes_thresh,
                 sizeof(uint32_t), "remove low covg supernode threshold", path);
      safe_fread(fh, &clean_kmers_thresh,
                 sizeof(uint32_t), "remove low covg kmer threshold", path);

      bytes_read += 4*sizeof(uint8_t) + 2*sizeof(uint32_t);

      // Fix for old versions with negative thresholds
      if(h->version <= 6) {
        if(!err_cleaning->cleaned_snodes && clean_snodes_thresh == (uint32_t)-1)
          clean_snodes_thresh = 0;
        if(!err_cleaning->cleaned_kmers && clean_kmers_thresh == (uint32_t)-1)
          clean_kmers_thresh = 0;
      }

      // Sanity checks
      if(!err_cleaning->cleaned_snodes && clean_snodes_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for supernodes "
             "when no cleaning was performed [path: %s]", path);
        clean_snodes_thresh = 0;
      }

      if(!err_cleaning->cleaned_kmers && clean_kmers_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for nodes "
             "when no cleaning was performed [path: %s]", path);
        clean_kmers_thresh = 0;
      }

      err_cleaning->clean_snodes_thresh = clean_snodes_thresh;
      err_cleaning->clean_kmers_thresh = clean_kmers_thresh;

      // Read cleaned against name
      uint32_t len;
      safe_fread(fh, &len, sizeof(uint32_t), "graph name length", path);
      if(len > 10000) die("Very big sample name. Length: %u", len);

      StrBuf *sbuf = &err_cleaning->intersection_name;
      strbuf_ensure_capacity(sbuf, len);

      safe_fread(fh, sbuf->b, len, "cleaned against graph name", path);
      sbuf->b[len] = '\0';

      bytes_read += sizeof(uint32_t) + len;

      size_t len2 = strlen(sbuf->b);
      sbuf->end = len2;
      sbuf->b[len2] = '\0';

      // Check sample length is as long as we were told
      if(len != len2)
      {
        warn("Sample %zu name has length %u but is only %zu chars long "
             "(premature '\\0') [path: %s]\n", i, len, len2, path);
      }
    }
  }

  // Read magic word at the end of header 'CORTEX'
  safe_fread(fh, magic_word, strlen("CORTEX"), "magic word (end)", path);
  if(strcmp(magic_word, "CORTEX") != 0)
  {
    die("Magic word doesn't match '%s' (end): '%s' [path: %s]\n",
        "CORTEX", magic_word, path);
  }
  bytes_read += strlen("CORTEX");

  return bytes_read;
}

// Only print errors once - these are externally visible
bool greader_zero_covg_error = false, greader_missing_covg_error = false;

size_t graph_file_read_kmer(FILE *fh, const GraphFileHeader *h, const char *path,
                            BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  size_t i, num_bytes_read;
  char kstr[MAX_KMER_SIZE+1];

  num_bytes_read = fread(bkmer->b, 1, sizeof(BinaryKmer), fh);

  if(num_bytes_read == 0) return 0;
  if(num_bytes_read != sizeof(uint64_t)*h->num_of_bitfields)
    die("Unexpected end of file: %s", path);

  safe_fread(fh, covgs, h->num_of_cols * sizeof(uint32_t), "Coverages", path);
  safe_fread(fh, edges, h->num_of_cols * sizeof(uint8_t), "Edges", path);
  num_bytes_read += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint8_t));

  // Check top word of each kmer
  if(binary_kmer_oversized(*bkmer, h->kmer_size))
    die("Oversized kmer in path [kmer: %u]: %s", h->kmer_size, path);

  // Check covg is not 0 for all colours
  for(i = 0; i < h->num_of_cols && covgs[i] == 0; i++) {}
  if(i == h->num_of_cols && !greader_zero_covg_error) {
    binary_kmer_to_str(*bkmer, h->kmer_size, kstr);
    warn("Kmer has zero covg in all colours [kmer: %s; path: %s]", kstr, path);
    greader_zero_covg_error = true;
  }

  // Check edges => coverage
  for(i = 0; i < h->num_of_cols && (!edges[i] || covgs[i]); i++) {}
  if(i < h->num_of_cols && !greader_missing_covg_error) {
    binary_kmer_to_str(*bkmer, h->kmer_size, kstr);
    warn("Kmer has edges but no coverage [kmer: %s; path: %s]", kstr, path);
    greader_missing_covg_error = true;
  }

  return num_bytes_read;
}

// Print some output
static void graph_loading_print_status(const GraphFileReader *file)
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
  size_t i, ncols_used = file_filter_into_ncols(&file->fltr), fromcol, intocol;
  FileFilter *fltr = &file->fltr;
  GraphFileHeader *hdr = &file->hdr;

  ctx_assert(file_filter_num(fltr) > 0);

  // Print status
  graph_loading_print_status(file);

  if(!file_filter_isstdin(fltr) && fseek(file->fh, file->hdr_size, SEEK_SET) != 0)
    die("fseek failed: %s", strerror(errno));

  // Check we can load this graph file into db_graph (kmer size + num colours)
  if(hdr->kmer_size != graph->kmer_size)
  {
    die("Graph has different kmer size [kmer_size: %u vs %zu; path: %s]",
        hdr->kmer_size, graph->kmer_size, fltr->path.b);
  }

  if(ncols_used > graph->num_of_cols)
  {
    die("Program has not assigned enough colours! "
        "[colours in graph: %zu vs file: %zu; path: %s]",
        graph->num_of_cols, ncols_used, fltr->path.b);
  }

  for(i = 0; i < file_filter_num(fltr); i++) {
    fromcol = file_filter_fromcol(fltr, i);
    intocol = file_filter_intocol(fltr, i);
    graph_info_merge(ginfo+intocol, hdr->ginfo+fromcol);
  }

  // Update number of colours loaded
  graph->num_of_cols_used = MAX2(graph->num_of_cols_used, ncols_used);

  // Read kmers, align colours to those they are updating
  //  e.g. covgs[i] -> colour i in the graph
  BinaryKmer bkmer;
  Covg covgs[ncols_used];
  Edges edges[ncols_used];

  size_t nkmers_parsed, num_of_kmers_loaded = 0;
  uint64_t num_of_kmers_already_loaded = graph->ht.num_kmers;

  for(nkmers_parsed = 0;
      graph_file_read_reset(file, ncols_used, &bkmer, covgs, edges);
      nkmers_parsed++)
  {
    // If kmer has no covg or edges -> don't load
    Covg keep_kmer = 0;
    for(i = 0; i < ncols_used; i++) keep_kmer |= covgs[i] | edges[i];
    if(keep_kmer == 0) continue;

    if(prefs.boolean_covgs)
      for(i = 0; i < ncols_used; i++)
        covgs[i] = covgs[i] > 0;

    // Fetch node in the de bruijn graph
    hkey_t hkey;

    if(prefs.must_exist_in_graph)
    {
      hkey = hash_table_find(&graph->ht, bkmer);
      if(hkey == HASH_NOT_FOUND) continue;

      // Edges union_edges = db_node_get_edges_union(graph, hkey);
      Edges union_edges = prefs.must_exist_in_edges[hkey];

      for(i = 0; i < ncols_used; i++) edges[i] &= union_edges;
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
      for(i = 0; i < ncols_used; i++) {
        db_node_or_col(graph, hkey, i, (covgs[i] || edges[i]));
      }
    }

    if(graph->col_covgs != NULL) {
      for(i = 0; i < ncols_used; i++)
        db_node_add_col_covg(graph, hkey, i, covgs[i]);
    }

    // Merge all edges into one colour
    if(graph->col_edges != NULL)
    {
      Edges *col_edges = &db_node_edges(graph, hkey, 0);

      if(graph->num_edge_cols == 1) {
        for(i = 0; i < ncols_used; i++)
          col_edges[0] |= edges[i];
      }
      else {
        for(i = 0; i < ncols_used; i++)
          col_edges[i] |= edges[i];
      }
    }

    num_of_kmers_loaded++;
  }

  if(file->num_of_kmers >= 0 && nkmers_parsed != (uint64_t)file->num_of_kmers)
  {
    warn("More kmers in graph than expected [expected: %zu; actual: %zu; "
         "path: %s]", (size_t)file->num_of_kmers, nkmers_parsed, fltr->path.b);
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

//
// Merging, filtering, combining graph files
//

// Load a kmer and write to a file one kmer at a time
// Optionally filter against the graph currently loaded
//   (i.e. only keep nodes and edges that are in the graph)
// Same functionality as graph_files_merge, but faster if dealing with only one
// input file. Reads in and dumps one kmer at a time
// parameters:
//   `only_load_if_in_edges`: Edges to mask edges with, 1 per hash table entry
size_t graph_stream_filter(const char *out_ctx_path, const GraphFileReader *file,
                           const dBGraph *db_graph, const GraphFileHeader *hdr,
                           const Edges *only_load_if_in_edges)
{
  const FileFilter *fltr = &file->fltr;
  bool only_load_if_in_graph = (only_load_if_in_edges != NULL);
  status("Filtering %s to %s with stream filter", fltr->path.b,
         futil_outpath_str(out_ctx_path));

  FILE *out = futil_fopen(out_ctx_path, "w");

  graph_loading_print_status(file);

  size_t i, nodes_dumped = 0, ncols = file_filter_into_ncols(fltr);

  graph_write_header(out, hdr);

  BinaryKmer bkmer;
  Covg covgs[ncols];
  Edges edges[ncols];

  while(graph_file_read_reset(file, ncols, &bkmer, covgs, edges))
  {
    // Collapse down colours
    Covg keep_kmer = 0;
    for(i = 0; i < ncols; i++) keep_kmer |= covgs[i] | edges[i];

    // If kmer has no covg or edges -> don't load
    if(keep_kmer)
    {
      if(only_load_if_in_graph)
      {
        hkey_t node = hash_table_find(&db_graph->ht, bkmer);

        if(node != HASH_NOT_FOUND) {
          Edges union_edges = only_load_if_in_edges[node];
          for(i = 0; i < ncols; i++) edges[i] &= union_edges;
        }
        else keep_kmer = 0;
      }

      if(keep_kmer) {
        graph_write_kmer(out, hdr->num_of_bitfields, hdr->num_of_cols,
                         bkmer, covgs, edges);
        nodes_dumped++;
      }
    }
  }

  fflush(out);
  fclose(out);

  graph_writer_print_status(nodes_dumped, hdr->num_of_cols,
                     out_ctx_path, CTX_GRAPH_FILEFORMAT);

  return nodes_dumped;
}

size_t graph_stream_filter_mkhdr(const char *out_ctx_path,
                                 GraphFileReader *file,
                                 const dBGraph *db_graph,
                                 const Edges *only_load_if_in_edges,
                                 const char *intersect_gname)
{
  ctx_assert(intersect_gname == NULL || db_graph->col_edges != NULL);
  ctx_assert(intersect_gname == NULL || only_load_if_in_edges != NULL);

  FileFilter *fltr = &file->fltr;
  size_t i, nodes_dumped, ncols = file_filter_into_ncols(fltr);

  GraphFileHeader outheader;
  memset(&outheader, 0, sizeof(outheader));

  outheader.version = CTX_GRAPH_FILEFORMAT;
  outheader.kmer_size = db_graph->kmer_size;
  outheader.num_of_bitfields = (db_graph->kmer_size*2+63)/64;
  outheader.num_of_cols = ncols;
  graph_header_alloc(&outheader, outheader.num_of_cols);

  uint32_t fromcol, intocol;
  for(i = 0; i < file_filter_num(fltr); i++) {
    fromcol = file_filter_fromcol(fltr, i);
    intocol = file_filter_intocol(fltr, i);
    GraphInfo *ginfo = &outheader.ginfo[intocol];
    graph_info_merge(ginfo, file->hdr.ginfo + fromcol);
    if(intersect_gname != NULL)
      graph_info_append_intersect(&ginfo->cleaning, intersect_gname);
  }

  nodes_dumped = graph_stream_filter(out_ctx_path, file,
                                     db_graph, &outheader,
                                     only_load_if_in_edges);
  graph_header_dealloc(&outheader);

  return nodes_dumped;
}

// `kmers_loaded`: means all kmers to dump have been loaded
// `colours_loaded`: means all kmer data have been loaded
// `only_load_if_in_edges`: Edges to mask edges with, 1 per hash table entry
//   should already be set
size_t graph_files_merge(const char *out_ctx_path,
                         GraphFileReader *files, size_t num_files,
                         bool kmers_loaded, bool colours_loaded,
                         const Edges *only_load_if_in_edges,
                         GraphFileHeader *hdr, dBGraph *db_graph)
{
  bool only_load_if_in_graph = (only_load_if_in_edges != NULL);
  ctx_assert(!only_load_if_in_graph || kmers_loaded);
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  ctx_assert(!colours_loaded || kmers_loaded);
  ctx_assert(hdr != NULL);

  size_t i, j, f, ncols, output_colours = 0;

  for(i = 0; i < num_files; i++) {
    ncols = file_filter_into_ncols(&files[i].fltr);
    output_colours = MAX2(output_colours, ncols);

    if(files[i].hdr.kmer_size != files[0].hdr.kmer_size) {
      die("Kmer-size mismatch %u vs %u [%s vs %s]",
          files[0].hdr.kmer_size, files[i].hdr.kmer_size,
          files[0].fltr.path.b, files[i].fltr.path.b);
    }
  }

  ctx_assert(output_colours <= hdr->num_of_cols);

  if(kmers_loaded && colours_loaded)
  {
    return graph_file_save(out_ctx_path, db_graph, hdr,
                           0, NULL, 0, output_colours);
  }
  else if(num_files == 1)
  {
    return graph_stream_filter(out_ctx_path, &files[0], db_graph, hdr,
                               only_load_if_in_edges);
  }

  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs prefs
    = {.db_graph = db_graph,
       .boolean_covgs = false,
       .must_exist_in_graph = only_load_if_in_graph,
       .must_exist_in_edges = only_load_if_in_edges,
       .empty_colours = false};

  if(output_colours <= db_graph->num_of_cols)
  {
    // Can load all files at once
    status("Loading and saving %zu colours at once", output_colours);

    if(!kmers_loaded) {
      for(i = 0; i < num_files; i++)
        graph_load(&files[i], prefs, &stats);
    }

    hash_table_print_stats(&db_graph->ht);
    graph_file_save(out_ctx_path, db_graph, hdr, 0, NULL, 0, output_colours);
  }
  else
  {
    ctx_assert2(strcmp(out_ctx_path,"-") != 0,
                "Cannot use STDOUT for output if not enough colours to load");

    // Have to load a few colours at a time then dump, rinse and repeat
    status("[mmap] Saving %zu colours, %zu colours at a time",
           output_colours, db_graph->num_of_cols);

    // Open file, write header
    FILE *fout = futil_fopen(out_ctx_path, "r+");

    size_t hdr_size = graph_write_header(fout, hdr);

    // Load all kmers into flat graph
    if(!kmers_loaded)
      graph_files_load_flat(files, num_files, prefs, &stats);

    // print file outline
    status("Generated merged hash table\n");
    hash_table_print_stats(&db_graph->ht);

    // Write empty file
    size_t file_len = hdr_size;
    file_len += graph_write_empty(db_graph, fout, output_colours);
    fflush(fout);

    // Open memory mapped file
    void *mmap_ptr = mmap(NULL, file_len, PROT_WRITE, MAP_SHARED, fileno(fout), 0);

    if(mmap_ptr == MAP_FAILED)
      die("Cannot memory map file: %s [%s]", out_ctx_path, strerror(errno));

    size_t num_kmer_cols = db_graph->ht.capacity * db_graph->num_of_cols;
    size_t firstcol, lastcol, fromcol, intocol;

    FileFilter origfltr;
    memset(&origfltr, 0, sizeof(origfltr));
    bool files_loaded = false;

    for(firstcol = 0; firstcol < output_colours; firstcol += db_graph->num_of_cols)
    {
      lastcol = MIN2(firstcol + db_graph->num_of_cols - 1, output_colours-1);

      // Wipe colour coverages and edges if needed
      if(firstcol == 0 || files_loaded) {
        status("Wiping colours");
        memset(db_graph->col_edges, 0, num_kmer_cols * sizeof(Edges));
        memset(db_graph->col_covgs, 0, num_kmer_cols * sizeof(Covg));
      }

      files_loaded = false;

      // Loop over files, loading only the colours currently required
      for(f = 0; f < num_files; f++)
      {
        FileFilter *fltr = &files[f].fltr;
        file_filter_copy(&origfltr, fltr);

        // Update filter to only load useful colours
        for(i = j = 0; i < file_filter_num(fltr); i++) {
          fromcol = file_filter_fromcol(fltr, i);
          intocol = file_filter_intocol(fltr, i);
          if(firstcol <= intocol && intocol <= lastcol) {
            fltr->filter.b[j++] = (Filter){.from = fromcol,
                                              .into = intocol-firstcol};
          }
        }
        file_filter_num(fltr) = j;

        if(file_filter_num(fltr) > 0)
          files_loaded |= (graph_load(&files[f], prefs, &stats) > 0);

        // Restore original filter
        file_filter_copy(fltr, &origfltr);
      }

      // if files_loaded, dump
      if(files_loaded) {
        if(db_graph->num_of_cols == 1)
          status("Dumping into colour %zu...\n", firstcol);
        else
          status("Dumping into colours %zu-%zu...\n", firstcol, lastcol);

        graph_update_mmap_kmers(db_graph, 0, lastcol-firstcol+1,
                                firstcol, output_colours,
                                mmap_ptr, hdr_size);
      }
    }

    if(munmap(mmap_ptr, file_len) == -1)
      die("Cannot release mmap file: %s [%s]", out_ctx_path, strerror(errno));

    fclose(fout);
    file_filter_close(&origfltr);

    // Force update of file timestamp so mmap doesn't mess with it
    futil_update_timestamp(out_ctx_path);

    // Print output status
    graph_writer_print_status(db_graph->ht.num_kmers, output_colours,
                              out_ctx_path, CTX_GRAPH_FILEFORMAT);
  }

  return db_graph->ht.num_kmers;
}

// if intersect_gname != NULL: only load kmers that are already in the hash table
//    and use string as name for cleaning against
// returns the number of kmers written
size_t graph_files_merge_mkhdr(const char *out_ctx_path,
                               GraphFileReader *files, size_t num_files,
                               bool kmers_loaded, bool colours_loaded,
                               const Edges *only_load_if_in_edges,
                               const char *intersect_gname, dBGraph *db_graph)
{
  size_t num_kmers;
  GraphFileHeader gheader;
  memset(&gheader, 0, sizeof(gheader));
  gheader.version = CTX_GRAPH_FILEFORMAT;
  gheader.kmer_size = db_graph->kmer_size;
  gheader.num_of_bitfields = (db_graph->kmer_size*2+63)/64;

  graph_reader_merge_headers(&gheader, files, num_files, intersect_gname);

  num_kmers = graph_files_merge(out_ctx_path, files, num_files,
                                kmers_loaded, colours_loaded,
                                only_load_if_in_edges,
                                &gheader, db_graph);

  graph_header_dealloc(&gheader);
  return num_kmers;
}

// Load all files into colour 0
void graph_files_load_flat(GraphFileReader *gfiles, size_t num_files,
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

