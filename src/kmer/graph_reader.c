#include "global.h"
#include "graph_file_filter.h"
#include "graph_format.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "range.h"

void graph_header_alloc(GraphFileHeader *h, size_t num_of_cols)
{
  size_t i, old_cap = h->capacity;

  if(h->capacity == 0)
    h->ginfo = calloc2(num_of_cols, sizeof(GraphInfo));
  else if(num_of_cols > h->capacity)
    h->ginfo = realloc2(h->ginfo, num_of_cols * sizeof(GraphInfo));

  for(i = old_cap; i < num_of_cols; i++)
    graph_info_alloc(h->ginfo + i);

  h->capacity = MAX2(old_cap, num_of_cols);
}

void graph_header_dealloc(GraphFileHeader *h)
{
  size_t i;
  if(h->capacity > 0) {
    for(i = 0; i < h->capacity; i++)
      graph_info_dealloc(h->ginfo + i);
    free(h->ginfo);
    h->capacity = 0;
  }
}

void graph_header_print(const GraphFileHeader *header)
{
  printf("HEADER\n");
  printf("  version: %u\n", header->version);
  printf("  kmer_size: %u\n", header->kmer_size);
  printf("  num_of_bitfields: %u\n", header->num_of_bitfields);
  printf("  num_of_cols: %u\n", header->num_of_cols);
  printf("  num_of_kmers: %zu\n", (size_t)header->num_of_kmers);
  printf("  [capacity: %zu]\n", header->capacity);
}

// Copy non-colour specific values
void graph_header_global_cpy(GraphFileHeader *dst, const GraphFileHeader *src)
{
  dst->version = src->version;
  dst->kmer_size = src->kmer_size;
  dst->num_of_bitfields = src->num_of_bitfields;
  dst->num_of_cols = src->num_of_cols;
  dst->num_of_kmers = src->num_of_kmers;
  graph_header_alloc(dst, src->num_of_cols);
  // size_t i;
  // for(i = 0; i < src->num_of_cols; i++)
  //   graph_info_cpy(&dst->ginfo[i], &src->ginfo[i]);
}

// Merge headers and set intersect name (if intersect_gname != NULL)
void graph_reader_merge_headers(GraphFileHeader *hdr,
                                const GraphFileReader *files, size_t num_files,
                                const char *intersect_gname)
{
  size_t i, j, ncols = 0, intocol, fromcol;

  for(i = 0; i < num_files; i++)
    ncols = MAX2(ncols, graph_file_usedcols(&files[i]));

  graph_header_global_cpy(hdr, &files[0].hdr);
  graph_header_alloc(hdr, ncols);
  hdr->num_of_cols = (uint32_t)ncols;

  for(i = 0; i < num_files; i++) {
    for(j = 0; j < files[i].fltr.ncols; j++) {
      intocol = graph_file_intocol(&files[i], j);
      fromcol = graph_file_fromcol(&files[i], j);
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

// Return number of bytes read or -1 if not valid
// Note: Doesn't set num_of_kmers if version < 7
int graph_file_read_header(FILE *fh, GraphFileHeader *h,
                           boolean fatal, const char *path)
{
  size_t i;
  int bytes_read = 0;

  char magic_word[7];
  magic_word[6] = '\0';

  SAFE_READ(fh, magic_word, strlen("CORTEX"), "Magic word", path, fatal);
  if(strcmp(magic_word, "CORTEX") != 0) {
    if(!fatal) return -1;
    die("Magic word doesn't match '%s' (start): %s", "CORTEX", path);
  }
  bytes_read += strlen("CORTEX");

  // Read version number, kmer_size, num bitfields, num colours
  SAFE_READ(fh, &h->version, sizeof(uint32_t), "graph version", path, fatal);
  SAFE_READ(fh, &h->kmer_size, sizeof(uint32_t), "kmer size", path, fatal);
  SAFE_READ(fh, &h->num_of_bitfields, sizeof(uint32_t), "num of bitfields", path, fatal);
  SAFE_READ(fh, &h->num_of_cols, sizeof(uint32_t), "number of colours", path, fatal);
  bytes_read += 4*sizeof(uint32_t);

  // Checks
  if(h->version > 7 || h->version < 4)
  {
    if(!fatal) return -1;
    die("Sorry, we only support graph file versions 4, 5, 6 & 7 "
        "[version: %u; path: %s]\n", h->version, path);
  }

  if(h->kmer_size % 2 == 0)
  {
    if(!fatal) return -1;
    die("kmer size is not an odd number [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3)
  {
    if(!fatal) return -1;
    die("kmer size is less than three [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_bitfields * 32 < h->kmer_size)
  {
    if(!fatal) return -1;
    die("Not enough bitfields for kmer size "
        "[kmer_size: %u; bitfields: %u; path: %s]\n",
        h->kmer_size, h->num_of_bitfields, path);
  }

  if((h->num_of_bitfields-1)*32 >= h->kmer_size) {
    if(!fatal) return -1;
    die("using more than the minimum number of bitfields [path: %s]\n", path);
  }

  if(h->num_of_cols == 0) {
    if(!fatal) return -1;
    die("number of colours is zero [path: %s]\n", path);
  }

  // graph_header_alloc will only alloc or realloc if it needs to
  graph_header_alloc(h, h->num_of_cols);

  if(h->version >= 7)
  {
    SAFE_READ(fh, &h->num_of_kmers, sizeof(uint64_t), "number of kmers", path, fatal);
    uint32_t tmp;
    SAFE_READ(fh, &tmp, sizeof(uint32_t), "number of shades", path, fatal);
    bytes_read += sizeof(uint64_t) + sizeof(uint32_t);

    if((tmp & 0x7) != 0) {
      warn("Number of shades is not a multiple of 8 [path: %s]", path);
    }
  }

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

      StrBuf *sbuf = &h->ginfo[i].sample_name;
      strbuf_ensure_capacity(sbuf, len);

      safe_fread(fh, sbuf->buff, len, "sample name", path);
      bytes_read += sizeof(uint32_t) + len;

      sbuf->buff[len] = '\0';

      size_t len2 = strlen(sbuf->buff);
      sbuf->len = len2;
      sbuf->buff[len2] = '\0';

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

      safe_fread(fh, &(err_cleaning->clean_snodes_thresh),
                 sizeof(uint32_t), "remove low covg supernode threshold", path);
      safe_fread(fh, &(err_cleaning->clean_kmers_thresh),
                 sizeof(uint32_t), "remove low covg kmer threshold", path);

      bytes_read += 4*sizeof(uint8_t) + 2*sizeof(uint32_t);

      // Fix for old versions with negative thresholds
      if(h->version <= 6) {
        if(!err_cleaning->cleaned_snodes &&
           err_cleaning->clean_snodes_thresh == (uint32_t)-1) {
          err_cleaning->clean_kmers_thresh = 0;
        } else if(!err_cleaning->cleaned_kmers &&
           err_cleaning->clean_kmers_thresh == (uint32_t)-1) {
          err_cleaning->clean_kmers_thresh = 0;
        }
      }

      // Sanity checks
      if(!err_cleaning->cleaned_snodes &&
         err_cleaning->clean_snodes_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for supernodes "
             "when no cleaning was performed [path: %s]", path);

        err_cleaning->clean_snodes_thresh = 0;
      }

      if(!err_cleaning->cleaned_kmers &&
         err_cleaning->clean_kmers_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for nodes "
             "when no cleaning was performed [path: %s]", path);

        err_cleaning->clean_kmers_thresh = 0;
      }

      // Read cleaned against name
      uint32_t len;
      safe_fread(fh, &len, sizeof(uint32_t), "graph name length", path);

      StrBuf *sbuf = &err_cleaning->intersection_name;
      strbuf_ensure_capacity(sbuf, len);

      safe_fread(fh, sbuf->buff, len, "cleaned against graph name", path);
      sbuf->buff[len] = '\0';

      bytes_read += sizeof(uint32_t) + len;

      size_t len2 = strlen(sbuf->buff);
      sbuf->len = len2;
      sbuf->buff[len2] = '\0';

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
    if(!fatal) return -1;
    die("Magic word doesn't match '%s' (end): '%s' [path: %s]\n",
        "CORTEX", magic_word, path);
  }
  bytes_read += strlen("CORTEX");

  // If we haven't set num_of_kmers set it now using file size
  if(h->version < 7) h->num_of_kmers = 0;
  /*
  if(h->version < 7)
  {
    off_t file_size = futil_get_file_size(path);
    size_t bytes_remaining = file_size - bytes_read;
    size_t shade_bytes = 0;

    // 2 * num_shade_bytes for shade + shade end data
    size_t num_bytes_per_kmer
      = sizeof(uint64_t) * NUM_BKMER_WORDS +
        sizeof(uint32_t) * h->num_of_cols + // coverage
        sizeof(uint8_t) * h->num_of_cols + // edges
        sizeof(uint8_t) * shade_bytes * 2; // shades

    h->num_of_kmers = bytes_remaining / num_bytes_per_kmer;

    if(num_bytes_per_kmer * h->num_of_kmers != bytes_remaining) {
      if(!fatal) return -1;
      die("Irregular size of graph file (corrupted?): %s", path);
    }
  }
  */

  return bytes_read;
}

size_t graph_file_read_kmer(FILE *fh, const GraphFileHeader *h, const char *path,
                            BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  size_t i, num_bytes_read;

  num_bytes_read = fread(bkmer->b, 1, sizeof(uint64_t)*h->num_of_bitfields, fh);

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
  if(i == h->num_of_cols)
    warn("Kmer has zero covg in all colours in path: %s", path);

  return num_bytes_read;
}

// Print some output
static void graph_loading_print_status(const GraphFileReader *file)
{
  const GraphFileHeader *hdr = &file->hdr;
  const FileFilter *fltr = &file->fltr;
  char nkmers_str[100], filesize_str[100];

  ulong_to_str(hdr->num_of_kmers, nkmers_str);
  bytes_to_str(fltr->file_size, 1, filesize_str);

  file_filter_status(fltr);
  status("  %s kmers, %s filesize", nkmers_str, filesize_str);
}

// if only_load_if_in_colour is >= 0, only kmers with coverage in existing
// colour only_load_if_in_colour will be loaded.
// We assume only_load_if_in_colour < load_first_colour_into
// if all_kmers_are_unique != 0 an error is thrown if a node already exists
// returns the number of colours in the graph
// If stats != NULL, updates:
//   stats->num_of_colours_loaded
//   stats->kmers_loaded
//   stats->total_bases_read
//   stats->ctx_files_loaded

size_t graph_load(GraphFileReader *file, const GraphLoadingPrefs prefs,
                  LoadingStats *stats)
{
  ctx_assert(!prefs.must_exist_in_graph || prefs.must_exist_in_edges != NULL);

  dBGraph *graph = prefs.db_graph;
  GraphInfo *ginfo = graph->ginfo;
  size_t i, load_ncols = graph_file_outncols(file), fromcol, intocol;
  FileFilter *fltr = &file->fltr;
  GraphFileHeader *hdr = &file->hdr;

  ctx_assert(load_ncols > 0);

  // Print status
  graph_loading_print_status(file);

  fseek(fltr->fh, file->hdr_size, SEEK_SET);

  // Check we can load this graph file into db_graph (kmer size + num colours)
  if(hdr->kmer_size != graph->kmer_size)
  {
    die("Graph has different kmer size [kmer_size: %u vs %zu; path: %s]",
        hdr->kmer_size, graph->kmer_size, fltr->file_path.buff);
  }

  if(fltr->intocol + load_ncols > graph->num_of_cols &&
     (graph->col_covgs != NULL || graph->col_edges != NULL))
  {
    die("Program has not assigned enough colours! "
        "[colours in graph: %zu vs graph: %zu, load into: %zu; path: %s]",
        load_ncols, graph->num_of_cols, fltr->intocol, fltr->file_path.buff);
  }

  for(i = 0; i < fltr->ncols; i++) {
    fromcol = file_filter_fromcol(fltr, i);
    intocol = file_filter_intocol(fltr, i);
    graph_info_merge(ginfo+intocol, hdr->ginfo+fromcol);
  }

  // Update number of colours loaded
  graph->num_of_cols_used = MAX2(graph->num_of_cols_used,
                                 file_filter_usedcols(fltr));

  // Read kmers
  BinaryKmer bkmer;
  Covg covgs[load_ncols];
  Edges edges[load_ncols];

  size_t nkmers_parsed, num_of_kmers_loaded = 0;
  uint64_t num_of_kmers_already_loaded = graph->ht.num_kmers;

  status("[CtxLoad] Reading into %zu colour%s from file with %u...",
         load_ncols, load_ncols != 1 ? "s" : "", hdr->num_of_cols);

  for(nkmers_parsed = 0; graph_file_read(file, &bkmer, covgs, edges); nkmers_parsed++)
  {
    // If kmer has no covg or edges -> don't load
    Covg keep_kmer = 0;
    for(i = 0; i < load_ncols; i++) keep_kmer |= covgs[i] | edges[i];
    if(keep_kmer == 0) continue;

    if(prefs.boolean_covgs)
      for(i = 0; i < load_ncols; i++)
        covgs[i] = covgs[i] > 0;

    // Fetch node in the de bruijn graph
    hkey_t node;

    if(prefs.must_exist_in_graph)
    {
      node = hash_table_find(&graph->ht, bkmer);
      if(node == HASH_NOT_FOUND) continue;

      // Edges union_edges = db_node_get_edges_union(graph, node);
      Edges union_edges = prefs.must_exist_in_edges[node];

      for(i = 0; i < load_ncols; i++) edges[i] &= union_edges;
    }
    else
    {
      boolean found;
      node = hash_table_find_or_insert(&graph->ht, bkmer, &found);

      if(prefs.empty_colours && found)
        die("Duplicate kmer loaded [cols:%zu:%zu]", fltr->intocol, load_ncols);
    }

    // Set presence in colours
    uint8_t has_col;
    if(graph->node_in_cols != NULL) {
      for(i = 0; i < load_ncols; i++) {
        has_col = (covgs[i] > 0 || edges[i] != 0);
        intocol = graph_file_intocol(file,i);
        db_node_cpy_col(graph, node, intocol, has_col);
      }
    }

    if(graph->col_covgs != NULL) {
      for(i = 0; i < load_ncols; i++)
        db_node_add_col_covg(graph, node, graph_file_intocol(file,i), covgs[i]);
    }

    // This may be an invalid pointer (if num_edge_cols == 0)
    Edges *col_edges = graph->col_edges + node * graph->num_edge_cols;

    // Merge all edges into one colour
    if(graph->num_edge_cols == 1) {
      for(i = 0; i < load_ncols; i++)
        col_edges[0] |= edges[i];
    }
    else if(graph->num_edge_cols > 0) {
      for(i = 0; i < load_ncols; i++)
        col_edges[graph_file_intocol(file,i)] |= edges[i];
    }

    num_of_kmers_loaded++;
  }

  if(nkmers_parsed != hdr->num_of_kmers)
  {
    warn("More kmers in graph than expected [expected: %zu; actual: %zu; "
         "path: %s]", (size_t)hdr->num_of_kmers, nkmers_parsed,
        fltr->file_path.buff);
  }

  if(stats != NULL)
  {
    stats->num_of_colours_loaded += load_ncols;
    stats->kmers_loaded += num_of_kmers_loaded;
    stats->num_kmers += graph->ht.num_kmers - num_of_kmers_already_loaded;
    for(i = 0; i < load_ncols; i++)
      stats->total_bases_read += hdr->ginfo[i].total_sequence;
    stats->ctx_files_loaded++;
  }

  char parsed_nkmers_str[100], loaded_nkmers_str[100];
  double loaded_nkmers_pct = (100.0 * num_of_kmers_loaded) / nkmers_parsed;
  ulong_to_str(num_of_kmers_loaded, loaded_nkmers_str);
  ulong_to_str(nkmers_parsed, parsed_nkmers_str);
  status("Loaded %s / %s (%.2f%%) of kmers parsed",
         loaded_nkmers_str, parsed_nkmers_str, loaded_nkmers_pct);

  return num_of_kmers_loaded;
}

size_t graph_load_colour(GraphFileReader *file,
                         const GraphLoadingPrefs prefs,
                         LoadingStats *stats,
                         size_t colour_idx, size_t intocol)
{
  size_t *tmpcols, tmpncols, tmpinto, newcol;
  boolean tmpflatten;
  FileFilter *fltr = &file->fltr;

  ctx_assert(colour_idx < fltr->ncols);
  ctx_assert(intocol < prefs.db_graph->num_of_cols);

  // Copy current values
  tmpcols = fltr->cols; tmpncols = fltr->ncols; tmpinto = fltr->intocol;
  tmpflatten = fltr->flatten;

  // Set new values
  newcol = fltr->cols[colour_idx];
  fltr->cols = &newcol; fltr->ncols = 1;
  file_filter_update_intocol(fltr, intocol);

  // Load graph
  size_t kmers_loaded = graph_load(file, prefs, stats);

  // Restore values
  fltr->cols = tmpcols; fltr->ncols = tmpncols;
  file_filter_update_intocol(fltr, tmpinto);
  fltr->flatten = tmpflatten;

  return kmers_loaded;
}

//
// Merging, filtering, combining graph files
//

// Load a kmer and write to a file one kmer at a time
// Optionally filter a against the graph currently loaded
//   (i.e. only keep nodes and edges that are in the graph)
// Same functionality as graph_files_merge, but faster if dealing with only one
// input file. Reads in and dumps one kmer at a time
// parameters:
//   `flatten`: if true merge colours into one
//   `only_load_if_in_edges`: Edges to mask edges with, 1 per hash table entry
size_t graph_stream_filter(const char *out_ctx_path, const GraphFileReader *file,
                           const dBGraph *db_graph, const GraphFileHeader *hdr,
                           const Edges *only_load_if_in_edges)
{
  const FileFilter *fltr = &file->fltr;
  boolean only_load_if_in_graph = (only_load_if_in_edges != NULL);
  status("Filtering %s to %s with stream filter", fltr->file_path.buff, out_ctx_path);

  FILE *out;
  if((out = fopen(out_ctx_path, "w")) == NULL)
    die("Cannot open output path: %s", out_ctx_path);
  setvbuf(out, NULL, _IOFBF, CTX_BUF_SIZE);

  graph_loading_print_status(file);

  size_t i, nodes_dumped = 0, ncols = graph_file_outncols(file);

  graph_write_header(out, hdr);

  BinaryKmer bkmer;
  Covg kmercovgs[fltr->intocol+ncols], *covgs = kmercovgs+fltr->intocol;
  Edges kmeredges[fltr->intocol+ncols], *edges = kmeredges+fltr->intocol;

  memset(kmercovgs, 0, sizeof(Covg)*(fltr->intocol+ncols));
  memset(kmeredges, 0, sizeof(Edges)*(fltr->intocol+ncols));

  while(graph_file_read(file, &bkmer, covgs, edges))
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
        graph_write_kmer(out, hdr, bkmer.b, kmercovgs, kmeredges);
        nodes_dumped++;
      }
    }
  }

  fclose(out);

  graph_write_status(nodes_dumped, hdr->num_of_cols,
                     out_ctx_path, CTX_GRAPH_FILEFORMAT);

  return nodes_dumped;
}

size_t graph_stream_filter_mkhdr(const char *out_ctx_path, GraphFileReader *file,
                                 const dBGraph *db_graph,
                                 const Edges *only_load_if_in_edges,
                                 const char *intersect_gname)
{
  ctx_assert(intersect_gname == NULL || db_graph->col_edges != NULL);
  ctx_assert(intersect_gname == NULL || only_load_if_in_edges != NULL);

  size_t i, nodes_dumped, ncols = graph_file_outncols(file);
  GraphFileHeader outheader = INIT_GRAPH_FILE_HDR;
  FileFilter *fltr = &file->fltr;

  graph_header_global_cpy(&outheader, &file->hdr);
  outheader.num_of_cols = (uint32_t)(fltr->intocol + ncols);
  graph_header_alloc(&outheader, outheader.num_of_cols);

  for(i = 0; i < fltr->ncols; i++) {
    GraphInfo *ginfo = &outheader.ginfo[graph_file_intocol(file, i)];
    graph_info_merge(ginfo, file->hdr.ginfo + fltr->cols[i]);
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
                         boolean kmers_loaded, boolean colours_loaded,
                         const Edges *only_load_if_in_edges,
                         GraphFileHeader *hdr, dBGraph *db_graph)
{
  boolean only_load_if_in_graph = (only_load_if_in_edges != NULL);
  ctx_assert(!only_load_if_in_graph || kmers_loaded);
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  ctx_assert(!colours_loaded || kmers_loaded);
  ctx_assert(hdr != NULL);

  size_t i, ncols, output_colours = 0;

  for(i = 0; i < num_files; i++) {
    ncols = graph_file_usedcols(&files[i]);
    output_colours = MAX2(output_colours, ncols);

    if(files[i].hdr.kmer_size != files[0].hdr.kmer_size) {
      die("Kmer-size mismatch %u vs %u [%s vs %s]",
          files[0].hdr.kmer_size, files[i].hdr.kmer_size,
          files[0].fltr.file_path.buff, files[i].fltr.file_path.buff);
    }
  }

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

  LoadingStats stats;
  loading_stats_init(&stats);

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
    // Have to load a few colours at a time then dump, rinse and repeat
    status("Saving %zu colours, %zu colours at a time",
           output_colours, db_graph->num_of_cols);

    // Open file, write header
    FILE *fout = fopen(out_ctx_path, "w");
    if(fout == NULL) die("Cannot open output ctx file: %s", out_ctx_path);
    setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);
    size_t header_size = graph_write_header(fout, hdr);

    // Load all kmers into flat graph
    if(!kmers_loaded)
    {
      size_t tmpinto; boolean tmpflatten;
      for(i = 0; i < num_files; i++)
      {
        tmpinto = files[i].fltr.intocol; tmpflatten = files[i].fltr.flatten;
        // files[i].fltr.intocol = 0;
        file_filter_update_intocol(&files[i].fltr, 0);
        files[i].fltr.flatten = true;
        graph_load(&files[i], prefs, &stats);
        // files[i].fltr.intocol = tmpinto;
        file_filter_update_intocol(&files[i].fltr, tmpinto);
        files[i].fltr.flatten = tmpflatten;
      }
    }

    // print file outline
    status("Generated merged hash table\n");
    hash_table_print_stats(&db_graph->ht);
    graph_write_empty(db_graph, fout, output_colours);

    size_t num_kmer_cols = db_graph->ht.capacity * db_graph->num_of_cols;
    size_t firstcol, lastcol, file_lastcol;

    for(firstcol = 0; firstcol < output_colours; firstcol += db_graph->num_of_cols)
    {
      boolean loaded = false;
      lastcol = MIN2(firstcol + db_graph->num_of_cols - 1, output_colours-1);

      // Wipe colour coverages and edges
      memset(db_graph->col_edges, 0, num_kmer_cols * sizeof(Edges));
      memset(db_graph->col_covgs, 0, num_kmer_cols * sizeof(Covg));

      // Check if any files cover this colour range
      for(i = 0; i < num_files; i++) {
        file_lastcol = graph_file_usedcols(&files[i]) - 1;
        if(files[i].fltr.intocol <= lastcol && file_lastcol >= firstcol) break;
      }

      if(i == num_files) {
        // No hits
        size_t newstrt = SIZE_MAX;
        for(i = 0; i < num_files; i++)
          if(files[i].fltr.intocol > firstcol)
            newstrt = MIN2(newstrt, files[i].fltr.intocol);
        firstcol = newstrt;
        lastcol = MIN2(firstcol + db_graph->num_of_cols - 1, output_colours-1);
        ctx_assert(newstrt < SIZE_MAX);
      }

      for(i = 0; i < num_files; i++)
      {
        file_lastcol = graph_file_usedcols(&files[i]) - 1;
        if(files[i].fltr.intocol <= lastcol && file_lastcol >= firstcol)
        {
          // Backup intocol, ncols, cols
          size_t tmpinto = files[i].fltr.intocol, tmpncols = files[i].fltr.ncols;
          size_t *tmpcols = files[i].fltr.cols;

          // Modify file filter to only load limited colours
          if(files[i].fltr.flatten)
            files[i].fltr.intocol -= firstcol;
          else {
            if(files[i].fltr.intocol < firstcol) {
              files[i].fltr.cols += firstcol - files[i].fltr.intocol;
              files[i].fltr.ncols -= firstcol - files[i].fltr.intocol;
              // files[i].fltr.intocol = 0;
              file_filter_update_intocol(&files[i].fltr, 0);
            }
            else files[i].fltr.intocol -= firstcol;

            files[i].fltr.ncols = MIN2(db_graph->num_of_cols-files[i].fltr.intocol,
                                  files[i].fltr.ncols);
          }

          fseek(files[i].fltr.fh, files[i].hdr_size, SEEK_SET);
          graph_load(&files[i], prefs, &stats);
          loaded = true;

          // files[i].fltr.intocol = tmpinto;
          file_filter_update_intocol(&files[i].fltr, tmpinto);
          files[i].fltr.cols = tmpcols;
          files[i].fltr.ncols = tmpncols;
        }
      }

      // if loaded, dump
      if(loaded) {
        if(db_graph->num_of_cols == 1)
          status("Dumping into colour %zu...\n", firstcol);
        else
          status("Dumping into colours %zu-%zu...\n", firstcol, lastcol);
        fseek(fout, (long)header_size, SEEK_SET);
        graph_file_write_colours(db_graph, 0, firstcol, lastcol-firstcol+1,
                                 output_colours, fout);
      }
    }

    fclose(fout);

    graph_write_status(db_graph->ht.num_kmers, output_colours,
                       out_ctx_path, CTX_GRAPH_FILEFORMAT);
  }

  return db_graph->ht.num_kmers;
}

// if intersect_gname != NULL: only load kmers that are already in the hash table
//    and use string as name for cleaning against
// returns the number of kmers written
size_t graph_files_merge_mkhdr(const char *out_ctx_path,
                               GraphFileReader *files, size_t num_files,
                               boolean kmers_loaded, boolean colours_loaded,
                               const Edges *only_load_if_in_edges,
                               const char *intersect_gname, dBGraph *db_graph)
{
  size_t num_kmers;
  GraphFileHeader gheader = INIT_GRAPH_FILE_HDR;

  graph_reader_merge_headers(&gheader, files, num_files, intersect_gname);

  num_kmers = graph_files_merge(out_ctx_path, files, num_files,
                                 kmers_loaded, colours_loaded,
                                 only_load_if_in_edges,
                                 &gheader, db_graph);

  graph_header_dealloc(&gheader);
  return num_kmers;
}
