#include "global.h"
#include "graph_writer.h"
#include "graphs_load.h" // need to load, merge then write some graphs
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"

// Construct graph header
// Free with graph_header_free(hdr)
GraphFileHeader* graph_writer_mkhdr(const dBGraph *db_graph,
                                    const FileFilter *fltr,
                                    size_t filencols)
{
  size_t i, from, into;
  GraphFileHeader *hdr = ctx_calloc(1, sizeof(*hdr));
  hdr->version = CTX_GRAPH_FILEFORMAT;
  hdr->kmer_size = (uint32_t)db_graph->kmer_size;
  hdr->num_of_bitfields = NUM_BKMER_WORDS;
  hdr->num_of_cols = (uint32_t)filencols;
  graph_header_capacity(hdr, filencols);

  for(i = 0; i < file_filter_num(fltr); i++) {
    from = file_filter_fromcol(fltr,i);
    into = file_filter_intocol(fltr, i);
    graph_info_merge(&hdr->ginfo[into], &db_graph->ginfo[from]);
  }

  return hdr;
}

// Returns number of bytes written
static size_t write_error_cleaning_object(FILE *fh, const ErrorCleaning *cleaning)
{
  size_t written = 0, expwrite;
  written += fwrite(&(cleaning->cleaned_tips), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->cleaned_unitigs), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->cleaned_kmers), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->is_graph_intersection), 1, sizeof(uint8_t), fh);

  uint32_t clean_unitigs_thresh
    = cleaning->cleaned_unitigs ? cleaning->clean_unitigs_thresh : 0;

  uint32_t clean_kmers_thresh
    = cleaning->cleaned_kmers ? cleaning->clean_kmers_thresh : 0;

  written += fwrite(&clean_unitigs_thresh, 1, sizeof(uint32_t), fh);
  written += fwrite(&clean_kmers_thresh, 1, sizeof(uint32_t), fh);

  uint32_t len = (uint32_t)cleaning->intersection_name.end;
  const char *str = cleaning->intersection_name.b;
  written += fwrite(&len, 1, sizeof(uint32_t), fh);
  written += fwrite(str, 1, len, fh);

  expwrite = 4 + sizeof(uint32_t) * 2 + sizeof(uint32_t) + len;
  if(written != expwrite) die("Cannot write to file");

  return written;
}

// Returns number of bytes written
size_t graph_write_header(FILE *fh, const GraphFileHeader *h)
{
  size_t i, b = 0, act = 0, tmp;

  act += fwrite("CORTEX", 1, strlen("CORTEX"), fh);
  act += fwrite(&h->version, 1, sizeof(uint32_t), fh);
  act += fwrite(&h->kmer_size, 1, sizeof(uint32_t), fh);
  act += fwrite(&h->num_of_bitfields, 1, sizeof(uint32_t), fh);
  act += fwrite(&h->num_of_cols, 1, sizeof(uint32_t), fh);

  b += strlen("CORTEX") + sizeof(uint32_t) * 4;

  for(i = 0; i < h->num_of_cols; i++)
    act += fwrite(&h->ginfo[i].mean_read_length, 1, sizeof(uint32_t), fh);
  for(i = 0; i < h->num_of_cols; i++)
    act += fwrite(&h->ginfo[i].total_sequence, 1, sizeof(uint64_t), fh);

  b += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint64_t));

  if(h->version >= 6)
  {
    for(i = 0; i < h->num_of_cols; i++)
    {
      uint32_t len = (uint32_t)h->ginfo[i].sample_name.end;
      const char *buff = h->ginfo[i].sample_name.b;
      act += fwrite(&len, 1, sizeof(uint32_t), fh);
      act += fwrite(buff, 1, len, fh);
      b += sizeof(uint32_t) + len;
    }

    // This fwrite seems to be upsetting valgrind - it doesn't like long doubles
    for(i = 0; i < h->num_of_cols; i++)
      act += fwrite(&h->ginfo[i].seq_err, 1, sizeof(long double), fh);

    b += h->num_of_cols * sizeof(long double);

    for(i = 0; i < h->num_of_cols; i++) {
      tmp = write_error_cleaning_object(fh, &h->ginfo[i].cleaning);
      b += tmp; act += tmp;
    }
  }

  act += fwrite("CORTEX", 1, strlen("CORTEX"), fh);
  b += strlen("CORTEX");

  if(act != b) die("Cannot write file");

  return b;
}


// Write a single kmer with its edges and coverages for colours 0..ncols-1 to
// 0..ncols-1 in the file. `covgs` and `edges` should be of length `filencols`.
// Returns number of bytes written
size_t graph_write_kmer(FILE *fh, size_t filencols,
                        const BinaryKmer bkmer,
                        const Covg *covgs,
                        const Edges *edges)
{
  size_t m = 0, expm = BKMER_BYTES + 5*filencols;
  m += fwrite(bkmer.b, 1, BKMER_BYTES, fh);
  m += fwrite(covgs, 1, sizeof(uint32_t) * filencols, fh);
  m += fwrite(edges, 1, sizeof(uint8_t) * filencols, fh);
  if(m != expm) die("Cannot write to file (%zu, %zu)", m, expm);
  return m;
}

// Write kmer with no re-ordering of colours
static inline void graph_write_kmer_direct(hkey_t hkey,
                                           const GraphFileHeader *hdr,
                                           FILE *fh, const dBGraph *db_graph)
{
  graph_write_kmer(fh, hdr->num_of_cols,
                   db_graph->ht.table[hkey],
                   &db_node_covg(db_graph, hkey, 0),
                   &db_node_edges(db_graph, hkey, 0));
}


// Dump node: only print kmers with coverages in given colours
static void graph_write_kmer_indirect(hkey_t hkey, const GraphFileHeader *hdr,
                                      const FileFilter *fltr, FILE *fh,
                                      const dBGraph *db_graph,
                                      size_t *num_dumped)
{
  size_t i, into, from;
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  Covg covgs[hdr->num_of_cols], merge_covgs = 0;
  Edges edges[hdr->num_of_cols], merge_edges = 0;

  memset(covgs, 0, sizeof(Covg) * hdr->num_of_cols);
  memset(edges, 0, sizeof(Edges) * hdr->num_of_cols);

  Edges (*col_edges)[db_graph->num_of_cols]
    = (Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  Covg (*col_covgs)[db_graph->num_of_cols]
    = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  for(i = 0; i < file_filter_num(fltr); i++) {
    into = file_filter_intocol(fltr, i);
    from = file_filter_fromcol(fltr, i);
    SAFE_SUM_COVG(covgs[into], col_covgs[hkey][from]);
    edges[into] |= col_edges[hkey][from];
    merge_covgs |= covgs[into];
    merge_edges |= edges[into];
  }

  // Check this node has coverage in one of the specified colours
  if(merge_covgs > 0) {
    graph_write_kmer(fh, hdr->num_of_cols, bkmer, covgs, edges);
    (*num_dumped)++;
  }
}


// Dump all kmers with all colours to given file.
// Write kmer with no re-ordering of colours
// write the first `ncols` from the graph in memory to the file handle
// `sort_kmer` if true, sort kmers before writing. Uses extra memory.
// Returns num of kmers written
size_t graph_write_all_kmers_direct(FILE *fh, const dBGraph *db_graph,
                                    bool sort_kmers, const GraphFileHeader *hdr)
{
  if(sort_kmers) {
    HASH_ITERATE_SORTED(&db_graph->ht, graph_write_kmer_direct,
                        hdr, fh, db_graph);
  } else {
    HASH_ITERATE(&db_graph->ht, graph_write_kmer_direct,
                 hdr, fh, db_graph);
  }
  return db_graph->ht.num_kmers;
}

size_t graph_write_all_kmers_filtered(FILE *fh, const dBGraph *db_graph,
                                      bool sort_kmers, const GraphFileHeader *hdr,
                                      const FileFilter *fltr)
{
  size_t num_nodes_dumped = 0;
  if(sort_kmers) {
    HASH_ITERATE_SORTED(&db_graph->ht, graph_write_kmer_indirect,
                        hdr, fltr, fh, db_graph,
                        &num_nodes_dumped);
  } else {
    HASH_ITERATE(&db_graph->ht, graph_write_kmer_indirect,
                 hdr, fltr, fh, db_graph,
                 &num_nodes_dumped);
  }
  return num_nodes_dumped;
}

// Pass your own header
// If sort_kmers is true, save kmers in lexigraphical order
// returns number of nodes written out
uint64_t graph_writer_save(const char *path, const dBGraph *db_graph,
                           const GraphFileHeader *hdr, bool sort_kmers,
                           const FileFilter *fltr)
{
  ctx_assert(db_graph->col_edges != NULL);
  ctx_assert(db_graph->col_covgs != NULL);

  uint64_t n_nodes = 0;
  const char *out_name = futil_outpath_str(path);

  status("[graphwriter] Saving file to: %s", path);
  file_filter_status(fltr, true);

  FILE *fh = futil_fopen(path, "w");

  // Write header
  graph_write_header(fh, hdr);

  if(file_filter_into_direct(fltr,hdr->num_of_cols)) {
    n_nodes = graph_write_all_kmers_direct(fh, db_graph, sort_kmers, hdr);
  }
  else {
    n_nodes = graph_write_all_kmers_filtered(fh, db_graph, sort_kmers, hdr,
                                             fltr);
  }

  fclose(fh);

  graph_writer_print_status(n_nodes, hdr->num_of_cols,
                            out_name, hdr->version);

  return n_nodes;
}

// filencols must be <= db_graph->num_of_cols
// Return number of kmers written out
uint64_t graph_writer_save_mkhdr(const char *path, const dBGraph *db_graph,
                                 bool sort_kmers,  size_t filencols)
{
  ctx_assert(filencols <= db_graph->num_of_cols);

  // Create a filter
  FileFilter fltr;
  memset(&fltr, 0, sizeof(fltr));
  file_filter_create_direct(&fltr, db_graph->num_of_cols, filencols);

  // Construct graph header
  GraphFileHeader *hdr = graph_writer_mkhdr(db_graph, &fltr, filencols);
  uint64_t nkmers = graph_writer_save(path, db_graph, hdr, sort_kmers, &fltr);

  graph_header_free(hdr);
  file_filter_close(&fltr);
  return nkmers;
}

void graph_writer_print_status(uint64_t nkmers, size_t ncols,
                               const char *path, uint32_t version)
{
  char num_kmer_str[100];
  ulong_to_str(nkmers, num_kmer_str);

  status("[graphwriter] Dumped %s kmers in %zu colour%s into: %s (ver: %u)\n",
         num_kmer_str, ncols, util_plural_str(ncols),
         futil_outpath_str(path), version);
}

//
// Merging, filtering, combining graph files
//

//
// Write empty file then update
//
static inline void _dump_empty_bkmer(hkey_t hkey, const dBGraph *db_graph,
                                     char *buf, size_t mem, FILE *fh)
{
  size_t written;
  const BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  written = fwrite(bkmer.b, 1, sizeof(BinaryKmer), fh) +
            fwrite(buf, 1, mem, fh);
  if(written != mem+sizeof(BinaryKmer)) die("Couldn't write to file");
}

/*!
  Write kmers from the graph to a file. The file header should already have been
  written.
  @return Number of bytes written
 */
static size_t graph_write_empty(const dBGraph *db_graph, FILE *fh,
                                bool sort_kmers, size_t num_of_cols)
{
  size_t mem = num_of_cols * (sizeof(Covg)+sizeof(Edges));
  char buf[mem];
  memset(buf, 0, mem);
  if(sort_kmers) {
    HASH_ITERATE_SORTED(&db_graph->ht, _dump_empty_bkmer, db_graph, buf, mem, fh);
  } else {
    HASH_ITERATE(&db_graph->ht, _dump_empty_bkmer, db_graph, buf, mem, fh);
  }
  return db_graph->ht.num_kmers * (sizeof(BinaryKmer) + mem);
}

/*!
  Overwrite kmers in an existing file.
  @param unordered_kmers file kmer order does not match hash table order
  @param first_graphcol  first colour in the dBGraph to read from
  @param first_filecol   first colour in the file to write into
  @param ngraphcols      number of colours to write to file
  @param nfilecols       total number of colours in file
  @param hdrsize         file header size i.e. byte pos of first kmer in file
  @param fh              file handle opened that we can fwrite/fseek
  @param path            path to the file (for error messages)
 */
static void graph_writer_update_file_kmers(const dBGraph *db_graph,
                                           bool unordered_kmers,
                                           size_t first_graphcol, size_t ngraphcols,
                                           size_t first_filecol, size_t nfilecols,
                                           size_t hdrsize, FILE *fh, const char *path)
{
  ctx_assert(db_graph->col_edges != NULL);
  ctx_assert(db_graph->col_covgs != NULL);
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  ctx_assert(first_graphcol+ngraphcols <= db_graph->num_of_cols);

  // db_graph->ht.num_kmers is also the number of kmers in the file
  // We are just overwriting some of the coverages and edges for colours
  // first_filecol..(first_filecol+ngraphcols) not including last
  size_t max_block_size = 16 * ONE_MEGABYTE;
  size_t filekmersize = sizeof(BinaryKmer)+(sizeof(Edges)+sizeof(Covg))*nfilecols;
  size_t kmers_per_block = max_block_size / filekmersize;
  size_t block_size = kmers_per_block * filekmersize;
  ctx_assert(block_size > 0);

  size_t nkmers_printed = 0, nkmers, nbytes, end;
  uint8_t *mem = ctx_malloc(block_size), *memptr;
  hkey_t hkey = 0;
  const Covg *covgs = db_graph->col_covgs;
  const Edges *edges = db_graph->col_edges;
  BinaryKmer bkmer;

  if(fseek(fh, hdrsize, SEEK_SET) != 0) die("Cannot seek to file start: %s", path);

  // Read block of kmers from the file
  // iterate over the hash table figuring which ones they are
  while(nkmers_printed < db_graph->ht.num_kmers)
  {
    nkmers = MIN2(db_graph->ht.num_kmers - nkmers_printed, kmers_per_block);
    nbytes = nkmers*filekmersize;
    if(fread(mem, 1, nbytes, fh) != nbytes) die("Cannot read: %s", path);
    memptr = mem;
    // Loop over the kmers we read and update them
    for(end = nkmers_printed+nkmers; nkmers_printed < end; nkmers_printed++, hkey++) {
      if(unordered_kmers) { // need to find kmer
        memcpy(&bkmer.b, memptr, sizeof(BinaryKmer));
        hkey = hash_table_find(&db_graph->ht, bkmer);
        ctx_assert(hkey != HASH_NOT_FOUND);
      }
      else { // linear search of the hash table (it's fast!)
        while(!db_graph_node_assigned(db_graph, hkey)) hkey++;
      }
      covgs = db_graph->col_covgs + hkey * db_graph->num_of_cols;
      edges = db_graph->col_edges + hkey * db_graph->num_of_cols;
      memptr += sizeof(BinaryKmer);
      memcpy(memptr + first_filecol*sizeof(Covg), covgs, ngraphcols*sizeof(Covg));
      memptr += sizeof(Covg)*nfilecols;
      memcpy(memptr + first_filecol*sizeof(Edges), edges, ngraphcols*sizeof(Edges));
      memptr += sizeof(Edges)*nfilecols;
    }
    if(fseek(fh, -nbytes, SEEK_CUR) != 0) die("fseek failed: %s", path);
    if(fwrite(mem, 1, nbytes, fh) != nbytes) die("fwrite failed: %s", path);
  }

  ctx_free(mem);
}


// Load a kmer and write to a file one kmer at a time
// Optionally filter against the graph currently loaded
//   (i.e. only keep nodes and edges that are in the graph)
// Same functionality as graph_writer_merge, but faster if dealing with only one
// input file. Reads in and dumps one kmer at a time
// parameters:
//   `only_load_if_in_edges`: Edges to mask edges with, 1 per hash table entry
size_t graph_writer_stream(const char *out_ctx_path, GraphFileReader *file,
                           const dBGraph *db_graph, const GraphFileHeader *hdr,
                           const Edges *only_load_if_in_edges)
{
  const FileFilter *fltr = &file->fltr;
  bool only_load_if_in_graph = (only_load_if_in_edges != NULL);
  status("Filtering %s to %s with stream filter", fltr->path.b,
         futil_outpath_str(out_ctx_path));
  graph_loading_print_status(file);

  // seek to start of input file after header
  if(graph_file_fseek(file, file->hdr_size, SEEK_SET) != 0)
    die("fseek failed: %s", strerror(errno));

  FILE *out = futil_fopen(out_ctx_path, "w");
  graph_write_header(out, hdr);

  size_t i, nodes_dumped = 0, ncols = file_filter_into_ncols(fltr);

  BinaryKmer bkmer;
  Covg covgs[ncols];
  Edges edges[ncols];

  while(graph_file_read_reset(file, &bkmer, covgs, edges))
  {
    // Collapse down colours
    Covg keep_kmer = 0;
    for(i = 0; i < ncols; i++) keep_kmer |= covgs[i];

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
        graph_write_kmer(out, hdr->num_of_cols, bkmer, covgs, edges);
        nodes_dumped++;
      }
    }
  }

  fflush(out);
  fclose(out);

  graph_writer_print_status(nodes_dumped, hdr->num_of_cols,
                            out_ctx_path, hdr->version);

  return nodes_dumped;
}

size_t graph_writer_stream_mkhdr(const char *out_ctx_path,
                                 GraphFileReader *file,
                                 const dBGraph *db_graph,
                                 const Edges *only_load_if_in_edges,
                                 const char *intersect_gname)
{
  ctx_assert(intersect_gname == NULL || db_graph->col_edges != NULL);
  ctx_assert(intersect_gname == NULL || only_load_if_in_edges != NULL);

  size_t i, nodes_dumped;

  GraphFileHeader outhdr;
  memset(&outhdr, 0, sizeof(outhdr));
  graph_file_merge_header(&outhdr, file);

  for(i = 0; i < outhdr.num_of_cols; i++)
    if(intersect_gname != NULL)
      graph_info_append_intersect(&outhdr.ginfo[i].cleaning, intersect_gname);

  nodes_dumped = graph_writer_stream(out_ctx_path, file,
                                     db_graph, &outhdr,
                                     only_load_if_in_edges);
  graph_header_dealloc(&outhdr);

  return nodes_dumped;
}

// `kmers_loaded`: means all kmers to dump have been loaded
// `colours_loaded`: means all kmer data have been loaded
// `only_load_if_in_edges`: Edges to mask edges with, 1 per hash table entry
//   should already be set
size_t graph_writer_merge(const char *out_ctx_path,
                          GraphFileReader *files, size_t num_files,
                          bool kmers_loaded, bool colours_loaded,
                          const Edges *only_load_if_in_edges,
                          GraphFileHeader *hdr, bool sort_kmers,
                          dBGraph *db_graph)
{
  bool only_load_if_in_graph = (only_load_if_in_edges != NULL);
  ctx_assert(!only_load_if_in_graph || kmers_loaded);
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  ctx_assert(!colours_loaded || kmers_loaded);
  ctx_assert(hdr != NULL);

  size_t i, j, f, out_ncols = 0;

  for(i = 0; i < num_files; i++) {
    out_ncols = MAX2(out_ncols, file_filter_into_ncols(&files[i].fltr));
    if(files[i].hdr.kmer_size != files[0].hdr.kmer_size) {
      die("Kmer-size mismatch %u vs %u [%s vs %s]",
          files[0].hdr.kmer_size, files[i].hdr.kmer_size,
          files[0].fltr.path.b, files[i].fltr.path.b);
    }
  }

  ctx_assert(0 < out_ncols && out_ncols <= hdr->num_of_cols);

  GraphLoadingPrefs gprefs = graph_loading_prefs(db_graph);
  gprefs.must_exist_in_graph = only_load_if_in_graph;
  gprefs.must_exist_in_edges = only_load_if_in_edges;

  if(kmers_loaded && colours_loaded)
  {
    uint64_t nnodes;
    FileFilter fltr;
    memset(&fltr, 0, sizeof(fltr));
    file_filter_create_direct(&fltr, db_graph->num_of_cols, hdr->num_of_cols);
    nnodes = graph_writer_save(out_ctx_path, db_graph, hdr, sort_kmers, &fltr);
    file_filter_close(&fltr);
    return nnodes;
  }
  else if(num_files == 1 && !sort_kmers)
  {
    return graph_writer_stream(out_ctx_path, &files[0], db_graph, hdr,
                               only_load_if_in_edges);
  }
  else if(hdr->num_of_cols <= db_graph->num_of_cols)
  {
    // Can load all files at once
    status("Loading and saving %zu colours at once", (size_t)hdr->num_of_cols);

    if(!kmers_loaded) {
      for(i = 0; i < num_files; i++)
        graph_load(&files[i], gprefs, NULL);
      hash_table_print_stats(&db_graph->ht);
    }

    uint64_t nnodes;
    FileFilter fltr;
    memset(&fltr, 0, sizeof(fltr));
    file_filter_create_direct(&fltr, db_graph->num_of_cols, hdr->num_of_cols);
    nnodes = graph_writer_save(out_ctx_path, db_graph, hdr, sort_kmers, &fltr);
    file_filter_close(&fltr);
    return nnodes;
  }
  else
  {
    ctx_assert2(strcmp(out_ctx_path,"-") != 0,
                "Cannot use STDOUT for output if not enough colours to load");

    // Have to load a few colours at a time then dump, rinse and repeat
    status("[overwriting] Saving %zu colours, %zu colours at a time",
           out_ncols, db_graph->num_of_cols);

    // Open file, write header
    FILE *fout = futil_fopen(out_ctx_path, "r+");

    size_t hdr_size = graph_write_header(fout, hdr);

    // Load all kmers into flat graph
    if(!kmers_loaded)
      graphs_load_files_flat(files, num_files, gprefs, NULL);

    // print file outline
    status("Generated merged hash table\n");
    hash_table_print_stats(&db_graph->ht);

    // Write empty file
    size_t file_len = hdr_size;
    file_len += graph_write_empty(db_graph, fout, sort_kmers, out_ncols);
    fflush(fout);

    size_t num_kmer_cols = db_graph->ht.capacity * db_graph->num_of_cols;
    size_t firstcol, lastcol, fromcol, intocol;

    FileFilter origfltr;
    memset(&origfltr, 0, sizeof(origfltr));
    bool files_loaded = false;

    for(firstcol = 0; firstcol < out_ncols; firstcol += db_graph->num_of_cols)
    {
      lastcol = MIN2(firstcol + db_graph->num_of_cols - 1, out_ncols-1);

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
        file_filter_update(fltr);

        if(file_filter_num(fltr) > 0)
          files_loaded |= (graph_load(&files[f], gprefs, NULL) > 0);

        // Restore original filter
        file_filter_copy(fltr, &origfltr);
      }

      // if files_loaded, dump
      if(files_loaded) {
        if(db_graph->num_of_cols == 1)
          status("[graphwriter] Dumping into colour %zu...\n", firstcol);
        else
          status("[graphwriter] Dumping into colours %zu-%zu...\n", firstcol, lastcol);

        graph_writer_update_file_kmers(db_graph, sort_kmers,
                                       0, lastcol-firstcol+1,
                                       firstcol, out_ncols,
                                       hdr_size, fout, out_ctx_path);
      }
    }

    fclose(fout);
    file_filter_close(&origfltr);

    // Print output status
    graph_writer_print_status(db_graph->ht.num_kmers, hdr->num_of_cols,
                              out_ctx_path, hdr->version);
  }

  return db_graph->ht.num_kmers;
}

// if intersect_gname != NULL: only load kmers that are already in the hash table
//    and use string as name for cleaning against
// returns the number of kmers written
size_t graph_writer_merge_mkhdr(const char *out_ctx_path,
                               GraphFileReader *files, size_t num_files,
                               bool kmers_loaded, bool colours_loaded,
                               const Edges *only_load_if_in_edges,
                               const char *intersect_gname,
                               bool sort_kmers, dBGraph *db_graph)
{
  size_t i, num_kmers;
  GraphFileHeader hdr;
  memset(&hdr, 0, sizeof(hdr));

  for(i = 0; i < num_files; i++)
    graph_file_merge_header(&hdr, &files[i]);

  if(intersect_gname != NULL) {
    for(i = 0; i < hdr.num_of_cols; i++)
      if(graph_file_is_colour_loaded(i, files, num_files))
        graph_info_append_intersect(&hdr.ginfo[i].cleaning, intersect_gname);
  }

  num_kmers = graph_writer_merge(out_ctx_path, files, num_files,
                                 kmers_loaded, colours_loaded,
                                 only_load_if_in_edges,
                                 &hdr, sort_kmers, db_graph);

  graph_header_dealloc(&hdr);
  return num_kmers;
}
