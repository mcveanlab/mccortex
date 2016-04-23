#include "global.h"
#include "graph_writer.h"
#include "graphs_load.h" // need to load, merge then write some graphs
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"

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
size_t graph_write_empty(const dBGraph *db_graph, FILE *fh, size_t num_of_cols)
{
  size_t mem = num_of_cols * (sizeof(Covg)+sizeof(Edges));
  char buf[mem];
  memset(buf, 0, mem);
  HASH_ITERATE(&db_graph->ht, _dump_empty_bkmer, db_graph, buf, mem, fh);
  return db_graph->ht.num_kmers * (sizeof(BinaryKmer) + mem);
}

// Returns number of bytes written
static size_t write_error_cleaning_object(FILE *fh, const ErrorCleaning *cleaning)
{
  size_t written = 0, expwrite;
  written += fwrite(&(cleaning->cleaned_tips), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->cleaned_snodes), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->cleaned_kmers), 1, sizeof(uint8_t), fh);
  written += fwrite(&(cleaning->is_graph_intersection), 1, sizeof(uint8_t), fh);

  uint32_t clean_snodes_thresh
    = cleaning->cleaned_snodes ? cleaning->clean_snodes_thresh : 0;

  uint32_t clean_kmers_thresh
    = cleaning->cleaned_kmers ? cleaning->clean_kmers_thresh : 0;

  written += fwrite(&clean_snodes_thresh, 1, sizeof(uint32_t), fh);
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

// Returns number of bytes written
size_t graph_write_kmer(FILE *fh, size_t num_bkmer_words, size_t num_cols,
                        const BinaryKmer bkmer, const Covg *covgs,
                        const Edges *edges)
{
  size_t m = 0, expm = 8*num_bkmer_words + 5*num_cols;
  m += fwrite(bkmer.b, 1, sizeof(uint64_t) * num_bkmer_words, fh);
  m += fwrite(covgs, 1, sizeof(uint32_t) * num_cols, fh);
  m += fwrite(edges, 1, sizeof(uint8_t) * num_cols, fh);
  if(m != expm) die("Cannot write to file (%zu, %zu)", m, expm);
  return m;
}

static inline void graph_write_graph_kmer(hkey_t hkey, FILE *fh,
                                          const dBGraph *db_graph)
{
  graph_write_kmer(fh, NUM_BKMER_WORDS, db_graph->num_of_cols,
                   db_graph->ht.table[hkey],
                   &db_node_covg(db_graph, hkey, 0),
                   &db_node_edges(db_graph, hkey, 0));
}

// Dump all kmers with all colours to given file. Return num of kmers written
size_t graph_write_all_kmers(FILE *fh, const dBGraph *db_graph)
{
  HASH_ITERATE(&db_graph->ht, graph_write_graph_kmer, fh, db_graph);
  return db_graph->ht.num_kmers;
}


// only called by graph_writer_update_mmap_kmers()
static inline void _graph_write_update_kmer(hkey_t hkey,
                                           const dBGraph *db_graph,
                                           size_t first_graphcol, size_t ngraphcols,
                                           size_t first_filecol, size_t nfilecols,
                                           char **ptr, size_t filekmersize)
{
  const Covg *covgs = &db_node_covg(db_graph, hkey, first_graphcol);
  const Edges *edges = &db_node_edges(db_graph, hkey, first_graphcol);

  void *covgs_out = *ptr+sizeof(BinaryKmer) + sizeof(Covg)*first_filecol;
  void *edges_out = *ptr+sizeof(BinaryKmer) + sizeof(Covg)*nfilecols +
                     sizeof(Edges)*first_filecol;

  memcpy(covgs_out, covgs, ngraphcols*sizeof(Covg));
  memcpy(edges_out, edges, ngraphcols*sizeof(Edges));

  *ptr += filekmersize;
}

/*!
  Overwrite kmers in an existing file.
  @param first_graphcol first colour in the dBGraph to read from
  @param first_filecol first colour in the file to write into
  @param ngraphcols Number of colours to write to file
  @param nfilecols Total number of colours in file
  @param mmap_ptr Memory mapped file pointer
  @param hdrsize Size of file header i.e. byte pos of first kmer in file
 */
void graph_writer_update_mmap_kmers(const dBGraph *db_graph,
                                    size_t first_graphcol, size_t ngraphcols,
                                    size_t first_filecol, size_t nfilecols,
                                    char *mmap_ptr, size_t hdrsize)
{
  ctx_assert(db_graph->col_edges != NULL);
  ctx_assert(db_graph->col_covgs != NULL);
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  ctx_assert(first_graphcol+ngraphcols <= db_graph->num_of_cols);

  char *ptr = mmap_ptr + hdrsize;
  size_t filekmersize = sizeof(BinaryKmer)+(sizeof(Edges)+sizeof(Covg))*nfilecols;

  HASH_ITERATE(&db_graph->ht, _graph_write_update_kmer,
               db_graph, first_graphcol, ngraphcols, first_filecol, nfilecols,
               &ptr, filekmersize);
}

/*!
  Overwrite kmers in an existing file.
  @param first_graphcol  first colour in the dBGraph to read from
  @param first_filecol   first colour in the file to write into
  @param ngraphcols      number of colours to write to file
  @param nfilecols       total number of colours in file
  @param hdrsize         file header size i.e. byte pos of first kmer in file
  @param fh              file handle opened that we can fwrite/fseek
  @param path            path to the file (for error messages)
 */
void graph_writer_update_file_kmers(const dBGraph *db_graph,
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

  if(fseek(fh, hdrsize, SEEK_SET) != 0) die("Cannot seek to file start: %s", path);

  // Read block of kmers from the file
  // iterate over the hash table figuring which ones they are
  while(nkmers_printed < db_graph->ht.num_kmers)
  {
    nkmers = MIN2(db_graph->ht.num_kmers - nkmers_printed, kmers_per_block);
    nbytes = nkmers*filekmersize;
    if(fread(mem, 1, nbytes, fh) != nbytes) die("Cannot read: %s", path);
    memptr = mem;
    for(end = nkmers_printed+nkmers; nkmers_printed < end; hkey++) {
      if(db_graph_node_assigned(db_graph, hkey)) {
        covgs = db_graph->col_covgs + hkey * db_graph->num_of_cols;
        edges = db_graph->col_edges + hkey * db_graph->num_of_cols;
        memptr += sizeof(BinaryKmer);
        memcpy(memptr + first_filecol*sizeof(Covg), covgs, ngraphcols*sizeof(Covg));
        memptr += sizeof(Covg)*nfilecols;
        memcpy(memptr + first_filecol*sizeof(Edges), edges, ngraphcols*sizeof(Edges));
        memptr += sizeof(Edges)*nfilecols;
        nkmers_printed++;
      }
    }
    if(fseek(fh, -nbytes, SEEK_CUR) != 0) die("fseek failed: %s", path);
    if(fwrite(mem, 1, nbytes, fh) != nbytes) die("fwrite failed: %s", path);
  }

  ctx_free(mem);
}

// Dump node: only print kmers with coverages in given colours
static void graph_write_node(hkey_t hkey, const dBGraph *db_graph,
                             FILE *fout, const GraphFileHeader *hdr,
                             size_t intocol, const Colour *colours,
                             size_t start_col, size_t num_of_cols,
                             uint64_t *num_dumped)
{
  ctx_assert(num_of_cols > 0);
  ctx_assert(intocol+num_of_cols <= hdr->num_of_cols);
  size_t i = 0;

  // Check this node has coverage in one of the specified colours
  if(colours != NULL)
    while(i < num_of_cols && db_node_get_covg(db_graph,hkey,colours[i]) == 0) i++;
  else
    while(i < num_of_cols && db_node_get_covg(db_graph,hkey,start_col+i) == 0) i++;

  if(i == num_of_cols) return;

  BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  Covg covg_store[hdr->num_of_cols], *covgs = covg_store + intocol;
  Edges edge_store[hdr->num_of_cols], *edges = edge_store + intocol;

  memset(covg_store, 0, sizeof(Covg) * hdr->num_of_cols);
  memset(edge_store, 0, sizeof(Edges) * hdr->num_of_cols);

  Edges (*col_edges)[db_graph->num_of_cols]
    = (Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  Covg (*col_covgs)[db_graph->num_of_cols]
    = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  if(colours != NULL) {
    for(i = 0; i < num_of_cols; i++) {
      covgs[i] = col_covgs[hkey][colours[i]];
      edges[i] = col_edges[hkey][colours[i]];
    }
  }
  else {
    memcpy(covgs, col_covgs[hkey]+start_col, num_of_cols*sizeof(Covg));
    memcpy(edges, col_edges[hkey]+start_col, num_of_cols*sizeof(Edges));
  }

  graph_write_kmer(fout, hdr->num_of_bitfields, hdr->num_of_cols,
                   bkmer, covg_store, edge_store);

  (*num_dumped)++;
}

// Returns true if we are dumping the graph 'as-is', without dropping or
// re-arranging colours
static bool saving_graph_as_is(const Colour *cols, Colour start_col,
                               size_t num_of_cols, size_t num_graph_cols)
{
  size_t i;

  if(cols != NULL) {
    for(i = 0; i < num_of_cols; i++)
      if(cols[i] != i)
        return false;
  }
  else if (start_col != 0)
    return false;

  return (num_of_cols == num_graph_cols);
}

// start_col is ignored unless colours is NULL
uint64_t graph_writer_save(const char *path, const dBGraph *db_graph,
                           const GraphFileHeader *header, size_t intocol,
                           const Colour *colours, Colour start_col,
                           size_t num_of_cols)
{
  // Cannot specify both colours array and start_col
  ctx_assert(colours == NULL || start_col == 0);
  ctx_assert(db_graph->col_edges != NULL);
  ctx_assert(db_graph->col_covgs != NULL);
  ctx_assert(num_of_cols > 0);
  ctx_assert(colours || start_col + num_of_cols <= db_graph->num_of_cols);
  ctx_assert(intocol + num_of_cols <= header->num_of_cols);

  size_t i;
  uint64_t num_nodes_dumped = 0;
  const char *out_name = futil_outpath_str(path);

  if(colours != NULL) {
    if(num_of_cols == 1)
      status("Dumping graph colour %zu into: %s", colours[0], out_name);
    else {
      timestamp();
      message("Dumping graph colours %zu", colours[0]);
      for(i = 1; i < num_of_cols; i++) message(",%zu", colours[i]);
      message(" into: %s\n", out_name);
    }
  }
  else if(num_of_cols == 1)
    status("Dumping graph colour %zu into: %s", start_col, out_name);
  else {
    status("Dumping graph colours %zu-%zu into: %s", start_col,
           start_col+num_of_cols-1, out_name);
  }

  status("[graph_writer_save] Writing colours %zu-%zu of %zu into: %s",
         intocol, intocol+num_of_cols-1, (size_t)header->num_of_cols,
         futil_outpath_str(path));

  FILE *fout = futil_fopen(path, "w");

  // Write header
  graph_write_header(fout, header);

  if(saving_graph_as_is(colours, start_col, num_of_cols, db_graph->num_of_cols)) {
    num_nodes_dumped = graph_write_all_kmers(fout, db_graph);
  }
  else {
    HASH_ITERATE(&db_graph->ht, graph_write_node,
                 db_graph, fout, header, intocol, colours, start_col, num_of_cols,
                 &num_nodes_dumped);
  }

  fclose(fout);
  // if(strcmp(path,"-") != 0) fclose(fout);

  graph_writer_print_status(num_nodes_dumped, num_of_cols,
                            out_name, header->version);

  return num_nodes_dumped;
}

uint64_t graph_writer_save_mkhdr(const char *path, const dBGraph *db_graph,
                                 uint32_t version,
                                 const Colour *colours, Colour start_col,
                                 size_t num_of_cols)
{
  // Construct graph header
  GraphInfo hdr_ginfo[num_of_cols];
  GraphFileHeader header = {.version = version,
                            .kmer_size = (uint32_t)db_graph->kmer_size,
                            .num_of_bitfields = NUM_BKMER_WORDS,
                            .num_of_cols = (uint32_t)num_of_cols,
                            .capacity = 0};

  size_t i;
  GraphInfo *ginfo = db_graph->ginfo;
  for(i = 0; i < num_of_cols; i++)
    hdr_ginfo[i] = ginfo[colours != NULL ? colours[i] : i];

  header.ginfo = hdr_ginfo;
  return graph_writer_save(path, db_graph, &header, 0,
                           colours, start_col, num_of_cols);
}

void graph_writer_print_status(uint64_t nkmers, size_t ncols,
                               const char *path, uint32_t version)
{
  char num_kmer_str[100];
  ulong_to_str(nkmers, num_kmer_str);

  status("Dumped %s kmers in %zu colour%s into: %s (format version: %u)\n",
         num_kmer_str, ncols, util_plural_str(ncols),
         futil_outpath_str(path), version);
}

//
// Merging, filtering, combining graph files
//

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
    return graph_writer_save(out_ctx_path, db_graph, hdr,
                             0, NULL, 0, output_colours);
  }
  else if(num_files == 1)
  {
    return graph_writer_stream(out_ctx_path, &files[0], db_graph, hdr,
                               only_load_if_in_edges);
  }

  GraphLoadingPrefs gprefs = graph_loading_prefs(db_graph);
  gprefs.must_exist_in_graph = only_load_if_in_graph;
  gprefs.must_exist_in_edges = only_load_if_in_edges;

  if(output_colours <= db_graph->num_of_cols)
  {
    // Can load all files at once
    status("Loading and saving %zu colours at once", output_colours);

    if(!kmers_loaded) {
      for(i = 0; i < num_files; i++)
        graph_load(&files[i], gprefs, NULL);
      hash_table_print_stats(&db_graph->ht);
    }

    graph_writer_save(out_ctx_path, db_graph, hdr, 0, NULL, 0, output_colours);
  }
  else
  {
    ctx_assert2(strcmp(out_ctx_path,"-") != 0,
                "Cannot use STDOUT for output if not enough colours to load");

    // Have to load a few colours at a time then dump, rinse and repeat
    status("[overwriting] Saving %zu colours, %zu colours at a time",
           output_colours, db_graph->num_of_cols);

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
    file_len += graph_write_empty(db_graph, fout, output_colours);
    fflush(fout);

    // Open memory mapped file
    // void *mmap_ptr = mmap(NULL, file_len, PROT_WRITE, MAP_SHARED, fileno(fout), 0);

    // if(mmap_ptr == MAP_FAILED)
    //   die("Cannot memory map file: %s [%s]", out_ctx_path, strerror(errno));

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
        file_filter_update(fltr);

        if(file_filter_num(fltr) > 0)
          files_loaded |= (graph_load(&files[f], gprefs, NULL) > 0);

        // Restore original filter
        file_filter_copy(fltr, &origfltr);
      }

      // if files_loaded, dump
      if(files_loaded) {
        if(db_graph->num_of_cols == 1)
          status("Dumping into colour %zu...\n", firstcol);
        else
          status("Dumping into colours %zu-%zu...\n", firstcol, lastcol);

        // graph_writer_update_mmap_kmers(db_graph, 0, lastcol-firstcol+1,
        //                                firstcol, output_colours,
        //                                mmap_ptr, hdr_size);

        graph_writer_update_file_kmers(db_graph, 0, lastcol-firstcol+1,
                                       firstcol, output_colours,
                                       hdr_size, fout, out_ctx_path);
      }
    }

    // if(munmap(mmap_ptr, file_len) == -1)
    //   die("Cannot release mmap file: %s [%s]", out_ctx_path, strerror(errno));

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
size_t graph_writer_merge_mkhdr(const char *out_ctx_path,
                               GraphFileReader *files, size_t num_files,
                               bool kmers_loaded, bool colours_loaded,
                               const Edges *only_load_if_in_edges,
                               const char *intersect_gname, dBGraph *db_graph)
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
                                &hdr, db_graph);

  graph_header_dealloc(&hdr);
  return num_kmers;
}
