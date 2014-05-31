#include "global.h"
#include "graph_format.h"
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"

static inline void dump_empty_bkmer(hkey_t hkey, const dBGraph *db_graph,
                                    char *buf, size_t mem, FILE *fh)
{
  size_t written;
  const BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  written = fwrite(&bkmer, 1, sizeof(BinaryKmer), fh) +
            fwrite(buf, 1, mem, fh);
  if(written != mem+sizeof(BinaryKmer)) die("Couldn't write to file");
}

void graph_write_empty(const dBGraph *db_graph, FILE *fh, size_t num_of_cols)
{
  size_t mem = num_of_cols * (sizeof(Covg)+sizeof(Edges));
  char buf[mem];
  memset(buf, 0, mem);
  HASH_ITERATE(&db_graph->ht, dump_empty_bkmer, db_graph, buf, mem, fh);
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

  uint32_t len = (uint32_t)cleaning->intersection_name.len;
  const char *str = cleaning->intersection_name.buff;
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
      uint32_t len = (uint32_t)h->ginfo[i].sample_name.len;
      const char *buff = h->ginfo[i].sample_name.buff;
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
  if(m != expm) die("Cannot write to file");
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

// Dump all kmers with all colours to given file
void graph_write_all_kmers(FILE *fh, const dBGraph *db_graph)
{
  HASH_ITERATE(&db_graph->ht, graph_write_graph_kmer, fh, db_graph);
}

static inline void overwrite_kmer_colours(hkey_t node,
                                          const dBGraph *db_graph,
                                          Colour graphcol, Colour intocol,
                                          size_t write_ncols, size_t file_ncols,
                                          FILE *fh)
{
  const Edges (*col_edges)[db_graph->num_of_cols]
    = (const Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  const Covg (*col_covgs)[db_graph->num_of_cols]
    = (const Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  size_t skip_cols = file_ncols - (intocol + write_ncols);
  const Covg *covg = col_covgs[node] + graphcol;
  const Edges *edges = col_edges[node] + graphcol;

  long skip0 = sizeof(BinaryKmer) + intocol * sizeof(Covg);
  long skip1 = skip_cols * sizeof(Covg) + intocol * sizeof(Edges);
  long skip2 = skip_cols * sizeof(Edges);

  bool success = (fseek(fh, skip0, SEEK_CUR) == 0) &&
                 (fwrite(covg, sizeof(Covg), write_ncols, fh) == write_ncols) &&
                 (fseek(fh, skip1, SEEK_CUR) == 0) &&
                 (fwrite(edges, sizeof(Edges), write_ncols, fh) == write_ncols) &&
                 (fseek(fh, skip2, SEEK_CUR) == 0);

  if(!success) die("Overwrite failed");
}

void graph_file_write_colours(const dBGraph *db_graph,
                             Colour graphcol, Colour intocol,
                             size_t write_ncols, size_t file_ncols,
                             FILE *fh)
{
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  HASH_ITERATE(&db_graph->ht, overwrite_kmer_colours,
                db_graph, graphcol, intocol, write_ncols, file_ncols, fh);
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

  return (num_of_cols == num_graph_cols && (cols == NULL || start_col == 0));
}

// This function will dump valid binaries by not printing edges to nodes that
// are not themselves printed
// graph info is REQUIRED!
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
uint64_t graph_file_save(const char *path, const dBGraph *db_graph,
                         const GraphFileHeader *header, size_t intocol,
                         const Colour *colours, Colour start_col,
                         size_t num_of_cols)
{
  // Cannot specify both colours array and start_col
  ctx_assert(colours == NULL || start_col == 0);
  ctx_assert(db_graph->col_edges != NULL);
  ctx_assert(db_graph->col_covgs != NULL);
  ctx_assert(num_of_cols > 0);
  ctx_assert(intocol + num_of_cols <= header->num_of_cols);

  size_t i;
  uint64_t num_nodes_dumped = 0;
  FILE *fout;
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

  status("Writing colours %zu-%zu of %zu", intocol, intocol+num_of_cols-1,
         (size_t)header->num_of_cols);

  if(strcmp(path,"-") == 0) fout = stdout;
  else if((fout = fopen(path, "w")) == NULL)
    die("Unable to open graph file to write: %s\n", path);

  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  // Write header
  graph_write_header(fout, header);

  if(saving_graph_as_is) {
    graph_write_all_kmers(fout, db_graph);
  }
  else {
    HASH_ITERATE(&db_graph->ht, graph_write_node,
                 db_graph, fout, header, intocol, colours, start_col, num_of_cols,
                 &num_nodes_dumped);
  }

  fflush(fout);
  fclose(fout);
  // if(strcmp(path,"-") != 0) fclose(fout);

  graph_write_status(num_nodes_dumped, num_of_cols, out_name, header->version);

  return num_nodes_dumped;
}

uint64_t graph_file_save_mkhdr(const char *path, const dBGraph *db_graph,
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
  return graph_file_save(path, db_graph, &header, 0,
                         colours, start_col, num_of_cols);
}

void graph_write_status(uint64_t nkmers, size_t ncols,
                        const char *path, uint32_t version)
{
  char num_kmer_str[100];
  ulong_to_str(nkmers, num_kmer_str);

  status("Dumped %s kmers in %zu colour%s into: %s (format version: %u)\n",
         num_kmer_str, ncols, util_plural_str(ncols),
         futil_outpath_str(path), version);
}
