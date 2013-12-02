#include "global.h"
#include "graph_format.h"
#include "db_graph.h"
#include "db_node.h"
#include "util.h"

static inline void dump_empty_bkmer(hkey_t node, const dBGraph *db_graph,
                                    char *buf, size_t mem, FILE *fh)
{
  const BinaryKmer bkmer = db_node_bkmer(db_graph, node);
  fwrite(&bkmer, sizeof(BinaryKmer), 1, fh);
  fwrite(buf, 1, mem, fh);
}

void graph_write_empty(const dBGraph *db_graph, FILE *fh, size_t num_of_cols)
{
  size_t mem = num_of_cols * (sizeof(Covg)+sizeof(Edges));
  char buf[mem];
  memset(buf, 0, mem);
  HASH_TRAVERSE(&db_graph->ht, dump_empty_bkmer, db_graph, buf, mem, fh);
}

// Returns number of bytes written
static size_t write_error_cleaning_object(FILE *fh, const ErrorCleaning *cleaning)
{
  fwrite(&(cleaning->tip_clipping), sizeof(uint8_t), 1, fh);
  fwrite(&(cleaning->remv_low_cov_sups), sizeof(uint8_t), 1, fh);
  fwrite(&(cleaning->remv_low_cov_nodes), sizeof(uint8_t), 1, fh);
  fwrite(&(cleaning->is_graph_intersection), sizeof(uint8_t), 1, fh);

  uint32_t supernodes_cleaning_thresh
    = cleaning->remv_low_cov_sups ? cleaning->remv_low_cov_sups_thresh : 0;

  uint32_t nodes_cleaning_thresh
    = cleaning->remv_low_cov_nodes ? cleaning->remv_low_cov_nodes_thresh : 0;

  fwrite(&supernodes_cleaning_thresh, sizeof(uint32_t), 1, fh);
  fwrite(&nodes_cleaning_thresh, sizeof(uint32_t), 1, fh);

  uint32_t len = cleaning->intersection_name.len;
  char *str = cleaning->intersection_name.buff;
  fwrite(&len, sizeof(uint32_t), 1, fh);
  fwrite(str, sizeof(uint8_t), len, fh);

  return 4 + sizeof(uint32_t) * 2 + sizeof(uint32_t) + len;
}

// Returns number of bytes written
size_t graph_write_header(FILE *fh, const GraphFileHeader *h)
{
  uint32_t i;
  size_t b = 0;

  fwrite("CORTEX", sizeof(char), strlen("CORTEX"), fh);
  fwrite(&h->version, sizeof(uint32_t), 1, fh);
  fwrite(&h->kmer_size, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_bitfields, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_cols, sizeof(uint32_t), 1, fh);

  b += strlen("CORTEX") + sizeof(uint32_t) * 4;

  if(h->version >= 7)
  {
    uint32_t num_of_shades = 0;
    fwrite(&(h->num_of_kmers), sizeof(uint64_t), 1, fh);
    fwrite(&num_of_shades, sizeof(uint32_t), 1, fh);
    b += sizeof(uint64_t) + sizeof(uint32_t);
  }

  for(i = 0; i < h->num_of_cols; i++)
    fwrite(&h->ginfo[i].mean_read_length, sizeof(uint32_t), 1, fh);
  for(i = 0; i < h->num_of_cols; i++)
    fwrite(&h->ginfo[i].total_sequence, sizeof(uint64_t), 1, fh);

  b += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint64_t));

  if(h->version >= 6)
  {
    for(i = 0; i < h->num_of_cols; i++)
    {
      uint32_t len = h->ginfo[i].sample_name.len;
      char *buff = h->ginfo[i].sample_name.buff;
      fwrite(&len, sizeof(uint32_t), 1, fh);
      fwrite(buff, sizeof(uint8_t), len, fh);
      b += sizeof(uint32_t) + len;
    }

    for(i = 0; i < h->num_of_cols; i++)
      fwrite(&h->ginfo[i].seq_err, sizeof(long double), 1, fh);

    b += h->num_of_cols * sizeof(long double);

    for(i = 0; i < h->num_of_cols; i++)
      b += write_error_cleaning_object(fh, &h->ginfo[i].cleaning);
  }

  fwrite("CORTEX", sizeof(uint8_t), strlen("CORTEX"), fh);
  b += strlen("CORTEX");

  return b;
}

// Returns number of bytes written
size_t graph_write_kmer(FILE *fh, const GraphFileHeader *h,
                        const uint64_t *bkmer, const Covg *covgs,
                        const Edges *edges)
{
  fwrite(bkmer, sizeof(uint64_t), h->num_of_bitfields, fh);
  fwrite(covgs, sizeof(uint32_t), h->num_of_cols, fh);
  fwrite(edges, sizeof(uint8_t),  h->num_of_cols, fh);
  return 8*h->num_of_bitfields + 5*h->num_of_cols;
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

  fseek(fh, sizeof(BinaryKmer) + intocol * sizeof(Covg), SEEK_CUR);
  fwrite(covg, sizeof(Covg), write_ncols, fh);
  fseek(fh, skip_cols * sizeof(Covg) + intocol * sizeof(Edges), SEEK_CUR);
  fwrite(edges, sizeof(Edges), write_ncols, fh);
  fseek(fh, skip_cols * sizeof(Edges), SEEK_CUR);
}

// DEV: remove this and graph_file_write_colour (+from header)
static inline void overwrite_kmer_colour(hkey_t node, const dBGraph *db_graph,
                                         Colour graphcol, Colour intocol,
                                         size_t file_ncols, FILE *fh)
{
  overwrite_kmer_colours(node, db_graph, graphcol, intocol, 1, file_ncols, fh);
}

// Dump a single colour into an existing binary
// FILE *fh must already point to the first bkmer
// if merge is true, read existing covg and edges and combine with outgoing
void graph_file_write_colour(const dBGraph *db_graph, Colour graphcol,
                             Colour intocol, size_t file_ncols, FILE *fh)
{
  HASH_TRAVERSE(&db_graph->ht, overwrite_kmer_colour,
                db_graph, graphcol, intocol, file_ncols, fh);
}

void graph_file_write_colours(const dBGraph *db_graph,
                             Colour graphcol, Colour intocol,
                             size_t write_ncols, size_t file_ncols,
                             FILE *fh)
{
  assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  HASH_TRAVERSE(&db_graph->ht, overwrite_kmer_colours,
                db_graph, graphcol, intocol, write_ncols, file_ncols, fh);
}

// Dump node: only print kmers with coverages in given colours
static void graph_write_node(hkey_t node, const dBGraph *db_graph,
                             FILE *fout, const GraphFileHeader *header,
                             size_t intocol, const Colour *colours,
                             size_t start_col, size_t num_of_cols,
                             uint64_t *num_dumped)
{
  size_t i = 0;

  // Check this node has coverage in one of the specified colours
  if(colours != NULL)
    while(i < num_of_cols && db_node_get_covg(db_graph,node,colours[i]) == 0) i++;
  else
    while(i < num_of_cols && db_node_get_covg(db_graph,node,start_col+i) == 0) i++;

  if(i == num_of_cols) return;

  BinaryKmer bkmer = db_node_bkmer(db_graph, node);
  Covg covg_store[header->num_of_cols], *covgs = covg_store + intocol;
  Edges edge_store[header->num_of_cols], *edges = edge_store + intocol;

  memset(covg_store, 0, sizeof(Covg) * header->num_of_cols);
  memset(edge_store, 0, sizeof(Edges) * header->num_of_cols);

  Edges (*col_edges)[db_graph->num_of_cols]
    = (Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  Covg (*col_covgs)[db_graph->num_of_cols]
    = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  if(colours != NULL) {
    for(i = 0; i < num_of_cols; i++) {
      covgs[i] = col_covgs[node][colours[i]];
      edges[i] = col_edges[node][colours[i]];
    }
  }
  else {
    memcpy(covgs, col_covgs[node]+start_col, num_of_cols*sizeof(Covg));
    memcpy(edges, col_edges[node]+start_col, num_of_cols*sizeof(Edges));
  }

  graph_write_kmer(fout, header, bkmer.b, covg_store, edge_store);

  (*num_dumped)++;
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
  assert(colours == NULL || start_col == 0);
  assert(db_graph->col_edges != NULL);
  assert(db_graph->col_covgs != NULL);
  assert(num_of_cols > 0);

  size_t i;
  uint64_t num_nodes_dumped = 0;
  FILE *fout;

  if(colours != NULL) {
    if(num_of_cols == 1)
      status("Dumping graph colour %zu into: %s", colours[0], path);
    else {
      timestamp(ctx_msg_out);
      message("Dumping graph colours %zu", colours[0]);
      for(i = 1; i < num_of_cols; i++) message(",%zu", colours[i]);
      message(" into: %s\n", path);
    }
  }
  else if(num_of_cols == 1)
    status("Dumping graph colour %zu into: %s", start_col, path);
  else {
    status("Dumping graph colours %zu-%zu into: %s", start_col,
           start_col+num_of_cols-1, path);
  }

  if((fout = fopen(path, "w")) == NULL)
    die("Unable to open dump binary file to write: %s\n", path);

  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  // Write header
  graph_write_header(fout, header);

  HASH_TRAVERSE(&db_graph->ht, graph_write_node, db_graph, fout, header,
                intocol, colours, start_col, num_of_cols, &num_nodes_dumped);

  if(header->version >= 7 && num_nodes_dumped != db_graph->ht.unique_kmers)
  {
    // Need to go back and update number of kmers dumped
    long pos = strlen("CORTEX") * sizeof(uint8_t) + 4 * sizeof(uint32_t);

    if(fseek(fout, pos, SEEK_SET) == 0) {
      fwrite(&num_nodes_dumped, sizeof(uint64_t), 1, fout);
    } else {
      warn("Couldn't update number of kmers in file [binary: %s]", path);
    }
  }

  fclose(fout);

  graph_write_status(num_nodes_dumped, num_of_cols, path, header->version);

  return num_nodes_dumped;
}

uint64_t graph_file_save_mkhdr(const char *path, const dBGraph *db_graph,
                         uint32_t version,
                         const Colour *colours, Colour start_col,
                         uint32_t num_of_cols)
{
  // Construct binary header
  GraphInfo hdr_ginfo[num_of_cols];
  GraphFileHeader header = {.version = version,
                            .kmer_size = db_graph->kmer_size,
                            .num_of_bitfields = NUM_BKMER_WORDS,
                            .num_of_cols = num_of_cols,
                            .num_of_kmers = db_graph->ht.unique_kmers,
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
         num_kmer_str, ncols, ncols != 1 ? "s" : "",
         path, version);
}
