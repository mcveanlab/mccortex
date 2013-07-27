#include "global.h"
#include "binary_format.h"
#include "db_graph.h"
#include "db_node.h"

const int CURR_CTX_VERSION = 6;
const char CTX_MAGIC_WORD[7] = "CORTEX";

static inline void dump_empty_bkmer(hkey_t node, dBGraph *db_graph,
                                    char *buf, size_t mem, FILE *fh)
{
  // printf("dump_empty_bkmer\n");
  fwrite(db_node_bkmer(db_graph, node), sizeof(BinaryKmer), 1, fh);
  fwrite(buf, 1, mem, fh);

  char bkmerstr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, bkmerstr);
}

void dump_empty_binary(dBGraph *db_graph, FILE *fh, uint32_t num_of_cols)
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
  fwrite(&(cleaning->cleaned_against_another_graph), sizeof(uint8_t), 1, fh);

  uint32_t supernodes_cleaning_thresh
    = cleaning->remv_low_cov_sups ? cleaning->remv_low_cov_sups_thresh : 0;

  uint32_t nodes_cleaning_thresh
    = cleaning->remv_low_cov_nodes ? cleaning->remv_low_cov_nodes_thresh : 0;

  fwrite(&supernodes_cleaning_thresh, sizeof(uint32_t), 1, fh);
  fwrite(&nodes_cleaning_thresh, sizeof(uint32_t), 1, fh);

  uint32_t len = cleaning->cleaned_against_graph_name.len;
  char *str = cleaning->cleaned_against_graph_name.buff;
  fwrite(&len, sizeof(uint32_t), 1, fh);
  fwrite(str, sizeof(uint8_t), len, fh);

  return 4 + sizeof(uint32_t) * 2 + sizeof(uint32_t) + len;
}

// Returns number of bytes written
size_t binary_write_header(FILE *fh, const BinaryFileHeader *h)
{
  uint32_t i;
  size_t b = 0;

  fwrite(CTX_MAGIC_WORD, sizeof(char), strlen(CTX_MAGIC_WORD), fh);
  fwrite(&h->version, sizeof(uint32_t), 1, fh);
  fwrite(&h->kmer_size, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_bitfields, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_cols, sizeof(uint32_t), 1, fh);

  b += strlen(CTX_MAGIC_WORD) + sizeof(uint32_t) * 4;

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

  fwrite(CTX_MAGIC_WORD, sizeof(uint8_t), strlen(CTX_MAGIC_WORD), fh);
  b += strlen(CTX_MAGIC_WORD);

  return b;
}

// Returns number of bytes written
size_t binary_write_kmer(FILE *fh, const BinaryFileHeader *h,
                         const uint64_t *bkmer, const Covg *covgs,
                         const Edges *edges)
{
  fwrite(bkmer, sizeof(uint64_t), h->num_of_bitfields, fh);
  fwrite(covgs, sizeof(uint32_t), h->num_of_cols, fh);
  fwrite(edges, sizeof(uint8_t),  h->num_of_cols, fh);
  return 8*h->num_of_bitfields + 5*h->num_of_cols;
}

static inline void overwrite_kmer_colour(hkey_t node, dBGraph *db_graph,
                                         Colour graphcol, Colour intocol,
                                         uint32_t num_of_cols, FILE *fh)
{
  Edges (*col_edges)[db_graph->num_of_cols]
    = (Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  Covg (*col_covgs)[db_graph->num_of_cols]
    = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  uint32_t skip_cols = num_of_cols - intocol - 1;
  Covg covg = col_covgs[node][graphcol];
  Edges edges = col_edges[node][graphcol];

  fseek(fh, sizeof(BinaryKmer) + intocol * sizeof(Covg), SEEK_CUR);
  fwrite(&covg, sizeof(Covg), 1, fh);
  fseek(fh, skip_cols * sizeof(Covg), SEEK_CUR);

  fseek(fh, intocol * sizeof(Edges), SEEK_CUR);
  fwrite(&edges, sizeof(Edges), 1, fh);
  fseek(fh, skip_cols * sizeof(Edges), SEEK_CUR);
}

// Dump a single colour into an existing binary
// FILE *fh must already point to the first bkmer
// if merge is true, read existing covg and edges and combine with outgoing
void binary_dump_colour(dBGraph *db_graph, Colour graphcol,
                        Colour intocol, uint32_t num_of_cols, FILE *fh)
{
  HASH_TRAVERSE(&db_graph->ht, overwrite_kmer_colour,
                db_graph, graphcol, intocol, num_of_cols, fh);
}

// Dump node: only print kmers with coverages in given colours
static void binary_dump_node_colours(hkey_t node, const dBGraph *db_graph,
                                     FILE *fout, const BinaryFileHeader *header,
                                     const Colour *colours, uint32_t start_col,
                                     uint64_t *num_of_nodes_dumped)
{
  uint32_t i = 0;

  // Check this node has coverage in one of the specified colours
  if(colours != NULL) {
    while(i < header->num_of_cols &&
          db_node_get_covg(db_graph,node,colours[i]) == 0) i++;
  }
  else {
    while(i < header->num_of_cols &&
          db_node_get_covg(db_graph,node,start_col+i) == 0) i++;
  }
  if(i == header->num_of_cols) return;

  ConstBinaryKmerPtr bkmer = db_node_bkmer(db_graph, node);
  Covg covgs[header->num_of_cols];
  Edges edges[header->num_of_cols];

  Edges (*col_edges)[db_graph->num_of_cols]
    = (Edges (*)[db_graph->num_of_cols])db_graph->col_edges;
  Covg (*col_covgs)[db_graph->num_of_cols]
    = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  if(colours != NULL)
  {
    uint32_t col;
    for(i = 0; i < header->num_of_cols; i++)
    {
      col = colours[i];
      covgs[i] = col_covgs[node][col];
      edges[i] = col_edges[node][col];
    }
  }
  else
  {
    memcpy(covgs, col_covgs[node]+start_col, header->num_of_cols*sizeof(Covg));
    memcpy(edges, col_edges[node]+start_col, header->num_of_cols*sizeof(Edges));
  }

  binary_write_kmer(fout, header, bkmer, covgs, edges);

  (*num_of_nodes_dumped)++;
}

// This function will dump valid binaries by not printing edges to nodes that
// are not themselves printed
// graph info is REQUIRED!
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
uint64_t binary_dump_graph(const char *path, dBGraph *db_graph,
                           uint32_t version,
                           const Colour *colours, Colour start_col,
                           uint32_t num_of_cols)
{
  assert(db_graph->col_edges != NULL);
  assert(db_graph->col_covgs != NULL);

  FILE *fout = fopen(path, "w");

  if(fout == NULL) {
    die("Unable to open dump binary file to write: %s\n", path);
  }

  GraphInfo *ginfo = db_graph->ginfo;

  // Construct binary header
  BinaryFileHeader header = {.version = version,
                             .kmer_size = db_graph->kmer_size,
                             .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                             .num_of_cols = num_of_cols,
                             .num_of_kmers = db_graph->ht.unique_kmers};

  GraphInfo header_ginfo[num_of_cols];

  uint32_t i, col;
  for(i = 0; i < num_of_cols; i++) {
    col = colours != NULL ? colours[i] : i;
    header_ginfo[i] = ginfo[col];
  }

  header.ginfo = header_ginfo;

  // Write header
  binary_write_header(fout, &header);

  uint64_t num_of_nodes_dumped = 0;

  HASH_TRAVERSE(&(db_graph->ht), binary_dump_node_colours,
                db_graph, fout, &header,
                colours, start_col, &num_of_nodes_dumped);

  if(version >= 7 && num_of_nodes_dumped != db_graph->ht.unique_kmers)
  {
    // Need to go back and update number of kmers dumped
    long pos = strlen(CTX_MAGIC_WORD) * sizeof(uint8_t) + 4 * sizeof(uint32_t);

    if(fseek(fout, pos, SEEK_SET) == 0) {
      fwrite(&num_of_nodes_dumped, sizeof(uint64_t), 1, fout);
    } else {
      warn("Couldn't update number of kmers in file [binary: %s]", path);
    }
  }

  fclose(fout);

  return num_of_nodes_dumped;
}
