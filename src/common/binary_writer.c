#include "global.h"
#include "binary_format.h"
#include "db_graph.h"

const int CURR_CTX_VERSION = 7;
const char CTX_MAGIC_WORD[7] = "CORTEX";

static void write_error_cleaning_object(FILE *fh, const ErrorCleaning *cleaning)
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

  uint32_t len = cleaning->cleaned_against_graph_name->len;
  char *str = cleaning->cleaned_against_graph_name->buff;
  fwrite(&len, sizeof(uint32_t), 1, fh);
  fwrite(str, sizeof(uint8_t), len, fh);
}

void binary_write_header(FILE *fh, const BinaryFileHeader *h)
{
  fwrite(CTX_MAGIC_WORD, sizeof(char), strlen(CTX_MAGIC_WORD), fh);
  fwrite(&h->version, sizeof(uint32_t), 1, fh);
  fwrite(&h->kmer_size, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_bitfields, sizeof(uint32_t), 1, fh);
  fwrite(&h->num_of_colours, sizeof(uint32_t), 1, fh);

  if(h->version >= 7)
  {
    fwrite(&(h->num_of_kmers), sizeof(uint64_t), 1, fh);
    fwrite(&(h->num_of_shades), sizeof(uint32_t), 1, fh);
  }

  fwrite(h->mean_read_lengths, sizeof(uint32_t), h->num_of_colours, fh);
  fwrite(h->total_seq_loaded, sizeof(uint64_t), h->num_of_colours, fh);

  if(h->version >= 6)
  {
    uint32_t i;

    for(i = 0; i < h->num_of_colours; i++)
    {
      uint32_t len = h->sample_names[i]->len;
      char *buff = h->sample_names[i]->buff;
      fwrite(&len, sizeof(uint32_t), 1, fh);
      fwrite(buff, sizeof(uint8_t), len, fh);
    }

    fwrite(h->seq_err_rates, sizeof(long double), h->num_of_colours, fh);

    for(i = 0; i < h->num_of_colours; i++) {
      write_error_cleaning_object(fh, h->err_cleaning[i]);
    }
  }

  fwrite(CTX_MAGIC_WORD, sizeof(uint8_t), strlen(CTX_MAGIC_WORD), fh);
}

void binary_write_kmer(FILE *fh, const BinaryFileHeader *h,
                       const Colour *colours, Colour start_col,
                       const uint64_t *bkmer, const Covg *covgs,
                       const Edges *edges)
{
  fwrite(bkmer, sizeof(uint64_t), h->num_of_bitfields, fh);

  if(colours == NULL) {
    fwrite(covgs+start_col, sizeof(uint32_t), h->num_of_colours, fh);
  } else {
    size_t i;
    for(i = 0; i < h->num_of_colours; i++)
      fwrite(covgs+colours[i], sizeof(uint32_t), 1, fh);
  }

  fwrite(edges, sizeof(uint8_t),  h->num_of_colours, fh);
}

// Dev: what if dumping a single colour -- do we dump nodes with no covg?
// Dump node: only print colour coverages for given colours
static void binary_dump_node_colours(hkey_t node, const dBGraph *db_graph,
                                     FILE *fout, const BinaryFileHeader *header,
                                     const Colour *colours, uint32_t start_col,
                                     uint64_t *num_of_nodes_dumped)
{
  uint32_t i = 0;

  if(header->num_of_colours < NUM_OF_COLOURS)
  {
    // Check this node has coverage in one of the specified colours
    if(colours != NULL) {
      while(i < header->num_of_colours && db_graph->covgs[colours[i]] == 0) i++;
    }
    else {
      while(i < header->num_of_colours && db_graph->covgs[start_col+i] == 0) i++;
    }
    if(i == header->num_of_colours) return;
  }

  BinaryKmerPtr bkmer = db_graph_bkmer(db_graph, node);

  if(colours == NULL)
  {
    binary_write_kmer(fout, header, colours, start_col, bkmer,
                      db_graph->covgs[node], db_graph->col_edges[node]);
  }
  else
  {
    Covg covgs[header->num_of_colours];
    Edges edges[header->num_of_colours];
    uint32_t col;

    for(i = 0; i < header->num_of_colours; i++)
    {
      col = colours[i];
      covgs[i] = db_graph->covgs[node][col];
      edges[i] = db_graph->col_edges[node][col];
    }

    binary_write_kmer(fout, header, colours, start_col,
                      bkmer, covgs, edges);
  }

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
                           uint32_t num_of_cols, uint32_t num_of_shades)
{
  assert(db_graph->col_edges != NULL);
  assert(db_graph->covgs != NULL);

  FILE *fout = fopen(path, "w");

  if(fout == NULL) {
    die("Unable to open dump binary file to write: %s\n", path);
  }

  GraphInfo *ginfo = &db_graph->ginfo;

  // Construct binary header
  BinaryFileHeader header = {.version = version,
                             .kmer_size = db_graph->kmer_size,
                             .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                             .num_of_colours = num_of_cols,
                             .num_of_shades = num_of_shades,
                             .num_of_kmers = db_graph->ht.unique_kmers};

  uint64_t total_seq_loaded[num_of_cols];
  uint32_t mean_read_lengths[num_of_cols];
  StrBuf *sample_names[num_of_cols];
  long double seq_err_rates[num_of_cols];
  ErrorCleaning *err_cleaning[num_of_cols];

  if(colours != NULL)
  {
    uint32_t i, col;
    for(i = 0; i < num_of_cols; i++) {
      col = colours[i];
      total_seq_loaded[i] = ginfo->total_sequence[col];
      sample_names[i] = ginfo->sample_names[col];
      mean_read_lengths[i] = ginfo->mean_read_length[col];
      seq_err_rates[i] = ginfo->seq_err[col];
      err_cleaning[i] = ginfo->cleaning + col;
    }

    header.total_seq_loaded = total_seq_loaded;
    header.mean_read_lengths = mean_read_lengths;
    header.sample_names = sample_names;
    header.seq_err_rates = seq_err_rates;
  }
  else {
    header.total_seq_loaded = ginfo->total_sequence + start_col;
    header.sample_names = ginfo->sample_names + start_col;
    header.mean_read_lengths = ginfo->mean_read_length + start_col;
    header.seq_err_rates = ginfo->seq_err + start_col;

    uint32_t i;
    for(i = 0; i < num_of_cols; i++)
      err_cleaning[i] = ginfo->cleaning + start_col + i;
  }

  header.err_cleaning = err_cleaning;

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
