#include "global.h"
#include "binary_format.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "range.h"


// offset for loading a binary
// e.g. 2:in.ctx means load in.ctx into colours 2..
boolean binary_get_offset(char *path, char ** ptr)
{
  char *str = path;
  for(str = path; *str >= '0' && *str <= '9'; str++);
  if(str > path && *str == ':') { *ptr = str; return true; }
  else return false;
}

uint32_t binary_get_num_colours(const char *path, uint32_t max_col)
{
  const char *ptr = strchr(path, ':');
  return (ptr != NULL ? range_get_num(ptr+1, max_col) : max_col+1);
}

void binary_parse_colour_array(const char *path, uint32_t *arr, uint32_t max_col)
{
  // Empty range is the same as :*
  const char *ptr = strchr(path, ':');
  range_parse_array(ptr == NULL ? "" : ptr+1, arr, max_col);
}

uint32_t binary_load_colour(const char *path,
                            SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                            uint32_t colour)
{
  size_t len = strlen(path);
  char new_path[len+61];
  memcpy(new_path, path, (len+1) * sizeof(char));

  char *end = strchr(new_path, ':');
  if(end == NULL) end = new_path+len;
  sprintf(end, ":%u", colour);

  return binary_load(new_path, prefs, stats, NULL);
}

static int skip_header(FILE *fh, uint32_t *kmer_size_ptr,
                       uint32_t *num_of_colours_ptr,
                       uint64_t *num_of_kmers_ptr, boolean *kmer_count_set,
                       size_t *bytes_per_kmer)
{
  // Estimate number of kmers
  size_t i, skip, read, total = 0;

  *kmer_count_set = false;

  char magic_word[7];
  magic_word[6] = '\0';

  uint32_t version, kmer_size, num_of_bitfields, num_of_colours;
  uint32_t num_of_shades = 0;
  uint64_t num_of_kmers;

  if(fread(magic_word, sizeof(char), 6, fh) != 6 ||
     fread(&version, sizeof(uint32_t), 1, fh) != 1 ||
     fread(&kmer_size, sizeof(uint32_t), 1, fh) != 1 ||
     fread(&num_of_bitfields, sizeof(uint32_t), 1, fh) != 1 ||
     fread(&num_of_colours, sizeof(uint32_t), 1, fh) != 1)
  {
    return 0;
  }

  // Simple check on number of magic word, bitfields vs kmer_size etc.
  if(strcmp(magic_word, CTX_MAGIC_WORD) != 0 ||
     num_of_bitfields * 32 < kmer_size ||
     (num_of_bitfields - 1) * 32 >= kmer_size ||
     num_of_colours == 0)
  {
    return 0;
  }

  *kmer_size_ptr = kmer_size;
  *num_of_colours_ptr = num_of_colours;

  total += 6 + sizeof(uint32_t) * 4;

  if(version >= 7)
  {
    // Read number of kmers
    if(fread(&num_of_kmers, sizeof(uint64_t), 1, fh) != 1 ||
       fread(&num_of_shades, sizeof(uint32_t), 1, fh) != 1)
    {
      return 0;
    }

    *num_of_kmers_ptr = num_of_kmers;
    *kmer_count_set = true;
    total += sizeof(uint64_t) + sizeof(uint32_t);
  }

  // Skip mean read length and total seq per colour
  skip = (sizeof(uint32_t)+sizeof(uint64_t))*num_of_colours;
  read = stream_skip(fh, skip);
  total += read;
  if(read != skip) return 0;

  if(version >= 6)
  {
    // Sample names
    for(i = 0; i < num_of_colours; i++)
    {
      if(fread(&skip, sizeof(uint32_t), 1, fh) != 1) return 0;
      total += sizeof(uint32_t);
      if((read = stream_skip(fh, skip)) != skip) return 0;
      total += read;
    }

    // Sequencing error rates
    skip = num_of_colours * sizeof(long double);
    if((read = stream_skip(fh, skip)) != skip) return 0;
    total += read;

    // Sample cleaning
    for(i = 0; i < num_of_colours; i++)
    {
      skip = sizeof(uint8_t)*4+sizeof(uint32_t)*2;
      if((read = stream_skip(fh, skip)) != skip) return 0;
      total += read;

      if(fread(&skip, sizeof(uint32_t), 1, fh) != 1) return 0;
      total += sizeof(uint32_t);
      if((read = stream_skip(fh, skip)) != skip) return 0;
      total += read;
    }
  }

  // 'CORTEX' end of header
  if(fread(magic_word, sizeof(char), 6, fh) != 6 ||
     strcmp(magic_word, CTX_MAGIC_WORD) != 0)
  {
    return 0;
  }
  total += 6;

  size_t shade_bytes = num_of_shades / 8;

  *bytes_per_kmer = sizeof(uint64_t) * num_of_bitfields +
                    sizeof(uint32_t) * num_of_colours + // coverage
                    sizeof(uint8_t) * num_of_colours + // edges
                    sizeof(uint8_t) * shade_bytes * 2; // shades

  return total;
}

// returns 0 if cannot read, 1 otherwise
char binary_probe(const char *ctx_path, boolean *valid_ctx,
                  uint32_t *kmer_size_ptr, uint32_t *num_of_colours_ptr,
                  uint32_t *max_col_index, uint64_t *num_of_kmers_ptr)
{
  *valid_ctx = 0;

  size_t pathlen = strlen(ctx_path);
  char path[pathlen+1];
  memcpy(path, ctx_path, pathlen+1);

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) return 0;

  uint32_t ctx_num_of_cols = 0;
  boolean kmer_count_set = false;
  size_t bytes_per_kmer = 0;

  size_t bytes_read = skip_header(fh, kmer_size_ptr, &ctx_num_of_cols,
                                  num_of_kmers_ptr, &kmer_count_set,
                                  &bytes_per_kmer);

  fclose(fh);

  // No reading errors, but not ctx binary
  if(bytes_read == 0) return 1;

  *max_col_index = ctx_num_of_cols-1;

  if(sep != NULL) *sep = ':';
  *num_of_colours_ptr = binary_get_num_colours(path, *max_col_index);
  if(sep != NULL) *sep = '\0';

  // Valid ctx binary
  *valid_ctx = 1;

  if(!kmer_count_set)
  {
    // Count number of kmers based on file size
    off_t fsize = get_file_size(path);
    size_t bytes_remaining = fsize - bytes_read;
    *num_of_kmers_ptr = (bytes_remaining / bytes_per_kmer);

    if(bytes_remaining % bytes_per_kmer != 0) {
      warn("Truncated ctx binary: %s "
           "[bytes per kmer: %zu remaining: %zu; fsize: %zu; header: %zu]",
           path, bytes_per_kmer, bytes_remaining, (size_t)fsize, bytes_read);
      *valid_ctx = 0;
    }
  }

  if(sep != NULL) *sep = ':';

  return 1;
}

void binary_header_alloc(BinaryFileHeader *h, size_t num_of_cols)
{
  size_t i;
  h->capacity = num_of_cols;
  if(num_of_cols == 0) return;
  h->ginfo = calloc2(h->capacity, sizeof(GraphInfo));
  for(i = 0; i < h->capacity; i++)
    graph_info_alloc(h->ginfo + i);
}

void binary_header_realloc(BinaryFileHeader *h, size_t num_of_cols)
{
  size_t i;
  if(num_of_cols < h->capacity) return;
  h->ginfo = realloc2(h->ginfo, num_of_cols * sizeof(GraphInfo));
  for(i = h->capacity; i < num_of_cols; i++)
    graph_info_alloc(h->ginfo + i);
  h->capacity = num_of_cols;
}

void binary_header_dealloc(BinaryFileHeader *h)
{
  size_t i;
  for(i = 0; i < h->capacity; i++)
    graph_info_dealloc(h->ginfo + i);
  free(h->ginfo);
}

void binary_read_cpy_basic(BinaryFileHeader *dst, BinaryFileHeader *src)
{
  dst->version = src->version;
  dst->kmer_size = src->kmer_size;
  dst->num_of_bitfields = src->num_of_bitfields;
  dst->num_of_cols = src->num_of_cols;
  dst->num_of_kmers = src->num_of_kmers;
}

// Return number of bytes read
size_t binary_read_header(FILE *fh, BinaryFileHeader *h, const char *path)
{
  size_t i, bytes_read = 0;

  char magic_word[7];
  magic_word[6] = '\0';

  safe_fread(fh, magic_word, strlen(CTX_MAGIC_WORD), "Magic word", path);
  if(strcmp(magic_word, CTX_MAGIC_WORD) != 0) {
    die("Magic word doesn't match '%s' (start): %s", CTX_MAGIC_WORD, path);
  }
  bytes_read += strlen(CTX_MAGIC_WORD);

  // Read version number, kmer_size, num bitfields, num colours
  safe_fread(fh, &h->version, sizeof(uint32_t), "binary version", path);
  safe_fread(fh, &h->kmer_size, sizeof(uint32_t), "kmer size", path);
  safe_fread(fh, &h->num_of_bitfields, sizeof(uint32_t), "number of bitfields", path);
  safe_fread(fh, &h->num_of_cols, sizeof(uint32_t), "number of colours", path);
  bytes_read += 4*sizeof(uint32_t);

  // Checks
  if(h->version > 7 || h->version < 4)
  {
    die("Sorry, we only support binary versions 4, 5, 6 & 7 "
        "[version: %i; binary: %s]\n", h->version, path);
  }

  if(h->kmer_size % 2 == 0)
  {
    die("kmer size is not an odd number [kmer_size: %i; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3)
  {
    die("kmer size is less than three [kmer_size: %i; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_bitfields * 32 < h->kmer_size)
  {
    die("Not enough bitfields for kmer size "
        "[kmer_size: %i; bitfields: %i; binary: %s]\n",
        h->kmer_size, h->num_of_bitfields, path);
  }

  if((h->num_of_bitfields-1)*32 >= h->kmer_size)
    die("using more than the minimum number of bitfields [binary: %s]\n", path);

  if(h->num_of_cols == 0)
    die("number of colours is zero [binary: %s]\n", path);

  if(h->capacity == 0)
    binary_header_alloc(h, h->num_of_cols);
  else if(h->num_of_cols > h->capacity)
    binary_header_realloc(h, h->num_of_cols);

  if(h->version >= 7)
  {
    safe_fread(fh, &h->num_of_kmers, sizeof(uint64_t), "number of kmers", path);
    uint32_t tmp;
    safe_fread(fh, &tmp, sizeof(uint32_t), "number of shades", path);
    bytes_read += sizeof(uint64_t) + sizeof(uint32_t);

    if((tmp & 0x7) != 0) {
      warn("Number of shades is not a multiple of 8 [binary: %s]", path);
    }
  }

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
             "(premature '\\0') [binary: %s]\n", i, len, len2, path);
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

      safe_fread(fh, &(err_cleaning->tip_clipping),
               sizeof(uint8_t), "tip cleaning", path);
      safe_fread(fh, &(err_cleaning->remv_low_cov_sups),
               sizeof(uint8_t), "remove low covg supernodes", path);
      safe_fread(fh, &(err_cleaning->remv_low_cov_nodes),
               sizeof(uint8_t), "remove low covg kmers", path);
      safe_fread(fh, &(err_cleaning->cleaned_against_another_graph),
               sizeof(uint8_t), "cleaned against graph", path);

      safe_fread(fh, &(err_cleaning->remv_low_cov_sups_thresh),
               sizeof(uint32_t), "remove low covg supernode threshold", path);
      safe_fread(fh, &(err_cleaning->remv_low_cov_nodes_thresh),
               sizeof(uint32_t), "remove low covg kmer threshold", path);

      bytes_read += 4*sizeof(uint8_t) + 2*sizeof(uint32_t);

      // Fix for old versions with negative thresholds
      if(h->version <= 6) {
        if(!err_cleaning->remv_low_cov_sups &&
           err_cleaning->remv_low_cov_sups_thresh == (uint32_t)-1) {
          err_cleaning->remv_low_cov_nodes_thresh = 0;
        } else if(!err_cleaning->remv_low_cov_nodes &&
           err_cleaning->remv_low_cov_nodes_thresh == (uint32_t)-1) {
          err_cleaning->remv_low_cov_nodes_thresh = 0;
        }
      }

      // Sanity checks
      if(!err_cleaning->remv_low_cov_sups &&
         err_cleaning->remv_low_cov_sups_thresh > 0)
      {
        warn("Binary header gives cleaning threshold for supernodes "
             "when no cleaning was performed [binary: %s]", path);

        err_cleaning->remv_low_cov_sups_thresh = 0;
      }

      if(!err_cleaning->remv_low_cov_nodes &&
         err_cleaning->remv_low_cov_nodes_thresh > 0)
      {
        warn("Binary header gives cleaning threshold for nodes "
             "when no cleaning was performed [binary: %s]", path);

        err_cleaning->remv_low_cov_nodes_thresh = 0;
      }

      // Read cleaned against name
      uint32_t len;
      safe_fread(fh, &len, sizeof(uint32_t), "graph name length", path);

      StrBuf *sbuf = &err_cleaning->cleaned_against_graph_name;
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
             "(premature '\\0') [binary: %s]\n", i, len, len2, path);
      }
    }
  }

  // Read magic word at the end of header 'CORTEX'
  safe_fread(fh, magic_word, strlen(CTX_MAGIC_WORD), "magic word (end)", path);
  if(strcmp(magic_word, CTX_MAGIC_WORD) != 0) {
    die("Magic word doesn't match '%s' (end): '%s' [binary: %s]\n",
        CTX_MAGIC_WORD, magic_word, path);
  }
  bytes_read += strlen(CTX_MAGIC_WORD);

  // If we haven't set num_of_kmers set it now using file size
  if(h->version < 7)
  {
    off_t file_size = get_file_size(path);
    size_t bytes_remaining = file_size - bytes_read;
    size_t shade_bytes = 0;

    // 2 * num_shade_bytes for shade + shade end data
    size_t num_bytes_per_kmer
      = sizeof(uint64_t) * NUM_BITFIELDS_IN_BKMER +
        sizeof(uint32_t) * h->num_of_cols + // coverage
        sizeof(uint8_t) * h->num_of_cols + // edges
        sizeof(uint8_t) * shade_bytes * 2; // shades

    h->num_of_kmers = bytes_remaining / num_bytes_per_kmer;

    if(num_bytes_per_kmer * h->num_of_kmers != bytes_remaining) {
      die("Irregular size of binary file (corrupted?): %s", path);
    }
  }

  return bytes_read;
}

size_t binary_read_kmer(FILE *fh, BinaryFileHeader *h, const char *path,
                        uint64_t *bkmer, Covg *covgs, Edges *edges)
{
  size_t i, num_bytes_read;

  num_bytes_read = fread(bkmer, 1, sizeof(uint64_t)*h->num_of_bitfields, fh);

  if(num_bytes_read == 0) return 0;
  if(num_bytes_read != sizeof(uint64_t)*h->num_of_bitfields)
    die("Unexpected end of file: %s", path);

  safe_fread(fh, covgs, h->num_of_cols * sizeof(uint32_t), "Coverages", path);
  safe_fread(fh, edges, h->num_of_cols * sizeof(uint8_t), "Edges", path);
  num_bytes_read += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint8_t));

  // Check top word of each kmer
  uint64_t top_word_mask = ~(uint64_t)0 << BKMER_TOP_BITS(h->kmer_size);
  if(bkmer[0] & top_word_mask) die("Oversized kmer in binary: %s", path);

  // Check covg is not 0 for all colours
  for(i = 0; i < h->num_of_cols && covgs[i] == 0; i++) {}
  if(i == h->num_of_cols)
    warn("Kmer has zero covg in all colours in binary: %s", path);

  return num_bytes_read;
}

// if only_load_if_in_colour is >= 0, only kmers with coverage in existing
// colour only_load_if_in_colour will be loaded.
// We assume only_load_if_in_colour < load_first_colour_into
// if all_kmers_are_unique != 0 an error is thrown if a node already exists
// returns the number of colours in the binary
// If stats != NULL, updates:
//   stats->num_of_colours_loaded
//   stats->kmers_loaded
//   stats->total_bases_read
//   stats->binaries_loaded
uint32_t binary_load(const char *ctx_path,
                     const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                     BinaryFileHeader *header_ptr)
{
  dBGraph *graph = prefs->db_graph;

  size_t len = strlen(ctx_path);
  char path[len+61];
  memcpy(path, ctx_path, (len+1) * sizeof(char));

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) die("Cannot open file: %s\n", path);

  BinaryFileHeader header_mem = {.capacity = 0};

  BinaryFileHeader *header = header_ptr != NULL ? header_ptr : &header_mem;
  binary_read_header(fh, header, path);

  if(sep != NULL) *sep = ':';
  uint32_t num_of_cols = binary_get_num_colours(path, header->num_of_cols-1);
  uint32_t load_colours[num_of_cols];
  binary_parse_colour_array(path, load_colours, header->num_of_cols-1);
  if(sep != NULL) *sep = '\0';

  uint32_t num_cols_loaded = prefs->merge_colours ? 1 : num_of_cols;

  uint32_t i;
  Colour col;

  message("Loading binary %s", path);
  if(sep != NULL) {
    message(" with colour filter: %u", load_colours[0]);
    for(i = 1; i < num_of_cols; i++) message(",%u", load_colours[i]);
  }
  if(prefs->merge_colours) message(" into colour %u\n", prefs->into_colour);
  else {
    message(" into colours %u-%u\n", prefs->into_colour,
            prefs->into_colour+num_cols_loaded-1);
  }

  GraphInfo *ginfo = graph->ginfo;

  // Checks for this binary with this executable (kmer size + num colours)
  if(header->kmer_size != graph->kmer_size)
  {
    die("binary has different kmer size [kmer_size: %u vs %u; binary: %s]",
        header->kmer_size, graph->kmer_size, path);
  }

  if(prefs->into_colour + num_cols_loaded > graph->num_of_cols &&
     (graph->col_covgs != NULL || graph->col_edges != NULL))
  {
    die("Program has not assigned enough colours! "
        "[colours in binary: %i vs graph: %i, load into: %i; binary: %s]",
        num_cols_loaded, graph->num_of_cols, prefs->into_colour, path);
  }

  if(prefs->merge_colours) {
    for(i = 0; i < num_of_cols; i++)
      graph_info_merge(ginfo+prefs->into_colour, header->ginfo+load_colours[i]);
  }
  else {
    for(i = 0; i < num_of_cols; i++)
      graph_info_merge(ginfo+prefs->into_colour+i, header->ginfo+load_colours[i]);
  }

  // Update number of colours loaded
  graph->num_of_cols_used = MAX2(graph->num_of_cols_used,
                                 prefs->into_colour + num_cols_loaded);

  // Read kmers
  BinaryKmer bkmer;
  Covg kmercovgs[header->num_of_cols], covgs[num_of_cols];
  Edges kmeredges[header->num_of_cols], edges[num_of_cols];

  size_t num_of_kmers_parsed, num_of_kmers_loaded = 0;
  uint64_t num_of_kmers_already_loaded = graph->ht.unique_kmers;

  for(num_of_kmers_parsed = 0; ; num_of_kmers_parsed++)
  {
    if(binary_read_kmer(fh, header, path, bkmer, kmercovgs, kmeredges) == 0) break;

    // Collapse down colours
    Covg keep_kmer = 0;
    if(prefs->merge_colours) {
      covgs[0] = 0;
      edges[0] = 0;
      for(i = 0; i < num_of_cols; i++) {
        covgs[0] += kmercovgs[load_colours[i]];
        edges[0] |= kmeredges[load_colours[i]];
        keep_kmer |= covgs[0] | edges[0];
      }
    }
    else {
      for(i = 0; i < num_of_cols; i++) {
        covgs[i] = kmercovgs[load_colours[i]];
        edges[i] = kmeredges[load_colours[i]];
        keep_kmer |= covgs[i] | edges[i];
      }
    }
    // If kmer has no covg or edges -> don't load
    if(keep_kmer == 0) continue;

    if(prefs->boolean_covgs)
      for(i = 0; i < num_of_cols; i++)
        covgs[i] = covgs[i] > 0;

    // Fetch node in the de bruijn graph
    hkey_t node;

    if(prefs->must_exist_in_graph)
    {
      node = hash_table_find(&graph->ht, bkmer);
      if(node == HASH_NOT_FOUND) continue;

      Edges union_edges
        = graph->edges != NULL ? db_node_edges(graph, node)
                                  : db_node_col_edges_union(graph, node);

      for(i = 0; i < num_of_cols; i++) edges[i] &= union_edges;
    }
    else
    {
      boolean found;
      node = hash_table_find_or_insert(&graph->ht, bkmer, &found);

      if(prefs->empty_colours && found)
      {
        die("Duplicate kmer loaded [cols:%u:%u]",
            prefs->into_colour, num_of_cols);
      }
    }

    // Set presence in colours
    if(graph->node_in_cols != NULL) {
      for(i = 0; i < num_cols_loaded; i++)
        if(covgs[i] > 0 || edges[i] != 0)
          db_node_set_col(graph, node, prefs->into_colour+i);
    }

    if(graph->col_covgs != NULL) {
      for(i = 0; i < num_cols_loaded; i++)
        db_node_add_col_covg(graph, node, prefs->into_colour+i, covgs[i]);
    }

    if(graph->col_edges != NULL) {
      Edges *col_edges = graph->col_edges + node * graph->num_of_cols;
      for(i = 0; i < num_cols_loaded; i++) {
        col = prefs->into_colour+i;
        col_edges[col] |= edges[i];
      }
    }

    // Merge all edges into one colour
    if(graph->edges != NULL) {
      for(i = 0; i < num_cols_loaded; i++)
        graph->edges[node] |= edges[i];
    }

    num_of_kmers_loaded++;
  }

  if(num_of_kmers_parsed > header->num_of_kmers)
  {
    warn("More kmers in binary than expected [expected: %zu; actual: %zu; "
         "path: %s]", (size_t)header->num_of_kmers, (size_t)num_of_kmers_parsed,
        path);
  }

  fclose(fh);

  if(stats != NULL)
  {
    stats->num_of_colours_loaded += num_cols_loaded;
    stats->kmers_loaded += num_of_kmers_loaded;
    stats->unique_kmers += graph->ht.unique_kmers - num_of_kmers_already_loaded;
    for(i = 0; i < num_cols_loaded; i++)
      stats->total_bases_read += header->ginfo[i].total_sequence;
    stats->binaries_loaded++;
  }

  if(header_ptr == NULL)
    binary_header_dealloc(header);

  return num_cols_loaded;
}

//
// Merging, filtering, combining binary graph files
//

// Filter a binary against the graph
// (only keep nodes and edges that are in the graph)
// Same functionality as binaries_merge, but faster if dealing with only one
// input file
size_t binary_filter_graph(const char *out_ctx_path, char *in_ctx_path,
                           boolean flatten, const dBGraph *db_graph)
{
  assert(db_graph->edges != NULL || db_graph->col_edges != NULL);

  size_t i, nodes_dumped = 0;
  FILE *in, *out;
  char *tmp_path, *endptr;
  uint32_t offset = 0;

  if(binary_get_offset(in_ctx_path, &tmp_path))
  {
    offset = strtoul(in_ctx_path, &endptr, 10);
    in_ctx_path = tmp_path+1;
  }

  char *split = strchr(in_ctx_path, ':');
  if(split != NULL) *split = '\0';

  if((in = fopen(in_ctx_path, "r")) == NULL)
    die("Cannot open input path: %s", in_ctx_path);
  if((out = fopen(out_ctx_path, "w")) == NULL)
    die("Cannot open output path: %s", out_ctx_path);

  BinaryFileHeader inheader = {.capacity = 0};
  binary_read_header(in, &inheader, in_ctx_path);

  if(split != NULL) *split = ':';
  uint32_t num_of_cols = binary_get_num_colours(in_ctx_path, inheader.num_of_cols-1);
  uint32_t load_colours[num_of_cols];
  binary_parse_colour_array(in_ctx_path, load_colours, inheader.num_of_cols-1);

  message("Loading binary %s", in_ctx_path);
  if(split != NULL) {
    message(" with colour filter: %u", load_colours[0]);
    for(i = 1; i < num_of_cols; i++) message(",%u", load_colours[i]);
  }
  if(flatten || num_of_cols == 1) message(" into colour %u\n", offset);
  else message(" into colours %u-%u\n", offset, offset+num_of_cols-1);

  BinaryFileHeader outheader = {.capacity = 0};
  binary_read_cpy_basic(&outheader, &inheader);
  outheader.num_of_cols = offset+num_of_cols;
  binary_header_alloc(&outheader, outheader.num_of_cols);

  for(i = 0; i < num_of_cols; i++)
    graph_info_merge(outheader.ginfo+offset+i, inheader.ginfo + load_colours[i]);

  binary_write_header(out, &outheader);

  BinaryKmer bkmer;
  Covg kmercovgs[inheader.num_of_cols], covgs[offset+num_of_cols];
  Edges kmeredges[inheader.num_of_cols], edges[offset+num_of_cols];
  memset(covgs, 0, (offset+num_of_cols)*sizeof(Covg));
  memset(edges, 0, (offset+num_of_cols)*sizeof(Edges));

  while(binary_read_kmer(in, &inheader, in_ctx_path, bkmer, kmercovgs, kmeredges))
  {
    // Collapse down colours
    Covg keep_kmer = 0;
    if(flatten) {
      covgs[0] = 0;
      edges[0] = 0;
      for(i = 0; i < num_of_cols; i++) {
        covgs[offset] += kmercovgs[load_colours[i]];
        edges[offset] |= kmeredges[load_colours[i]];
        keep_kmer |= covgs[offset] | edges[offset];
      }
    }
    else {
      for(i = 0; i < num_of_cols; i++) {
        covgs[offset+i] = kmercovgs[load_colours[i]];
        edges[offset+i] = kmeredges[load_colours[i]];
        keep_kmer |= covgs[offset+i] | edges[offset+i];
      }
    }
    // If kmer has no covg or edges -> don't load
    if(keep_kmer == 0) continue;

    hkey_t node = hash_table_find(&db_graph->ht, bkmer);
    Edges union_edges;

    if(node != HASH_NOT_FOUND)
    {
      union_edges = (db_graph->edges != NULL) ? db_node_edges(db_graph, node)
                                              : db_node_col_edges_union(db_graph,
                                                                        node);

      for(i = 0; i < num_of_cols; i++) edges[offset+i] &= union_edges;
      binary_write_kmer(out, &outheader, bkmer, covgs, edges);
      nodes_dumped++;
    }
  }

  fclose(in);
  fclose(out);

  message("Dumped %zu kmers in %u colour%s into: %s\n",
          nodes_dumped, outheader.num_of_cols,
          outheader.num_of_cols != 1 ? "s" : "", out_ctx_path);

  binary_header_dealloc(&inheader);
  binary_header_dealloc(&outheader);

  return nodes_dumped;
}



// ctx_num_cols and ctx_max_cols are the numbers returned from binary_probe
// if merge pool colour 0 from each binary into colour 0, 1 -> 1 etc.
// if flatten, pool all colours into colour 0
// if intersect only load kmers that are already in the hash table
size_t binaries_merge(char *out_ctx_path, char **binary_paths,
                      size_t num_binaries,
                      uint32_t ctx_num_cols[num_binaries],
                      uint32_t ctx_max_cols[num_binaries],
                      boolean merge, boolean flatten, boolean intersect,
                      dBGraph *db_graph)
{
  assert(!merge || !flatten);
  assert(intersect || db_graph->ht.unique_kmers == 0);

  if(num_binaries == 1 && intersect)
    return binary_filter_graph(out_ctx_path, binary_paths[0], flatten, db_graph);

  uint32_t i, offsets[num_binaries];
  uint32_t max_cols = 0, sum_cols = 0;

  // Check all binaries are valid binaries with matching kmer size
  char *ptr;
  char *endptr;

  for(i = 0; i < num_binaries; i++)
  {
    if(binary_get_offset(binary_paths[i], &ptr)) {
      offsets[i] = strtoul(binary_paths[i], &endptr, 10);
      binary_paths[i] = ptr+1;
    } else {
      offsets[i] = 0;
    }

    max_cols = MAX2(offsets[i] + ctx_num_cols[i], max_cols);
    sum_cols += ctx_num_cols[i];
  }

  uint32_t output_colours;

  if(flatten) output_colours = 1;
  else if(merge) output_colours = max_cols;
  else output_colours = sum_cols;

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs
    = {.into_colour = 0, .merge_colours = true,
       .boolean_covgs = false,
       .load_seq = true,
       .quality_cutoff = 0, .ascii_fq_offset = 0,
       .homopolymer_cutoff = 0,
       .remove_dups_se = false, .remove_dups_pe = false,
       .load_binaries = true,
       .must_exist_in_graph = intersect, .empty_colours = false,
       .update_ginfo = true,
       .db_graph = db_graph};
  //

  if(output_colours == 1)
  {
    // e.g. flatten
    for(i = 0; i < num_binaries; i++)
      binary_load(binary_paths[i], &prefs, stats, NULL);

    hash_table_print_stats(&db_graph->ht);
    binary_dump_graph(out_ctx_path, db_graph, CURR_CTX_VERSION, NULL, 0, 1);
  }
  else
  {
    uint32_t load_colours[num_binaries][max_cols];
    for(i = 0; i < num_binaries; i++)
      binary_parse_colour_array(binary_paths[i], load_colours[i], ctx_max_cols[i]);

    // Construct binary header
    BinaryFileHeader tmpheader;
    BinaryFileHeader output_header = {.version = CURR_CTX_VERSION,
                                      .kmer_size = db_graph->kmer_size,
                                      .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                                      .num_of_cols = output_colours,
                                      .num_of_kmers = db_graph->ht.unique_kmers};

    binary_header_alloc(&tmpheader, max_cols);
    binary_header_alloc(&output_header, output_colours);

    Colour j, output_colour = 0;
    for(i = 0; i < num_binaries; i++)
    {
      binary_load(binary_paths[i], &prefs, stats, &tmpheader);
      if(merge) output_colour = 0;
      output_colour += offsets[i];
      for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
      {
        graph_info_merge(output_header.ginfo + output_colour,
                         tmpheader.ginfo + load_colours[i][j]);
      }
    }

    FILE *fh = fopen(out_ctx_path, "w");
    if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

    size_t header_size = binary_write_header(fh, &output_header);

    // Free header resources
    binary_header_dealloc(&tmpheader);
    binary_header_dealloc(&output_header);

    // print file outline
    message("Generated merged hash table\n\n");
    hash_table_print_stats(&db_graph->ht);
    dump_empty_binary(db_graph, fh, output_colours);

    if(merge)
    {
      for(output_colour = 0; output_colour < output_colours; output_colour++)
      {
        memset(db_graph->col_edges, 0, db_graph->ht.capacity * sizeof(Edges));
        memset(db_graph->col_covgs, 0, db_graph->ht.capacity * sizeof(Covg));

        boolean data_loaded_in_col = false;
        for(i = 0; i < num_binaries; i++)
        {
          if(output_colour >= offsets[i] &&
             output_colour < offsets[i] + ctx_num_cols[i])
          {
            uint32_t ctx_col = output_colour - offsets[i];
            binary_load_colour(binary_paths[i], &prefs, stats,
                               load_colours[i][ctx_col]);
            data_loaded_in_col = true;
          }
        }
        if(data_loaded_in_col) {
          message("Dumping into colour %u...\n\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          binary_dump_colour(db_graph, 0, output_colour, output_colours, fh);
        }
      }
    }
    else
    {
      output_colour = 0;
      for(i = 0; i < num_binaries; i++)
      {
        output_colour += offsets[i];
        for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
        {
          memset(db_graph->col_edges, 0, db_graph->ht.capacity * sizeof(Edges));
          memset(db_graph->col_covgs, 0, db_graph->ht.capacity * sizeof(Covg));

          binary_load_colour(binary_paths[i], &prefs, stats,
                             load_colours[i][j]);

          message("Dumping into colour %u...\n\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          binary_dump_colour(db_graph, 0, output_colour, output_colours, fh);
        }
      }
    }

    fclose(fh);
  }

  seq_loading_stats_free(stats);

  message("Dumped %zu kmers in %u colour%s into: %s\n",
          (size_t)db_graph->ht.unique_kmers, output_colours,
          output_colours != 1 ? "s" : "", out_ctx_path);

  return db_graph->ht.unique_kmers;
}

