#include "global.h"
#include "binary_format.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "range.h"

uint32_t binary_get_num_colours(const char *path, uint32_t ctx_cols)
{
  const char *ptr = strchr(path, ':');
  return (ptr != NULL ? range_get_num(ptr+1, ctx_cols-1) : ctx_cols);
}

void binary_parse_colour_array(const char *path, uint32_t *arr, uint32_t ctx_cols)
{
  // Empty range is the same as :*
  const char *ptr = strchr(path, ':');
  range_parse_array(ptr == NULL ? "" : ptr+1, arr, ctx_cols-1);
}

uint32_t binary_colour_array_alloc(char *ctx_path, uint32_t **result)
{
  boolean is_binary = false;
  uint32_t kmer_size, ctx_num_of_cols;
  uint64_t num_kmers;

  char *sep = strchr(ctx_path, ':');
  if(sep != NULL) *sep = '\0';

  if(!binary_probe(ctx_path, &is_binary, &kmer_size, &ctx_num_of_cols, &num_kmers))
    die("Cannot read input binary file: %s", ctx_path);
  else if(!is_binary)
    die("Input binary file isn't valid: %s", ctx_path);

  uint32_t num_cols, *cols;

  if(sep != NULL) *sep = ':';

  num_cols = binary_get_num_colours(ctx_path, ctx_num_of_cols);
  cols = malloc(num_cols * sizeof(uint32_t));
  binary_parse_colour_array(ctx_path, cols, ctx_num_of_cols);
  *result = cols;
  return num_cols;
}

uint32_t binary_load_colour(const char *path, dBGraph *db_graph,
                            SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                            uint32_t colour)
{
  size_t len = strlen(path);
  char new_path[len+61];
  memcpy(new_path, path, (len+1) * sizeof(char));

  char *end = strchr(new_path, ':');
  if(end == NULL) end = new_path+len;
  sprintf(end, ":%u", colour);

  return binary_load(new_path, db_graph, prefs, stats);
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
                  uint64_t *num_of_kmers_ptr)
{
  *valid_ctx = 0;

  size_t pathlen = strlen(ctx_path);
  char path[pathlen+1];
  memcpy(path, ctx_path, pathlen+1);

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) return 0;

  boolean kmer_count_set = false;
  size_t bytes_per_kmer = 0;
  uint32_t ctx_num_of_cols = 0;

  size_t bytes_read = skip_header(fh, kmer_size_ptr, &ctx_num_of_cols,
                                  num_of_kmers_ptr, &kmer_count_set,
                                  &bytes_per_kmer);

  fclose(fh);

  // No reading errors, but not ctx binary
  if(bytes_read == 0) return 1;

  if(sep != NULL) *sep = ':';
  *num_of_colours_ptr = binary_get_num_colours(path, ctx_num_of_cols);
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
           "[bytes per kmer: %zu remaining: %zu; fsize: %lli; header: %zu]",
           path, bytes_per_kmer, bytes_remaining, fsize, bytes_read);
      *valid_ctx = 0;
    }
  }

  if(sep != NULL) *sep = ':';

  return 1;
}

static void binary_header_alloc(BinaryFileHeader *h)
{
  size_t i;
  h->ginfo = calloc(h->num_of_cols, sizeof(GraphInfo));

  for(i = 0; i < h->num_of_cols; i++)
    graph_info_alloc(h->ginfo + i);
}

void binary_header_destroy(BinaryFileHeader *h)
{
  size_t i;
  for(i = 0; i < h->num_of_cols; i++)
    graph_info_dealloc(h->ginfo + i);

  free(h->ginfo);
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

  binary_header_alloc(h);

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
    die("Kmer has zero covg in all colours in binary: %s", path);

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
uint32_t binary_load(const char *ctx_path, dBGraph *graph,
                     const SeqLoadingPrefs *prefs, SeqLoadingStats *stats)
{
  assert(prefs->must_exist_in_colour == -1 || graph->node_in_cols != NULL);

  size_t len = strlen(ctx_path);
  char path[len+61];
  memcpy(path, ctx_path, (len+1) * sizeof(char));

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) die("Cannot open file: %s\n", path);

  BinaryFileHeader header;
  binary_read_header(fh, &header, path);

  if(sep != NULL) *sep = ':';
  uint32_t num_of_cols = binary_get_num_colours(path, header.num_of_cols);
  uint32_t load_colours[num_of_cols];
  binary_parse_colour_array(path, load_colours, header.num_of_cols);
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
            prefs->into_colour + num_cols_loaded);
  }

  GraphInfo *ginfo = graph->ginfo;

  // Checks for this binary with this executable (kmer size + num colours)
  if(header.kmer_size != graph->kmer_size)
  {
    die("binary has different kmer size [kmer_size: %u vs %u; binary: %s]",
        header.kmer_size, graph->kmer_size, path);
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
      graph_info_merge(ginfo+prefs->into_colour, header.ginfo+load_colours[i]);
  }
  else {
    for(i = 0; i < num_of_cols; i++)
      graph_info_merge(ginfo+prefs->into_colour+i, header.ginfo+load_colours[i]);
  }

  // Update number of colours loaded
  graph->num_of_cols_used = MAX2(graph->num_of_cols_used,
                                 prefs->into_colour + num_cols_loaded);

  // Read kmers
  BinaryKmer bkmer;
  Covg kmercovgs[header.num_of_cols], covgs[num_of_cols];
  Edges kmeredges[header.num_of_cols], edges[num_of_cols];

  size_t num_of_kmers_parsed, num_of_kmers_loaded = 0;
  uint64_t num_of_kmers_already_loaded = graph->ht.unique_kmers;

  for(num_of_kmers_parsed = 0; ; num_of_kmers_parsed++)
  {
    if(binary_read_kmer(fh, &header, path, bkmer, kmercovgs, kmeredges) == 0) break;

    memset(covgs, 0, sizeof(Covg) * num_of_cols);
    memset(edges, 0, sizeof(Edges) * num_of_cols);

    // Collapse down colours
    if(prefs->merge_colours) {
      for(i = 0; i < num_of_cols; i++) {
        covgs[0] += kmercovgs[load_colours[i]];
        edges[0] |= kmeredges[load_colours[i]];
      }
    }
    else {
      uint32_t keep_kmer = 0;
      for(i = 0; i < num_of_cols; i++) {
        covgs[i] = kmercovgs[load_colours[i]];
        edges[i] = kmeredges[load_colours[i]];
        keep_kmer = covgs[i] | edges[i];
      }
      // If kmer has no covg or edges -> don't load
      if(keep_kmer == 0) continue;
    }

    // Fetch node in the de bruijn graph
    hkey_t node;

    if(prefs->must_exist_in_colour >= 0)
    {
      node = hash_table_find(&graph->ht, bkmer);

      if(node != HASH_NOT_FOUND &&
         !db_node_has_col(graph, node, prefs->must_exist_in_colour))
      {
        node = HASH_NOT_FOUND;
      }
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

    if(node != HASH_NOT_FOUND)
    {
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

      if(graph->col_edges != NULL)
      {
        // For each colour take the union of edges
        Edges *col_edges = graph->col_edges + node * graph->num_of_cols;

        if(prefs->must_exist_in_colour >= 0) {
          for(i = 0; i < num_cols_loaded; i++) {
            col = prefs->into_colour+i;
            col_edges[col] |= edges[i] & col_edges[prefs->must_exist_in_colour];
          }
        }
        else {
          for(i = 0; i < num_cols_loaded; i++) {
            col = prefs->into_colour+i;
            col_edges[col] |= edges[i];
          }
        }
      }
      
      if(graph->edges != NULL)
      {
        // Merge all edges into one colour
        for(i = 0; i < num_cols_loaded; i++) {
          graph->edges[node] |= edges[i];
        }
      }

      num_of_kmers_loaded++;
    }
  }

  if(num_of_kmers_parsed > header.num_of_kmers)
  {
    warn("More kmers in binary than expected [expected: %zu; actual: %zu; "
         "path: %s]", (size_t)header.num_of_kmers, (size_t)num_of_kmers_parsed,
        path);
  }

  fclose(fh);

  if(stats != NULL)
  {
    stats->num_of_colours_loaded += num_cols_loaded;
    stats->kmers_loaded += num_of_kmers_loaded;
    stats->unique_kmers += graph->ht.unique_kmers - num_of_kmers_already_loaded;
    for(i = 0; i < num_cols_loaded; i++)
      stats->total_bases_read += header.ginfo[i].total_sequence;
    stats->binaries_loaded++;
  }

  binary_header_destroy(&header);

  return num_cols_loaded;
}
