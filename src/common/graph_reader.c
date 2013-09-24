#include "global.h"
#include "graph_format.h"
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

uint32_t graph_file_get_ncols(const char *path, uint32_t max_col)
{
  const char *ptr = strchr(path, ':');
  return (ptr != NULL ? range_get_num(ptr+1, max_col) : max_col+1);
}

void graph_file_parse_colours(const char *path, uint32_t *arr, uint32_t max_col)
{
  // Empty range is the same as :*
  const char *ptr = strchr(path, ':');
  range_parse_array(ptr == NULL ? "" : ptr+1, arr, max_col);
}

uint32_t graph_load_colour(const char *path, const SeqLoadingPrefs *prefs,
                           SeqLoadingStats *stats, uint32_t colour)
{
  size_t len = strlen(path);
  char new_path[len+61];
  memcpy(new_path, path, (len+1) * sizeof(char));

  char *end = strchr(new_path, ':');
  if(end == NULL) end = new_path+len;
  sprintf(end, ":%u", colour);

  return graph_load(new_path, prefs, stats, NULL);
}

void graph_header_alloc(GraphFileHeader *h, size_t num_of_cols)
{
  size_t i, old_cap = h->capacity;

  if(h->capacity == 0)
    h->ginfo = calloc2(num_of_cols, sizeof(GraphInfo));
  else
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

void graph_header_cpy(GraphFileHeader *dst, const GraphFileHeader *src)
{
  size_t i;
  dst->version = src->version;
  dst->kmer_size = src->kmer_size;
  dst->num_of_bitfields = src->num_of_bitfields;
  dst->num_of_cols = src->num_of_cols;
  dst->max_col = src->max_col;
  dst->num_of_kmers = src->num_of_kmers;
  graph_header_alloc(dst, dst->num_of_cols);
  for(i = 0; i < dst->num_of_cols; i++)
    graph_info_cpy(&dst->ginfo[i], &src->ginfo[i]);
}

// Return number of bytes read or -1 if not valid
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
  SAFE_READ(fh, &h->version, sizeof(uint32_t), "binary version", path, fatal);
  SAFE_READ(fh, &h->kmer_size, sizeof(uint32_t), "kmer size", path, fatal);
  SAFE_READ(fh, &h->num_of_bitfields, sizeof(uint32_t), "num of bitfields", path, fatal);
  SAFE_READ(fh, &h->num_of_cols, sizeof(uint32_t), "number of colours", path, fatal);
  bytes_read += 4*sizeof(uint32_t);

  // Checks
  if(h->version > 7 || h->version < 4)
  {
    if(!fatal) return -1;
    die("Sorry, we only support binary versions 4, 5, 6 & 7 "
        "[version: %u; binary: %s]\n", h->version, path);
  }

  if(h->kmer_size % 2 == 0)
  {
    if(!fatal) return -1;
    die("kmer size is not an odd number [kmer_size: %u; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3)
  {
    if(!fatal) return -1;
    die("kmer size is less than three [kmer_size: %u; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_bitfields * 32 < h->kmer_size)
  {
    if(!fatal) return -1;
    die("Not enough bitfields for kmer size "
        "[kmer_size: %u; bitfields: %u; binary: %s]\n",
        h->kmer_size, h->num_of_bitfields, path);
  }

  if((h->num_of_bitfields-1)*32 >= h->kmer_size) {
    if(!fatal) return -1;
    die("using more than the minimum number of bitfields [binary: %s]\n", path);
  }

  if(h->num_of_cols == 0) {
    if(!fatal) return -1;
    die("number of colours is zero [binary: %s]\n", path);
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
      warn("Number of shades is not a multiple of 8 [binary: %s]", path);
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
      safe_fread(fh, &(err_cleaning->is_graph_intersection),
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
             "(premature '\\0') [binary: %s]\n", i, len, len2, path);
      }
    }
  }

  // Read magic word at the end of header 'CORTEX'
  safe_fread(fh, magic_word, strlen("CORTEX"), "magic word (end)", path);
  if(strcmp(magic_word, "CORTEX") != 0)
  {
    if(!fatal) return -1;
    die("Magic word doesn't match '%s' (end): '%s' [binary: %s]\n",
        "CORTEX", magic_word, path);
  }
  bytes_read += strlen("CORTEX");

  // If we haven't set num_of_kmers set it now using file size
  if(h->version < 7)
  {
    off_t file_size = get_file_size(path);
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
      die("Irregular size of binary file (corrupted?): %s", path);
    }
  }

  return bytes_read;
}

size_t graph_file_read_kmer(FILE *fh, GraphFileHeader *h, const char *path,
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

// returns 0 if cannot read, 1 otherwise
boolean graph_file_probe(const char* ctx_path, boolean *valid_ctx,
                         GraphFileHeader *gheader)
{
  *valid_ctx = false;

  size_t pathlen = strlen(ctx_path);
  char path[pathlen+1];
  memcpy(path, ctx_path, pathlen+1);

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) return false;
  setvbuf(fh, NULL, _IOFBF, CTX_BUF_SIZE);
  int hret = graph_file_read_header(fh, gheader, false, path);
  fclose(fh);

  // Get num of colours from path
  if(hret > 0)
  {
    *valid_ctx = true;
    if(sep != NULL) {
      *sep = ':';
      gheader->num_of_cols = graph_file_get_ncols(path, gheader->num_of_cols-1);
    }
    gheader->max_col = gheader->num_of_cols-1;
  }

  return true;
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
uint32_t graph_load(const char *ctx_path,
                     const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                     GraphFileHeader *header_ptr)
{
  dBGraph *graph = prefs->db_graph;

  size_t len = strlen(ctx_path);
  char path[len+61];
  memcpy(path, ctx_path, (len+1) * sizeof(char));

  char *sep = strchr(path, ':');
  if(sep != NULL) *sep = '\0';

  FILE* fh = fopen(path, "r");
  if(fh == NULL) die("Cannot open file: %s\n", path);
  setvbuf(fh, NULL, _IOFBF, CTX_BUF_SIZE);

  GraphFileHeader header_mem = {.capacity = 0};

  GraphFileHeader *header = header_ptr != NULL ? header_ptr : &header_mem;
  graph_file_read_header(fh, header, true, path);

  if(sep != NULL) *sep = ':';
  size_t num_of_cols = graph_file_get_ncols(path, header->num_of_cols-1);
  uint32_t load_colours[num_of_cols];
  graph_file_parse_colours(path, load_colours, header->num_of_cols-1);
  if(sep != NULL) *sep = '\0';

  size_t num_cols_loaded = prefs->merge_colours ? 1 : num_of_cols;

  size_t i;
  Colour col;

  timestamp(ctx_msg_out);
  message(" Loading binary %s", path);
  if(sep != NULL) {
    message(" with colour filter: %u", load_colours[0]);
    for(i = 1; i < num_of_cols; i++) message(",%u", load_colours[i]);
  }
  if(prefs->merge_colours || num_cols_loaded == 1)
    message(" into colour %zu\n", prefs->into_colour);
  else {
    message(" into colours %zu-%zu\n", prefs->into_colour,
            prefs->into_colour+num_cols_loaded-1);
  }

  GraphInfo *ginfo = graph->ginfo;

  // Checks for this binary with this executable (kmer size + num colours)
  if(header->kmer_size != graph->kmer_size)
  {
    die("binary has different kmer size [kmer_size: %u vs %zu; binary: %s]",
        header->kmer_size, graph->kmer_size, path);
  }

  if(prefs->into_colour + num_cols_loaded > graph->num_of_cols &&
     (graph->col_covgs != NULL || graph->col_edges != NULL))
  {
    die("Program has not assigned enough colours! "
        "[colours in binary: %zu vs graph: %zu, load into: %zu; binary: %s]",
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
    if(graph_file_read_kmer(fh, header, path, bkmer.b, kmercovgs, kmeredges) == 0) break;

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

      Edges union_edges = db_node_col_edges_union(graph, node);

      for(i = 0; i < num_of_cols; i++) edges[i] &= union_edges;
    }
    else
    {
      boolean found;
      node = hash_table_find_or_insert(&graph->ht, bkmer, &found);

      if(prefs->empty_colours && found)
      {
        die("Duplicate kmer loaded [cols:%zu:%zu]",
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

    // This may be an invalid pointer (if num_edge_cols == 0)
    Edges *col_edges = graph->col_edges + node * graph->num_edge_cols;

    // Merge all edges into one colour
    if(graph->num_edge_cols == 1) {
      for(i = 0; i < num_cols_loaded; i++)
        col_edges[0] |= edges[i];
    }
    else if(graph->num_edge_cols > 0) {
      for(i = 0; i < num_cols_loaded; i++) {
        col = prefs->into_colour+i;
        col_edges[col] |= edges[i];
      }
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
    graph_header_dealloc(header);

  return num_cols_loaded;
}

//
// Merging, filtering, combining binary graph files
//

// Filter a binary against the graph
// (only keep nodes and edges that are in the graph)
// Same functionality as graph_files_merge, but faster if dealing with only one
// input file
// reads in and dumps one kmer at a time
size_t binary_filter_graph(const char *out_ctx_path, char *in_ctx_path,
                           boolean flatten, const dBGraph *db_graph,
                           const char *intersect_gname)
{
  assert(db_graph->col_edges != NULL);

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

  setvbuf(in, NULL, _IOFBF, CTX_BUF_SIZE);

  GraphFileHeader inheader = {.capacity = 0};
  graph_file_read_header(in, &inheader, true, in_ctx_path);

  if(split != NULL) *split = ':';
  uint32_t num_of_cols = graph_file_get_ncols(in_ctx_path, inheader.num_of_cols-1);
  uint32_t load_colours[num_of_cols];
  graph_file_parse_colours(in_ctx_path, load_colours, inheader.num_of_cols-1);

  timestamp(ctx_msg_out);
  message(" Loading binary %s", in_ctx_path);
  if(split != NULL) {
    message(" with colour filter: %u", load_colours[0]);
    for(i = 1; i < num_of_cols; i++) message(",%u", load_colours[i]);
  }
  if(flatten || num_of_cols == 1) message(" into colour %u\n", offset);
  else message(" into colours %u-%u\n", offset, offset+num_of_cols-1);

  GraphFileHeader outheader = {.capacity = 0};
  graph_header_cpy(&outheader, &inheader);
  outheader.num_of_cols = offset+num_of_cols;
  graph_header_alloc(&outheader, outheader.num_of_cols);

  for(i = 0; i < num_of_cols; i++) {
    GraphInfo *ginfo_out = outheader.ginfo+offset+i;
    graph_info_merge(ginfo_out, inheader.ginfo + load_colours[i]);
    if(intersect_gname != NULL)
      graph_info_set_intersect(&ginfo_out->cleaning, intersect_gname);
  }

  graph_write_header(out, &outheader);

  BinaryKmer bkmer;
  Covg kmercovgs[inheader.num_of_cols], covgs[offset+num_of_cols];
  Edges kmeredges[inheader.num_of_cols], edges[offset+num_of_cols];
  memset(covgs, 0, (offset+num_of_cols)*sizeof(Covg));
  memset(edges, 0, (offset+num_of_cols)*sizeof(Edges));

  while(graph_file_read_kmer(in, &inheader, in_ctx_path, bkmer.b, kmercovgs, kmeredges))
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
      union_edges = db_node_col_edges_union(db_graph, node);

      for(i = 0; i < num_of_cols; i++) edges[offset+i] &= union_edges;
      graph_write_kmer(out, &outheader, bkmer.b, covgs, edges);
      nodes_dumped++;
    }
  }

  fclose(in);
  fclose(out);

  graph_write_status(db_graph->ht.unique_kmers, outheader.num_of_cols,
                     out_ctx_path, CTX_GRAPH_FILEFORMAT);

  graph_header_dealloc(&inheader);
  graph_header_dealloc(&outheader);

  return nodes_dumped;
}



// ctx_num_cols and ctx_max_cols are the numbers returned from binary_probe
// if merge: pool colour 0 from each binary into colour 0, 1 -> 1 etc.
// if flatten: pool all colours into colour 0
// if intersect_gname != NULL: only load kmers that are already in the hash table
//    and use string as name for cleaning against
size_t graph_files_merge(char *out_ctx_path, char **binary_paths,
                         size_t num_binaries,
                         uint32_t ctx_num_cols[num_binaries],
                         uint32_t ctx_max_num_cols[num_binaries],
                         boolean merge, boolean flatten,
                         const char *intersect_gname, dBGraph *db_graph)
{
  assert(!merge || !flatten);
  assert(intersect_gname || db_graph->ht.unique_kmers == 0);

  if(num_binaries == 1 && intersect_gname) {
    return binary_filter_graph(out_ctx_path, binary_paths[0], flatten,
                               db_graph, intersect_gname);
  }

  uint32_t i, offsets[num_binaries];
  uint32_t max_num_cols = 0, sum_cols = 0;

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

    max_num_cols = MAX2(offsets[i] + ctx_num_cols[i], max_num_cols);
    sum_cols += ctx_num_cols[i];
  }

  uint32_t output_colours;

  if(flatten) output_colours = 1;
  else if(merge) output_colours = max_num_cols;
  else output_colours = sum_cols;

  // printf("output_colours: %u sum_cols %u\n", output_colours, sum_cols);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs
    = {.into_colour = 0, .db_graph = db_graph,
       .merge_colours = true,
       .boolean_covgs = false,
       .must_exist_in_graph = intersect_gname != NULL, .empty_colours = false};
  //

  if(output_colours == 1)
  {
    // e.g. flatten
    for(i = 0; i < num_binaries; i++)
      graph_load(binary_paths[i], &prefs, stats, NULL);

    if(intersect_gname) {
      graph_info_set_intersect(&db_graph->ginfo[0].cleaning, intersect_gname);
    }

    hash_table_print_stats(&db_graph->ht);
    graph_file_save(out_ctx_path, db_graph, CTX_GRAPH_FILEFORMAT, NULL, 0, 1);
  }
  else
  {
    uint32_t load_colours[num_binaries][max_num_cols];
    for(i = 0; i < num_binaries; i++)
      graph_file_parse_colours(binary_paths[i], load_colours[i], ctx_max_num_cols[i]);

    // Construct binary header
    GraphFileHeader tmpheader = {.capacity = 0};
    GraphFileHeader output_header = {.version = CTX_GRAPH_FILEFORMAT,
                                     .kmer_size = db_graph->kmer_size,
                                     .num_of_bitfields = NUM_BKMER_WORDS,
                                     .num_of_cols = output_colours,
                                     .num_of_kmers = db_graph->ht.unique_kmers,
                                     .capacity = 0};

    graph_header_alloc(&tmpheader, max_num_cols);
    graph_header_alloc(&output_header, output_colours);

    Colour j, output_colour = 0;
    for(i = 0; i < num_binaries; i++)
    {
      graph_load(binary_paths[i], &prefs, stats, &tmpheader);
      if(merge) output_colour = 0;
      output_colour += offsets[i];
      for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
      {
        graph_info_merge(output_header.ginfo + output_colour,
                         tmpheader.ginfo + load_colours[i][j]);
      }
    }

    if(intersect_gname) {
      for(i = 0; i < output_colours; i++) {
        db_graph->ginfo[i].cleaning.is_graph_intersection = true;
        strbuf_set(&db_graph->ginfo[i].cleaning.intersection_name,
                   intersect_gname);
      }
    }

    FILE *fh = fopen(out_ctx_path, "w");
    if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

    size_t header_size = graph_write_header(fh, &output_header);

    // Free header resources
    graph_header_dealloc(&tmpheader);
    graph_header_dealloc(&output_header);

    // print file outline
    status("Generated merged hash table\n");
    hash_table_print_stats(&db_graph->ht);
    graph_write_empty(db_graph, fh, output_colours);

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
            graph_load_colour(binary_paths[i], &prefs, stats,
                               load_colours[i][ctx_col]);
            data_loaded_in_col = true;
          }
        }
        if(data_loaded_in_col) {
          status("Dumping into colour %zu...\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          graph_file_write_colour(db_graph, 0, output_colour, output_colours, fh);
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

          graph_load_colour(binary_paths[i], &prefs, stats,
                             load_colours[i][j]);

          status("Dumping into colour %zu...\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          graph_file_write_colour(db_graph, 0, output_colour, output_colours, fh);
        }
      }
    }

    fclose(fh);

    graph_write_status(db_graph->ht.unique_kmers, output_colours,
                       out_ctx_path, CTX_GRAPH_FILEFORMAT);
  }

  seq_loading_stats_free(stats);

  return db_graph->ht.unique_kmers;
}

