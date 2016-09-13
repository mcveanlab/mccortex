#include "global.h"
#include "graph_file_reader.h"
#include "db_node.h"
#include "cmd.h"
#include "file_util.h"

int graph_file_fseek(GraphFileReader *file, off_t offset, int whence)
{
  if(file_filter_isstdin(&file->fltr)) die("Cannot fseek on STDIN");
  if(graph_file_is_buffered(file))
    return fseek_buf(file->fh, offset, whence, &file->strm);
  else
    return fseek(file->fh, offset, whence);
}

off_t graph_file_ftell(GraphFileReader *file)
{
  if(graph_file_is_buffered(file))
    return ftell_buf(file->fh, &file->strm);
  else
    return ftell(file->fh);
}

size_t gfr_fread_bytes(GraphFileReader *file, void *ptr, size_t size)
{
  size_t n;
  if(graph_file_is_buffered(file))
    n = fread_buf(file->fh, ptr, size, &file->strm);
  else
    n = fread2(file->fh, ptr, size);
  // check for error
  if(ferror(file->fh))
    die("File error: %s [%s]", strerror(errno), file_filter_path(&file->fltr));
  return n;
}

// Read an element from the graph file
#define _gfread(gfile,ptr,size,desc) \
do { \
  size_t _n = gfr_fread_bytes(gfile, ptr, size); \
  const char *_path = file_filter_path(&(gfile)->fltr); \
  if(_n != (size)) { \
    die("Couldn't read '%s': expected %zu; recieved: %zu; [file: %s]\n",\
        (desc), (size_t)(size), _n, _path); \
  } \
} while(0)

void graph_file_merge_header(GraphFileHeader *hdr, const GraphFileReader *file)
{
  size_t i, fromcol, intocol;
  hdr->version = CTX_GRAPH_FILEFORMAT;
  hdr->num_of_bitfields = file->hdr.num_of_bitfields;
  hdr->kmer_size = file->hdr.kmer_size;
  hdr->num_of_cols = MAX2(hdr->num_of_cols, file_filter_into_ncols(&file->fltr));
  graph_header_capacity(hdr, hdr->num_of_cols);

  for(i = 0; i < file_filter_num(&file->fltr); i++) {
    fromcol = file_filter_fromcol(&file->fltr, i);
    intocol = file_filter_intocol(&file->fltr, i);
    graph_info_merge(hdr->ginfo + intocol, file->hdr.ginfo + fromcol);
  }
}

// Return number of bytes read or die() with error
size_t graph_file_read_header(GraphFileReader *file)
{
  size_t i;
  int bytes_read = 0;
  GraphFileHeader *h = &file->hdr;
  const char *path = file_filter_path(&file->fltr);

  char magic_word[7];
  magic_word[6] = '\0';

  _gfread(file, magic_word, strlen("CORTEX"), "Magic word");
  if(strcmp(magic_word, "CORTEX") != 0) {
    die("Magic word doesn't match '%s' (start): %s", "CORTEX", path);
  }
  bytes_read += strlen("CORTEX");

  // Read version number, kmer_size, num bitfields, num colours
  _gfread(file, &h->version, sizeof(uint32_t), "graph version");
  _gfread(file, &h->kmer_size, sizeof(uint32_t), "kmer size");
  _gfread(file, &h->num_of_bitfields, sizeof(uint32_t), "num of bitfields");
  _gfread(file, &h->num_of_cols, sizeof(uint32_t), "number of colours");
  bytes_read += 4*sizeof(uint32_t);

  // Checks
  if(h->version > 7 || h->version < 4)
  {
    die("Sorry, we only support graph file versions 4, 5, 6 & 7 "
        "[version: %u; path: %s]\n", h->version, path);
  }

  if(h->kmer_size % 2 == 0)
  {
    die("kmer size is not an odd number [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3)
  {
    die("kmer size is less than three [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_bitfields * 32 < h->kmer_size)
  {
    die("Not enough bitfields for kmer size "
        "[kmer_size: %u; bitfields: %u; path: %s]\n",
        h->kmer_size, h->num_of_bitfields, path);
  }

  if((h->num_of_bitfields-1)*32 >= h->kmer_size) {
    die("using more than the minimum number of bitfields [path: %s]\n", path);
  }

  if(h->num_of_cols == 0)
    die("number of colours is zero [path: %s]\n", path);
  if(h->num_of_cols > 10000)
    die("Very high number of colours: %zu [path: %s]", (size_t)h->num_of_cols, path);

  // graph_header_capacity will only alloc or realloc if it needs to
  graph_header_capacity(h, h->num_of_cols);

  // Assume to be graph file now, any error is therefore fatal

  for(i = 0; i < h->num_of_cols; i++) {
    _gfread(file, &h->ginfo[i].mean_read_length, sizeof(uint32_t),
             "mean read length for each colour");
  }

  for(i = 0; i < h->num_of_cols; i++) {
    _gfread(file, &h->ginfo[i].total_sequence, sizeof(uint64_t),
             "total sequance loaded for each colour");
  }

  bytes_read += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint64_t));

  if(h->version >= 6)
  {
    for(i = 0; i < h->num_of_cols; i++)
    {
      uint32_t len;
      _gfread(file, &len, sizeof(uint32_t), "sample name length");
      if(len > 10000) die("Very big sample name. Length: %u", len);

      StrBuf *sbuf = &h->ginfo[i].sample_name;
      strbuf_ensure_capacity(sbuf, len);

      _gfread(file, sbuf->b, len, "sample name");
      bytes_read += sizeof(uint32_t) + len;

      sbuf->b[len] = '\0';

      size_t len2 = strlen(sbuf->b);
      sbuf->end = len2;
      sbuf->b[len2] = '\0';

      // Check sample length is as long as we were told
      if(len != len2)
      {
        warn("Sample %zu name has length %u but is only %zu chars long "
             "(premature '\\0') [path: %s]\n", i, len, len2, path);
      }
    }

    for(i = 0; i < h->num_of_cols; i++) {
      _gfread(file, &h->ginfo[i].seq_err, sizeof(long double), "seq error rates");
    }

    bytes_read += sizeof(long double) * h->num_of_cols;

    for(i = 0; i < h->num_of_cols; i++)
    {
      ErrorCleaning *err_cleaning = &h->ginfo[i].cleaning;

      _gfread(file, &(err_cleaning->cleaned_tips),
              sizeof(uint8_t), "tip cleaning");
      _gfread(file, &(err_cleaning->cleaned_unitigs),
              sizeof(uint8_t), "remove low covg unitig");
      _gfread(file, &(err_cleaning->cleaned_kmers),
              sizeof(uint8_t), "remove low covg kmers");
      _gfread(file, &(err_cleaning->is_graph_intersection),
              sizeof(uint8_t), "cleaned against graph");

      uint32_t clean_unitigs_thresh = 0, clean_kmers_thresh = 0;
      _gfread(file, &clean_unitigs_thresh,
              sizeof(uint32_t), "remove low covg unitig threshold");
      _gfread(file, &clean_kmers_thresh,
              sizeof(uint32_t), "remove low covg kmer threshold");

      bytes_read += 4*sizeof(uint8_t) + 2*sizeof(uint32_t);

      // Fix for old versions with negative thresholds
      if(h->version <= 6) {
        if(!err_cleaning->cleaned_unitigs && clean_unitigs_thresh == (uint32_t)-1)
          clean_unitigs_thresh = 0;
        if(!err_cleaning->cleaned_kmers && clean_kmers_thresh == (uint32_t)-1)
          clean_kmers_thresh = 0;
      }

      // Sanity checks
      if(!err_cleaning->cleaned_unitigs && clean_unitigs_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for unitig "
             "when no cleaning was performed [path: %s]", path);
        clean_unitigs_thresh = 0;
      }

      if(!err_cleaning->cleaned_kmers && clean_kmers_thresh > 0)
      {
        warn("Graph header gives cleaning threshold for nodes "
             "when no cleaning was performed [path: %s]", path);
        clean_kmers_thresh = 0;
      }

      err_cleaning->clean_unitigs_thresh = clean_unitigs_thresh;
      err_cleaning->clean_kmers_thresh = clean_kmers_thresh;

      // Read cleaned against name
      uint32_t len;
      _gfread(file, &len, sizeof(uint32_t), "graph name length");
      if(len > 10000) die("Very big sample name. Length: %u", len);

      StrBuf *sbuf = &err_cleaning->intersection_name;
      strbuf_ensure_capacity(sbuf, len);

      _gfread(file, sbuf->b, len, "cleaned against graph name");
      sbuf->b[len] = '\0';

      bytes_read += sizeof(uint32_t) + len;

      size_t len2 = strlen(sbuf->b);
      sbuf->end = len2;
      sbuf->b[len2] = '\0';

      // Check sample length is as long as we were told
      if(len != len2)
      {
        warn("Sample %zu name has length %u but is only %zu chars long "
             "(premature '\\0') [path: %s]\n", i, len, len2, path);
      }
    }
  }

  // Read magic word at the end of header 'CORTEX'
  _gfread(file, magic_word, strlen("CORTEX"), "magic word (end)");
  if(strcmp(magic_word, "CORTEX") != 0)
  {
    die("Magic word doesn't match '%s' (end): '%s' [path: %s]\n",
        "CORTEX", magic_word, path);
  }
  bytes_read += strlen("CORTEX");

  return bytes_read;
}

int graph_file_open(GraphFileReader *file, const char *path)
{
  return graph_file_open2(file, path, "r", true, 0);
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open2(GraphFileReader *file, const char *input, const char *mode,
                     bool usebuf, size_t into_offset)
{
  GraphFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;
  file_filter_open(fltr, input); // calls die() on error
  const char *path = fltr->path.b;

  // Reset reading errors
  file->error_zero_covg = false;
  file->error_missing_covg = false;

  // Stat will fail on streams, so file_size and num_of_kmers with both be -1
  struct stat st;
  file->file_size = -1;
  file->num_of_kmers = -1;

  if(strcmp(input,"-") != 0) {
    if(stat(path, &st) == 0) file->file_size = st.st_size;
    else warn("Couldn't get file size: %s", futil_outpath_str(path));
  }

  file->fh = futil_fopen(path, mode);
  if(usebuf) strm_buf_alloc(&file->strm, ONE_MEGABYTE);
  else memset(&file->strm, 0, sizeof(file->strm));
  file->hdr_size = graph_file_read_header(file);

  file_filter_set_cols(fltr, hdr->num_of_cols, into_offset);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(file->hdr.kmer_size, file->fltr.path.b);

  size_t bytes_per_kmer, bytes_remaining;

  // If reading from STDIN we don't know file size
  if(file->file_size != -1)
  {
    // File header checks
    // Get number of kmers
    bytes_per_kmer = sizeof(BinaryKmer) +
                     hdr->num_of_cols * (sizeof(Covg) + sizeof(Edges));
    bytes_remaining = (size_t)(file->file_size - file->hdr_size);
    file->num_of_kmers = (bytes_remaining / bytes_per_kmer);

    if(bytes_remaining % bytes_per_kmer != 0) {
      warn("Truncated graph file: %s [bytes per kmer: %zu "
           "remaining: %zu; fsize: %zu; header: %zu; nkmers: %zu]",
           path, bytes_per_kmer, bytes_remaining,
           (size_t)file->file_size, (size_t)file->hdr_size,
           (size_t)file->num_of_kmers);
    }
  }

  return 1;
}

// Close file
void graph_file_close(GraphFileReader *file)
{
  strm_buf_dealloc(&file->strm);
  if(file->fh) fclose(file->fh);
  file_filter_close(&file->fltr);
  graph_header_dealloc(&file->hdr);
  memset(file, 0, sizeof(*file));
}

size_t graph_file_read_raw(GraphFileReader *file,
                           BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  GraphFileHeader *h = &file->hdr;
  const char *path = file_filter_path(&file->fltr);

  size_t i;
  int num_bytes_read;
  char kstr[MAX_KMER_SIZE+1];

  num_bytes_read = gfr_fread_bytes(file, bkmer->b, sizeof(BinaryKmer));

  if(num_bytes_read == 0) return 0;
  if(num_bytes_read != (int)(sizeof(uint64_t)*h->num_of_bitfields))
    die("Unexpected end of file: %s", path);

  _gfread(file, covgs, h->num_of_cols * sizeof(uint32_t), "Coverages");
  _gfread(file, edges, h->num_of_cols * sizeof(uint8_t), "Edges");
  num_bytes_read += h->num_of_cols * (sizeof(uint32_t) + sizeof(uint8_t));

  // Check top word of each kmer
  if(binary_kmer_oversized(*bkmer, h->kmer_size))
    die("Oversized kmer in path [kmer: %u]: %s", h->kmer_size, path);

  // Check covg is not 0 for all colours
  for(i = 0; i < h->num_of_cols && covgs[i] == 0; i++) {}
  if(i == h->num_of_cols && !file->error_zero_covg) {
    binary_kmer_to_str(*bkmer, h->kmer_size, kstr);
    warn("Kmer has zero covg in all colours [kmer: %s; path: %s]", kstr, path);
    file->error_zero_covg = true;
  }

  // Check edges => coverage
  for(i = 0; i < h->num_of_cols && (!edges[i] || covgs[i]); i++) {}
  if(i < h->num_of_cols && !file->error_missing_covg) {
    binary_kmer_to_str(*bkmer, h->kmer_size, kstr);
    warn("Kmer has edges but no coverage [kmer: %s; path: %s]", kstr, path);
    file->error_missing_covg = true;
  }

  return num_bytes_read;
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// be sure to zero covgs, edges before reading in
bool graph_file_read(GraphFileReader *file,
                     BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  // status("Header colours: %u", file->hdr.num_of_cols);
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i, from, into;
  const FileFilter *fltr = &file->fltr;

  if(!graph_file_read_raw(file, bkmer, kmercovgs, kmeredges)) return false;

  for(i = 0; i < file_filter_num(fltr); i++) {
    from = file_filter_fromcol(fltr, i);
    into = file_filter_intocol(fltr, i);
    covgs[into] = SAFE_ADD_COVG(covgs[into], kmercovgs[from]);
    edges[into] |= kmeredges[from];
  }

  return true;
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// @ncols is file_filter_into_ncols(&file->fltr)
bool graph_file_read_reset(GraphFileReader *file,
                           BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  size_t ncols = file_filter_into_ncols(&file->fltr);
  memset(covgs, 0, ncols*sizeof(Covg));
  memset(edges, 0, ncols*sizeof(Edges));
  return graph_file_read(file, bkmer, covgs, edges);
}

// Read coverages and edges only (no kmer)
void graph_file_read_covgs_edges(GraphFileReader *f, Covg *covgs, Edges *edges)
{
  size_t from, into, i;
  const FileFilter *fltr = &f->fltr;
  Covg allcovgs[f->hdr.num_of_cols];
  Edges alledges[f->hdr.num_of_cols];
  gfr_fread_bytes(f, allcovgs, sizeof(allcovgs));
  gfr_fread_bytes(f, alledges, sizeof(alledges));
  if(covgs) {
    memset(covgs, 0, file_filter_into_ncols(fltr) * sizeof(Covg));
    for(i = 0; i < file_filter_num(fltr); i++) {
      from = file_filter_fromcol(fltr, i);
      into = file_filter_intocol(fltr, i);
      covgs[into] = SAFE_ADD_COVG(covgs[into], allcovgs[from]);
    }
  }
  if(edges) {
    memset(edges, 0, file_filter_into_ncols(fltr) * sizeof(Edges));
    for(i = 0; i < file_filter_num(fltr); i++) {
      from = file_filter_fromcol(fltr,i);
      into = file_filter_intocol(fltr, i);
      edges[into] |= alledges[from];
    }
  }
}


// Returns true if one or more files passed loads data into colour
bool graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                 size_t num_files)
{
  size_t i;
  for(i = 0; i < num_files; i++) {
    if(file_filter_iscolloaded(&files[i].fltr, colour))
      return true;
  }
  return false;
}

// if one of the files is reading from stdin, sum_kmers_ptr is set to 0
// `max_cols_ptr` is used to return the most colours being loaded from a single file
// returns the number of colours being loaded in total
size_t graph_files_open(char **graph_paths,
                        GraphFileReader *gfiles, size_t num_gfiles,
                        size_t *max_kmers_ptr, size_t *sum_kmers_ptr)
{
  size_t i, ctx_max_kmers = 0, ctx_sum_kmers = 0;
  bool ctx_uses_stdin = false;
  size_t ncols = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    memset(&gfiles[i], 0, sizeof(GraphFileReader));
    graph_file_open2(&gfiles[i], graph_paths[i], "r", true, ncols);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                      gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size);
    }

    ncols = MAX2(ncols, file_filter_into_ncols(&gfiles[i].fltr));

    ctx_max_kmers = MAX2(ctx_max_kmers, graph_file_nkmers(&gfiles[i]));
    ctx_sum_kmers += graph_file_nkmers(&gfiles[i]);
    ctx_uses_stdin |= file_filter_isstdin(&gfiles[i].fltr);
  }

  if(ctx_uses_stdin) ctx_sum_kmers = SIZE_MAX;

  *max_kmers_ptr = ctx_max_kmers;
  *sum_kmers_ptr = ctx_sum_kmers;

  return ncols;
}
