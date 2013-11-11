#include "global.h"
#include "range.h"
#include "graph_file_filter.h"
#include "graph_format.h"
#include "file_util.h"

const GraphFileHeader INIT_GRAPH_FILE_HDR = INIT_GRAPH_FILE_HDR_MACRO;
const GraphFileReader INIT_GRAPH_READER = INIT_GRAPH_READER_MACRO;

// Get pointers to start and end of actual path
// (\d+:)?path.ctx(:\d+(-\d+)?(,\d+(-\d+)?)*)?
void graph_file_filter_deconstruct(char *path, char **start, char **end)
{
  char *ptr, *c;
  ptr = *start = path;
  for(ptr = path; *ptr >= '0' && *ptr <= '9'; ptr++);
  if(ptr > path && *ptr == ':') { ptr++; *start = ptr; }
  // Count backwards to match /:[-,0123456789]*$/
  c = *end = path + strlen(path);
  while(c > (*start)+1) {
    c--;
    if(*c == ':') { *end = c; break; }
    else if(!(*c == ',' || *c == '-' || (*c >= '0' && *c <= '9'))) break;
  }
}

static void graph_file_capacity(GraphFileReader *file, size_t ncolscap)
{
  if(ncolscap == 0) return;
  else if(file->ncolscap == 0) {
    file->cols = malloc2(ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
  else if(file->ncolscap < ncolscap) {
    file->cols = realloc2(file->cols, ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, char *path, boolean fatal)
{
  size_t i;
  char *path_start, *path_end, path_lchar;

  if(file->fh != NULL) graph_file_close(file);

  graph_file_filter_deconstruct(path, &path_start, &path_end);
  file->intocol = (path_start == path ? 0 : atoi(path));

  path_lchar = *path_end;
  *path_end = '\0';

  file->file_size = get_file_size(path_start);
  if(file->file_size == -1) die("Cannot get file size: %s", path_start);

  if((file->fh = fopen(path_start, "r")) == NULL) {
    if(fatal) die("Cannot load binary file: %s", path_start);
    else return 0;
  }
  setvbuf(file->fh, NULL, _IOFBF, CTX_BUF_SIZE);

  file->hdr_size = graph_file_read_header(file->fh, &file->hdr, fatal, path_start);
  if(file->hdr_size == -1) return -1;

  // status("File: %s; hsize: %zu fsize: %zu", path_start,
  //        (size_t)file->file_size, (size_t)file->hdr_size);
  // graph_header_print(&file->hdr);

  // Get number of kmers
  size_t bytes_per_kmer = sizeof(BinaryKmer) +
                          file->hdr.num_of_cols * (sizeof(Covg) + sizeof(Edges));
  size_t bytes_remaining = file->file_size - file->hdr_size;
  size_t nkmers = (bytes_remaining / bytes_per_kmer);

  if(file->hdr.version > 6 && file->hdr.num_of_kmers != nkmers) {
    warn("File size and number of kmers do not match: %s [bytes per kmer: %zu "
         "remaining: %zu; fsize: %zu; header: %zu; expect: %zu; got: %zu]",
         path_start, bytes_per_kmer, bytes_remaining,
         (size_t)file->file_size, (size_t)file->hdr_size,
         (size_t)file->hdr.num_of_kmers, nkmers);
  }

  if(bytes_remaining % bytes_per_kmer != 0) {
    warn("Truncated graph file: %s [bytes per kmer: %zu "
         "remaining: %zu; fsize: %zu; header: %zu; nkmers: %zu]",
         path_start, bytes_per_kmer, bytes_remaining,
         (size_t)file->file_size, (size_t)file->hdr_size, nkmers);
  }
  file->hdr.num_of_kmers = nkmers;

  // Success
  if(file->path.buff == NULL) strbuf_alloc(&file->path, 1024);
  strbuf_set(&file->path, path_start);
  *path_end = path_lchar;

  if(*path_end == ':') {
    file->ncols = range_get_num(path_end+1, file->hdr.num_of_cols-1);
    graph_file_capacity(file, file->ncols);
    range_parse_array(path_end+1, file->cols, file->hdr.num_of_cols-1);
  }
  else {
    file->ncols = file->hdr.num_of_cols;
    graph_file_capacity(file, file->ncols);
    for(i = 0; i < file->hdr.num_of_cols; i++) file->cols[i] = i;
  }

  return 1;
}

// Close file
void graph_file_close(GraphFileReader *file)
{
  if(file->fh != NULL) { fclose(file->fh); file->fh = NULL; }
}

void graph_file_dealloc(GraphFileReader *file)
{
  graph_file_close(file);
  graph_header_dealloc(&file->hdr);
  if(file->ncolscap > 0) { free(file->cols); file->ncolscap = 0; }
  if(file->path.buff != NULL) { strbuf_dealloc(&file->path); }
  if(file->fh != NULL) { fclose(file->fh); }
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
boolean graph_file_read(const GraphFileReader *file,
                        BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  // status("Header colours: %u", file->hdr.num_of_cols);
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i;

  if(!graph_file_read_kmer(file->fh, &file->hdr, file->path.buff,
                           bkmer->b, kmercovgs, kmeredges)) return 0;

  // covgs += file->intocol;
  // edges += file->intocol;

  if(file->flatten) {
    covgs[0] = 0;
    edges[0] = 0;
    for(i = 0; i < file->ncols; i++) {
      covgs[0] += kmercovgs[file->cols[i]];
      edges[0] |= kmeredges[file->cols[i]];
    }
  }
  else {
    for(i = 0; i < file->ncols; i++) {
      covgs[i] = kmercovgs[file->cols[i]];
      edges[i] = kmeredges[file->cols[i]];
    }
  }

  return 1;
}

// Return true if all colours are being loaded once in their original order
boolean graph_file_no_filter(const GraphFileReader *file)
{
  size_t i;
  if(file->ncols != file->hdr.num_of_cols) return false;
  for(i = 0; i < file->ncols && file->cols[i] == i; i++);
  return (i == file->ncols);
}

// Print file filter description
void graph_file_status(const GraphFileReader *file)
{
  size_t i;
  timestamp(ctx_msg_out);
  message(" Loading graph file %s", file->path.buff);
  if(!graph_file_no_filter(file)) {
    message(" with colour filter: %u", file->cols[0]);
    for(i = 1; i < file->ncols; i++) message(",%u", file->cols[i]);
  }
  size_t into_ncols = graph_file_outncols(file);
  if(into_ncols == 1)
    message(" into colour %u\n", file->intocol);
  else
    message(" into colours %u-%zu\n", file->intocol, file->intocol+into_ncols-1);
}
