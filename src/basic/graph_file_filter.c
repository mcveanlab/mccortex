#include "global.h"
#include "range.h"
#include "graph_file_filter.h"
#include "graph_format.h"
#include "file_util.h"

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

// // (\d+:)?path.ctx(:\d+)?
// void graph_file_filter_alloc(GraphFileFilter *gff, char *path, GraphFileHeader *hdr)
// {
//   char *start = path, *end;
//   size_t i;

//   graph_file_filter_deconstruct(path, &start, &end);
//   gff->intocol = (start == path ? 0 : atoi(path));

//   if(*end == ':') {
//     gff->ncols = range_get_num(end+1, hdr->num_of_cols-1);
//     gff->cols = malloc2(gff->ncols * sizeof(*(gff->cols)));
//     range_parse_array(end+1, gff->cols, hdr->num_of_cols-1);
//     // printf("%u cols:", gff->ncols);
//     // for(i=0; i<gff->ncols; i++) {printf(" %u", gff->cols[i]);}
//     // printf("\n");
//   }
//   else {
//     gff->ncols = hdr->num_of_cols;
//     gff->cols = malloc2(gff->ncols * sizeof(*(gff->cols)));
//     for(i = 0; i < hdr->num_of_cols; i++) gff->cols[i] = i;
//   }

//   char end_store = *end;
//   *end = '\0';
//   gff->path = strdup(start);
//   *end = end_store;

//   gff->hdr = hdr;
//   gff->fh = NULL;
//   gff->flatten = false;
// }

// void graph_file_filter_dealloc(GraphFileFilter *gff)
// {
//   free(gff->path);
//   free(gff->cols);
// }

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, char *path, boolean fatal)
{
  size_t i;
  char *path_start, *path_end, path_lchar;

  graph_file_filter_deconstruct(path, &path_start, &path_end);
  file->intocol = (path_start == path ? 0 : atoi(path));

  path_lchar = *path_end;
  *path_end = '\0';

  file->file_size = get_file_size(path_start);
  if(file->file_size == -1) die("Cannot get file size: %s", path_start);

  if((file->fh = fopen(path_start, "r")) == NULL) return 0;
  setvbuf(file->fh, NULL, _IOFBF, CTX_BUF_SIZE);

  file->hdr_size = graph_file_read_header(file->fh, &file->hdr, fatal, path_start);
  if(file->hdr_size == -1) return -1;

  // Success
  file->path = strdup(path_start);
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
  graph_header_dealloc(&file->hdr);
  if(file->ncolscap > 0) { free(file->cols); file->ncolscap = 0; }
  if(file->path != NULL) { free(file->path); }
  if(file->fh != NULL) { fclose(file->fh); }
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
boolean graph_file_read(const GraphFileReader *file,
                        BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i;

  if(!graph_file_read_kmer(file->fh, &file->hdr, file->path,
                           bkmer->b, kmercovgs, kmeredges)) return 0;

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
