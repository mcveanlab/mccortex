#include "global.h"
#include "range.h"
#include "graph_file_filter.h"

void graph_file_filter_alloc(GraphFileFilter *gff, size_t ncolcap)
{
  gff->cols = malloc2(ncolcap * sizeof(*(gff->cols)));
  gff->ncols = 0;
  gff->ncolcap = ncolcap;
  gff->intocol = 0;
  gff->flatten = false;
}

void graph_file_filter_dealloc(GraphFileFilter *gff)
{
  free(gff->cols);
}

void graph_file_filter_capacity(GraphFileFilter *gff, size_t ncap)
{
  if(gff->ncolcap == 0) {
    gff->ncolcap = ROUNDUP2POW(ncap);
    gff->cols = malloc2(gff->ncolcap * sizeof(*(gff->cols)));
  }
  else if(gff->ncolcap < ncap) {
    gff->ncolcap = ROUNDUP2POW(ncap);
    gff->cols = realloc2(gff->cols, gff->ncolcap * sizeof(*(gff->cols)));
  }
}

// Get pointers to start and end of actual path
// (\d+:)?path.ctx(:\d+)?
void graph_file_filter_deconstruct(const char *path,
                                   const char **start, const char **end)
{
  const char *ptr = path;
  for(ptr = path; *ptr >= '0' && *ptr <= '9'; ptr++);
  if(ptr > path && *ptr == ':') { *start = ptr+1; ptr++; }
  ptr = strchr(ptr, ':');
  *end = (ptr == NULL ? path + strlen(path) : ptr);
}

// (\d+:)?path.ctx(:\d+)?
void graph_file_filter_parse(const char *path, GraphFileFilter *gff,
                             size_t file_ncols)
{
  const char *start, *end;
  size_t i;

  graph_file_filter_deconstruct(path, &start, &end);
  gff->intocol = start == path ? 0 : atoi(path);

  if(*end == ':') {
    gff->ncols = range_get_num(end+1, file_ncols-1);
    graph_file_filter_capacity(gff, gff->ncols);
    range_parse_array(end+1, gff->cols, file_ncols-1);
  }
  else {
    graph_file_filter_capacity(gff, gff->ncols);
    for(i = 0; i < file_ncols; i++) gff->cols[i] = i;
  }
}

// Returns true if kmers are non-zero in covg or edges
boolean graph_file_filter_load(const Covg *kmercovgs, const Edges *kmeredges,
                               Covg *covgs, Edges *edges,
                               const GraphFileFilter *gff)
{
  size_t i; Covg nonzero = 0;
  if(gff->flatten) {
    covgs[0] = 0;
    edges[0] = 0;
    for(i = 0; i < gff->ncols; i++) {
      covgs[0] += kmercovgs[gff->cols[i]];
      edges[0] |= kmeredges[gff->cols[i]];
      nonzero |= covgs[0] | edges[0];
    }
  }
  else {
    for(i = 0; i < gff->ncols; i++) {
      covgs[i] = kmercovgs[gff->cols[i]];
      edges[i] = kmeredges[gff->cols[i]];
      nonzero |= covgs[i] | edges[i];
    }
  }
  return (nonzero != 0);
}
