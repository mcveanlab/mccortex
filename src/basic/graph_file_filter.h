#ifndef GRAPH_FILE_FILTER_H_
#define GRAPH_FILE_FILTER_H_

#include "graph_typedef.h"

typedef struct
{
  uint32_t ncols, ncolcap, *cols, intocol;
  boolean flatten;
} GraphFileFilter;

#define INIT_GRAPH_FILE_FILTER {.ncolcap = 0}

#define graph_file_filter_outncols(gff) ((gff)->flatten ? 1 : (gff)->ncols)

void graph_file_filter_alloc(GraphFileFilter *gff, size_t ncolcap);
void graph_file_filter_dealloc(GraphFileFilter *gff);

void graph_file_filter_parse(const char *path, GraphFileFilter *gff,
                             size_t file_ncols);

#endif /* GRAPH_FILE_FILTER_H_ */
