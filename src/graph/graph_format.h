#ifndef GRAPH_FORMAT_H_
#define GRAPH_FORMAT_H_

// graph file format version
#define CTX_GRAPH_FILEFORMAT 6

#include "graph_info.h"

// Graph (.ctx)
typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols;
  GraphInfo *ginfo; // Cleaning info etc for each colour
  size_t capacity;
} GraphFileHeader;

void graph_header_capacity(GraphFileHeader *header, size_t num_of_cols);
void graph_header_dealloc(GraphFileHeader *header);
void graph_header_print(const GraphFileHeader *header);

static inline void graph_header_free(GraphFileHeader *hdr) {
  graph_header_dealloc(hdr);
  ctx_free(hdr);
}

#endif /* GRAPH_FORMAT_H_ */
