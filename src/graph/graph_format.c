#include "global.h"
#include "graph_format.h"

void graph_header_capacity(GraphFileHeader *h, size_t num_of_cols)
{
  size_t i;

  if(num_of_cols > h->capacity) {
    h->ginfo = ctx_recallocarray(h->ginfo, h->capacity, num_of_cols, sizeof(GraphInfo));
    for(i = h->capacity; i < num_of_cols; i++)
      graph_info_alloc(&h->ginfo[i]);
    h->capacity = num_of_cols;
  }
}

void graph_header_dealloc(GraphFileHeader *h)
{
  size_t i;
  for(i = 0; i < h->capacity; i++)
    graph_info_dealloc(&h->ginfo[i]);
  ctx_free(h->ginfo);
  memset(h, 0, sizeof(*h));
}

void graph_header_print(const GraphFileHeader *header)
{
  printf("HEADER\n");
  printf("  version: %u\n", header->version);
  printf("  kmer_size: %u\n", header->kmer_size);
  printf("  num_of_bitfields: %u\n", header->num_of_bitfields);
  printf("  num_of_cols: %u\n", header->num_of_cols);
  printf("  [capacity: %zu]\n", header->capacity);
}
