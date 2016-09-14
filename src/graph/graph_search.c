#include "global.h"
#include "graph_search.h"

struct GraphFileSearch {
  GraphFileReader *file;
  size_t nkmers, ncols;
  BinaryKmer *index;
  size_t blocksize,nblocks;
};

// #define INDEX_SIZE 4 /* debugging */
#define INDEX_SIZE 4*1024*1024 /* 4M */

GraphFileSearch *graph_search_new(GraphFileReader *file)
{
  if(file->num_of_kmers < 0) {
    warn("Cannot open GraphFileSearch with file stream");
    return NULL;
  }
  size_t i;
  GraphFileSearch *gs = ctx_calloc(sizeof(GraphFileSearch), 1);
  gs->file = file;
  gs->nkmers = file->num_of_kmers;
  gs->ncols = file->fltr.srcncols;
  gs->nblocks = MIN2(gs->nkmers, INDEX_SIZE);
  gs->blocksize = gs->nkmers / gs->nblocks;
  gs->nblocks = (gs->nkmers+gs->blocksize-1) / gs->blocksize;
  gs->index = ctx_calloc(sizeof(BinaryKmer), gs->nblocks+1); // sentinel
  memset(gs->index[gs->nblocks].b,0xff,BKMER_BYTES); // sentinel kmer
  status("[graph_search] on-disk-graph %zu cols %zu blocks %zu bsize %zu kmers"
         " building...", gs->ncols, gs->nblocks, gs->blocksize, gs->nkmers);
  graph_file_set_buffered(file, false); // Turn OFF buffered input
  for(i = 0; i < gs->nblocks; i++) {
    graph_file_fseek(file, graph_file_offset(file, i*gs->blocksize), SEEK_SET);
    graph_file_fread(file, &gs->index[i], sizeof(BinaryKmer));
  }
  graph_file_set_buffered(file, true); // Turn ON buffered input
  // check file is sorted
  for(i = 0; i+1 < gs->nblocks; i++)
    if(!binary_kmer_lt(gs->index[i],gs->index[i+1]))
      die("File is not sorted: %s", file_filter_path(&file->fltr));
  status("[graph_search] Index built.");
  return gs;
}

// We don't close the file
void graph_search_destroy(GraphFileSearch *gs)
{
  ctx_free(gs->index);
  ctx_free(gs);
}

// bkmers[n] must be a sentinel kmer (i.e. MAX_KMER)
static inline int binary_search_index(BinaryKmer bkey,
                                      const BinaryKmer *bkmers, size_t n)
{
  int l = 0, r = n, mid;
  while(l < r) {
    mid = (l+r)/2;
    if(binary_kmer_le(bkmers[mid],bkey)) {
      if(binary_kmer_lt(bkey,bkmers[mid+1])) return mid;
      else l = mid+1;
    }
    else r = mid;
  }
  return -1;
}

static inline bool search_file_sec(GraphFileSearch *gs, BinaryKmer bkey,
                                   size_t start, size_t end)
{
  const size_t hdrsize = gs->file->hdr_size;
  const size_t covg_edges_size = gs->ncols*(sizeof(Covg)+sizeof(Edges));
  const size_t entsize = sizeof(BinaryKmer)+covg_edges_size;
  const size_t max_block = 512; // linear search over small blocks
  size_t mid;
  BinaryKmer bmid;
  // Binary search
  while(start + max_block < end) {
    mid = (start+end) / 2;
    graph_file_fseek(gs->file, hdrsize+entsize*mid, SEEK_SET);
    graph_file_fread(gs->file, bmid.b, sizeof(BinaryKmer));
    if(binary_kmer_eq(bkey,bmid)) { fprintf(stderr,"index hit\n"); return true; }
    if(binary_kmer_lt(bkey,bmid)) end = mid;
    else start = mid + 1;
  }
  // Linear search
  // char bstr[MAX_KMER_SIZE+1];
  // const size_t kmer_size = gs->file->hdr.kmer_size;
  // fprintf(stderr, "Linear search over %zu-%zu\n", start, end);
  graph_file_fseek(gs->file, hdrsize+entsize*start, SEEK_SET);
  for(mid = start; mid < end; mid++) {
    graph_file_fread(gs->file, bmid.b, sizeof(BinaryKmer));
    // fprintf(stderr, "Read %zu: %s\n", mid, binary_kmer_to_str(bmid, kmer_size, bstr));
    if(binary_kmer_eq(bkey,bmid)) return true;
    if(binary_kmer_lt(bkey,bmid)) return false;
    graph_file_fseek(gs->file, covg_edges_size, SEEK_CUR); // skip over covgs+edges
  }
  return false;
}

bool graph_search_find(GraphFileSearch *gs, BinaryKmer bkey,
                            Covg *covgs, Edges *edges)
{
  // Binary search on the index
  long x = binary_search_index(bkey,gs->index,gs->nblocks);
  if(x < 0) return false;
  size_t blockstart = x*gs->blocksize;
  size_t blockend = (size_t)x < gs->nblocks ? blockstart+gs->blocksize : gs->nkmers;
  if(!search_file_sec(gs, bkey, blockstart, blockend)) return false;
  if(edges || covgs) graph_file_read_covgs_edges(gs->file, covgs, edges);
  return true;
}

void graph_search_fetch(GraphFileSearch *gs, size_t idx, BinaryKmer *bkey,
                             Covg *covgs, Edges *edges)
{
  const size_t entsize = sizeof(BinaryKmer)+gs->ncols*(sizeof(Covg)+sizeof(Edges));
  graph_file_fseek(gs->file, gs->file->hdr_size+entsize*idx, SEEK_SET);
  graph_file_fread(gs->file, bkey->b, sizeof(BinaryKmer));
  if(edges || covgs) graph_file_read_covgs_edges(gs->file, covgs, edges);
}

void graph_search_rand(GraphFileSearch *gs,
                            BinaryKmer *bkey, Covg *covgs, Edges *edges)
{
  size_t idx = (rand() / (double)RAND_MAX) * gs->nkmers;
  graph_search_fetch(gs, idx, bkey, covgs, edges);
}
