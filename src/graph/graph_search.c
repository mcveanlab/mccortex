#include "global.h"
#include "graph_search.h"

struct GraphFileSearch {
  GraphFileReader *file;
  size_t nkmers, ncols, entrysize; // nkmers in file, size of kmer entry in file
  BinaryKmer *index;
  size_t blocksize, nblocks;
  void *block; // read file into block to linear search
};

// #define INDEX_SIZE 4 /* debugging */
#define INDEX_SIZE 4*1024*1024 /* 4M */
#define MAX_LIN_SEARCH 512

/* with MAX_LIN_SEARCH of 512, 1 MiB allows 227 colours to be loaded */

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
  gs->ncols = file->hdr.num_of_cols;
  gs->entrysize = sizeof(BinaryKmer) + gs->ncols * (sizeof(Covg)+sizeof(Edges));
  gs->nblocks = MIN2(gs->nkmers, INDEX_SIZE);
  gs->blocksize = gs->nkmers / gs->nblocks;
  gs->nblocks = (gs->nkmers+gs->blocksize-1) / gs->blocksize;
  gs->index = ctx_calloc(sizeof(BinaryKmer), gs->nblocks+1); // sentinel
  gs->block = ctx_calloc(MAX_LIN_SEARCH * gs->entrysize, 1);
  memset(gs->index[gs->nblocks].b,0xff,BKMER_BYTES); // sentinel kmer
  status("[graph_search] on-disk-graph %zu cols %zu blocks %zu bsize %zu kmers"
         " building...", gs->ncols, gs->nblocks, gs->blocksize, gs->nkmers);
  graph_file_set_buffered(file, 0); // Turn OFF buffered input
  for(i = 0; i < gs->nblocks; i++) {
    graph_file_fseek(file, graph_file_offset(file, i*gs->blocksize), SEEK_SET);
    if(graph_file_fread(file, &gs->index[i], sizeof(BinaryKmer)) != sizeof(BinaryKmer))
      die("Cannot index graph: %s", file_filter_path(&gs->file->fltr));
  }
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
  ctx_free(gs->block);
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

// Return pointer to block of Covgs+Edges
static inline void* search_file_sec(GraphFileSearch *gs, BinaryKmer bkey,
                                    size_t start, size_t end)
{
  const size_t hdrsize = gs->file->hdr_size;
  size_t mid;
  BinaryKmer bmid;
  // Binary search
  while(start + MAX_LIN_SEARCH < end) {
    mid = (start+end) / 2;
    graph_file_fseek(gs->file, hdrsize+gs->entrysize*mid, SEEK_SET);
    if(graph_file_fread(gs->file, gs->block, gs->entrysize) != gs->entrysize)
      die("Cannot search graph from disk: %s", file_filter_path(&gs->file->fltr));
    memcpy(bmid.b, gs->block, sizeof(BinaryKmer)); // copy binary kmer
    if(binary_kmer_eq(bkey,bmid)) return gs->block;
    if(binary_kmer_lt(bkey,bmid)) end = mid;
    else start = mid + 1;
  }
  // Linear search
  size_t blockmem = gs->entrysize*(end-start);
  graph_file_fseek(gs->file, hdrsize+gs->entrysize*start, SEEK_SET);
  if(graph_file_fread(gs->file, gs->block, blockmem) != blockmem)
    die("Cannot search graph from disk: %s", file_filter_path(&gs->file->fltr));
  char *p, *endp = (char*)gs->block + blockmem;
  for(p = gs->block; p < endp; p += gs->entrysize)
  {
    memcpy(bmid.b, p, sizeof(BinaryKmer));
    if(binary_kmer_eq(bkey,bmid)) return p;
    if(binary_kmer_lt(bkey,bmid)) return NULL;
  }
  return NULL;
}

// Given an entry from a graph file, load edges and coverage
static inline void filter_covgs_edges(const FileFilter *fltr,
                                      Covg *covgs, Edges *edges,
                                      const void *ptr)
{
  size_t from, into, i;
  const char *allcovgs = (const char*)ptr + sizeof(BinaryKmer);
  const char *alledges = (const char*)allcovgs + fltr->srcncols*sizeof(Covg);
  Covg c;
  Edges e;
  memset(covgs, 0, file_filter_into_ncols(fltr) * sizeof(Covg));
  for(i = 0; i < file_filter_num(fltr); i++) {
    from = file_filter_fromcol(fltr, i);
    into = file_filter_intocol(fltr, i);
    memcpy(&c, allcovgs+sizeof(Covg)*from, sizeof(Covg));
    covgs[into] = SAFE_ADD_COVG(covgs[into], c);
  }
  memset(edges, 0, file_filter_into_ncols(fltr) * sizeof(Edges));
  for(i = 0; i < file_filter_num(fltr); i++) {
    from = file_filter_fromcol(fltr,i);
    into = file_filter_intocol(fltr, i);
    memcpy(&e, alledges+sizeof(Edges)*from, sizeof(Edges));
    edges[into] |= e;
  }
}

bool graph_search_find(GraphFileSearch *gs, BinaryKmer bkey,
                       Covg *covgs, Edges *edges)
{
  char *ptr;
  // Binary search on the index
  long x = binary_search_index(bkey,gs->index,gs->nblocks);
  if(x < 0) return false;
  size_t blockstart = x*gs->blocksize;
  size_t blockend = (size_t)x+1 < gs->nblocks ? blockstart+gs->blocksize : gs->nkmers;
  if((ptr = search_file_sec(gs, bkey, blockstart, blockend)) == NULL) return false;
  filter_covgs_edges(&gs->file->fltr, covgs, edges, ptr);
  return true;
}

void graph_search_fetch(GraphFileSearch *gs, size_t idx, BinaryKmer *bkey,
                        Covg *covgs, Edges *edges)
{
  graph_file_fseek(gs->file, gs->file->hdr_size+gs->entrysize*idx, SEEK_SET);
  // read one entry
  if(graph_file_fread(gs->file, gs->block, gs->entrysize) != gs->entrysize)
    die("Cannot search graph from disk: %s", file_filter_path(&gs->file->fltr));
  memcpy(bkey, gs->block, sizeof(BinaryKmer)); // copy binary kmer
  filter_covgs_edges(&gs->file->fltr, covgs, edges, gs->block);
}

void graph_search_rand(GraphFileSearch *gs,
                       BinaryKmer *bkey, Covg *covgs, Edges *edges)
{
  size_t idx = (rand() / (double)RAND_MAX) * gs->nkmers;
  graph_search_fetch(gs, idx, bkey, covgs, edges);
}
