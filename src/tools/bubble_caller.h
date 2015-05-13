#ifndef BUBBLE_CALLER_H_
#define BUBBLE_CALLER_H_

#include "db_graph.h"
#include "graph_cache.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "cmd.h"

#include "cJSON/cJSON.h"

#define BUBBLE_FORMAT_VERSION 2

typedef struct
{
  // Max lengths in kmers (not bases)
  const size_t max_allele_len, max_flank_len;
  // haploid colours are samples which are haploid (e.g. reference genomes)
  // these are used to filter out repeats (e.g. if haploid covg on both alleles)
  const size_t *haploid_cols, num_haploid;
} BubbleCallingPrefs;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(cache_stepptr_buf, GCacheStepPtrBuf, GCacheStep*);

typedef struct
{
  // Specific to this instance
  const size_t threadid, nthreads;

  // Temporary memory specific to this instance
  GraphCache cache;
  bool *const haploid_seen; // used to record which of the haploids we've seen
  GCacheStepPtrBuf spp_forward, spp_reverse;

  dBNodeBuffer flank5p, pathbuf;
  GraphWalker wlk;
  RepeatWalker rptwlk;

  StrBuf output_buf;

  // Shared data
  size_t *num_bubbles_ptr; // statistics - shared pointer
  const BubbleCallingPrefs prefs;
  const dBGraph *db_graph;
  gzFile gzout;
  pthread_mutex_t *const out_lock;
} BubbleCaller;

BubbleCaller* bubble_callers_new(size_t num_callers,
                                 BubbleCallingPrefs prefs,
                                 gzFile gzout,
                                 const dBGraph *db_graph);

void bubble_callers_destroy(BubbleCaller *callers, size_t num_callers);

// `fork_node` is a node with outdegree > 1
void find_bubbles(BubbleCaller *caller, dBNode fork_node);

// Load GCacheSteps into caller->spp_forward (if they traverse the snode forward)
// or caller->spp_reverse (if they traverse the snode in reverse)
void find_bubbles_ending_with(BubbleCaller *caller, GCacheSnode *snode);

// Run bubble caller, write output to gzout
// @param hdrs JSON headers of input files
// @param nhdrs number of JSON headers of input files
void invoke_bubble_caller(size_t num_of_threads, BubbleCallingPrefs prefs,
                          gzFile gzout, const char *out_path,
                          cJSON **hdrs, size_t nhdrs,
                          const dBGraph *db_graph);

#endif /* BUBBLE_CALLER_H_ */
