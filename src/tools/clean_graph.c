#include "global.h"
#include "util.h"
#include "supernode.h"
#include "prune_nodes.h"
#include "clean_graph.h"

#define DUMP_COVG_ARRSIZE 1000

// Define a vector of Covg
#include "objbuf_macro.h"
create_objbuf(covg_buf,CovgBuffer,Covg);

#define supernode_covg(covgs,len) supernode_covg_mean(covgs,len)
// #define supernode_covg(covgs,len) supernode_read_starts(covgs,len)

static inline size_t supernode_covg_mean(Covg *covgs, size_t len)
{
  ctx_assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len; // round to nearest integer
}

// Mark each node in the supernode as visited
static inline size_t fetch_supernode(hkey_t node,
                                     dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                     uint64_t *visited, const dBGraph *db_graph)
{
  size_t i;
  hkey_t hkey;

  db_node_buf_reset(nbuf);
  covg_buf_reset(cbuf);
  supernode_find(node, nbuf, db_graph);
  covg_buf_ensure_capacity(cbuf, nbuf->len);
  cbuf->len = nbuf->len;

  for(i = 0; i < nbuf->len; i++) {
    hkey = nbuf->data[i].key;
    bitset_set(visited, hkey);
    cbuf->data[i] = db_graph->col_covgs[hkey];
  }

  return supernode_covg(cbuf->data, cbuf->len);
}

static inline void supernode_clean(hkey_t node,
                                   dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                   uint64_t *visited, Covg covg_threshold,
                                   dBGraph *db_graph)
{
  size_t reads_arriving;

  if(!bitset_get(visited, node))
  {
    reads_arriving = fetch_supernode(node, nbuf, cbuf, visited, db_graph);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_thresh: ", stdout);
      db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
      printf(" len: %zu covg: %zu threshold: %u\n", nbuf->len,
             reads_arriving, covg_threshold);
    #endif

    if(reads_arriving < covg_threshold) {
      prune_supernode(nbuf->data, nbuf->len, db_graph);
    }
  }
}

static inline void clip_tip(hkey_t node, dBNodeBuffer *nbuf,
                            uint64_t *visited, size_t min_keep_len,
                            dBGraph *db_graph)
{
  size_t i;
  Edges first, last;
  int in, out;

  if(!bitset_get(visited, node))
  {
    db_node_buf_reset(nbuf);
    supernode_find(node, nbuf, db_graph);

    for(i = 0; i < nbuf->len; i++) bitset_set(visited, nbuf->data[i].key);

    first = db_node_get_edges_union(db_graph, nbuf->data[0].key);
    last = db_node_get_edges_union(db_graph, nbuf->data[nbuf->len-1].key);
    in = edges_get_indegree(first, nbuf->data[0].orient);
    out = edges_get_outdegree(last, nbuf->data[nbuf->len-1].orient);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_clip: ", stdout);
      db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
      fprintf(stdout, " len: %zu junc: %i\n", nbuf->len, in+out);
    #endif

    if(in+out <= 1 && nbuf->len < min_keep_len)
      prune_supernode(nbuf->data, nbuf->len, db_graph);
  }
}

// min_tip_len is the min tip length to KEEP
void cleaning_remove_tips(size_t min_tip_len, uint64_t *visited,
                          dBGraph *db_graph)
{
  ctx_assert(min_tip_len > 1);

  // Need to use _SAFE hash traverse since we remove elements in clip_tip()
  status("[cleaning] Clipping tips shorter than %zu...\n", min_tip_len);

  char remain_nkmers_str[100], removed_nkmers_str[100];
  size_t init_nkmers = db_graph->ht.num_kmers, remain_nkmers, removed_nkmers;

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 2048);

  HASH_ITERATE_SAFE(&db_graph->ht, clip_tip,
                    &nbuf, visited, min_tip_len, db_graph);

  db_node_buf_dealloc(&nbuf);

  remain_nkmers = db_graph->ht.num_kmers;
  removed_nkmers = init_nkmers - remain_nkmers;
  ulong_to_str(remain_nkmers, remain_nkmers_str);
  ulong_to_str(removed_nkmers, removed_nkmers_str);
  status("[cleaning] Remaining kmers: %s removed: %s (%.1f%%)",
         remain_nkmers_str, removed_nkmers_str,
         (100.0*removed_nkmers)/init_nkmers);
}

#define status_return(msg,thresh) do { status(msg); return (thresh); } while(0)

size_t cleaning_supernode_threshold(uint64_t *covgs, size_t len,
                                    double seq_depth, const dBGraph *db_graph)
{
  ctx_assert(len > 5);
  ctx_assert(db_graph->ht.num_kmers > 0);
  size_t i, d1len = len-2, d2len = len-3, f1, f2;
  double *tmp = malloc2((d1len+d2len) * sizeof(double));
  double *delta1 = tmp, *delta2 = tmp + d1len;

  // Get sequencing depth from coverage
  uint64_t covg_sum = 0, capacity = db_graph->ht.capacity * db_graph->num_of_cols;
  for(i = 0; i < capacity; i++) covg_sum += db_graph->col_covgs[i];
  double seq_depth_est = (double)covg_sum / db_graph->ht.num_kmers;

  status("[cleaning] Kmer depth before cleaning supernodes: %.2f", seq_depth_est);
  if(seq_depth <= 0) seq_depth = seq_depth_est;
  else status("[cleaning] Using sequence depth argument: %f", seq_depth);

  size_t fallback_thresh = (size_t)MAX2(1, (seq_depth+1)/2);

  // +1 to ensure covgs is never 0
  for(i = 0; i < d1len; i++) delta1[i] = (double)(covgs[i+1]+1) / (covgs[i+2]+1);

  d1len = i;
  d2len = d1len - 1;

  if(d1len <= 2) {
    status("(using fallback1)\n");
    free(tmp);
    return fallback_thresh;
  }

  // d2len is d1len-1
  for(i = 0; i < d2len; i++) delta2[i] = delta1[i] / delta1[i+1];

  for(f1 = 0; f1 < d1len && delta1[f1] >= 1; f1++);
  for(f2 = 0; f2 < d2len && delta2[f2] > 1; f2++);

  free(tmp);

  if(f1 < d1len && f1 < (seq_depth*0.75))
    status_return("[cleaning] (using f1)", f1+1);
  else if(f2 < d2len)
    status_return("[cleaning] (using f2)", f2+1);
  else
    status_return("[cleaning] (using fallback1)", fallback_thresh+1);
}

static inline void covg_histogram(hkey_t node,
                                  dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                  uint64_t *visited, uint64_t *covg_hist,
                                  size_t covg_arrlen, const dBGraph *db_graph)
{
  size_t reads_arriving;

  if(!bitset_get(visited, node))
  {
    reads_arriving = fetch_supernode(node, nbuf, cbuf, visited, db_graph);
    reads_arriving = MIN2(reads_arriving, covg_arrlen-1);
    covg_hist[reads_arriving]++;
  }
}


// If covg_threshold is zero, uses covg distribution to calculate
// Returns covg threshold used
// visited will be dirty after calling
// Returns 0 if `covg_threshold` == 0 and computed threshold also <= 1
//  => no supernodes can be removed
Covg cleaning_remove_supernodes(bool do_cleaning, Covg covg_threshold,
                                double seq_depth, const char *dump_covgs,
                                uint64_t *visited, dBGraph *db_graph)
{
  if(db_graph->ht.num_kmers == 0) return covg_threshold;

  dBNodeBuffer nbuf;
  CovgBuffer cbuf;
  db_node_buf_alloc(&nbuf, 2048);
  covg_buf_alloc(&cbuf, 2048);

  size_t visited_words = roundup_bits2words64(db_graph->ht.capacity);

  char remain_nkmers_str[100], removed_nkmers_str[100];
  size_t init_nkmers = db_graph->ht.num_kmers, remain_nkmers, removed_nkmers;

  // Estimate optimum cleaning threshold
  if(covg_threshold == 0 || dump_covgs != NULL)
  {
    status("[cleaning] Calculating supernode cleaning threshold...");

    // Get supernode coverages
    uint64_t *covg_hist = calloc2(DUMP_COVG_ARRSIZE, sizeof(uint64_t));

    HASH_ITERATE(&db_graph->ht, covg_histogram,
                 &nbuf, &cbuf, visited, covg_hist, DUMP_COVG_ARRSIZE, db_graph);

    if(dump_covgs != NULL)
      cleaning_dump_covg_histogram(dump_covgs, covg_hist, DUMP_COVG_ARRSIZE);

    if(do_cleaning)
      memset(visited, 0, visited_words * sizeof(uint64_t));

    // set threshold using histogram and genome size
    size_t threshold_est = cleaning_supernode_threshold(covg_hist,
                                                        DUMP_COVG_ARRSIZE,
                                                        seq_depth, db_graph);

    status("[cleaning] Recommended supernode cleaning threshold: < %zu",
           threshold_est);

    if(covg_threshold == 0)
      covg_threshold = (Covg)threshold_est;

    free(covg_hist);
  }

  if(do_cleaning)
  {
    // Remove low coverage supernodes
    status("[cleaning] Cleaning supernodes with coverage < %zu...",
           (size_t)covg_threshold);

    if(covg_threshold <= 1)
      warn("Supernode cleaning failed, cleaning with threshold of <= 1");
    else {
      HASH_ITERATE(&db_graph->ht, supernode_clean,
                   &nbuf, &cbuf, visited,
                   covg_threshold, db_graph);
    }

    remain_nkmers = db_graph->ht.num_kmers;
    removed_nkmers = init_nkmers - remain_nkmers;
    ulong_to_str(remain_nkmers, remain_nkmers_str);
    ulong_to_str(removed_nkmers, removed_nkmers_str);
    status("[cleaning] Remaining kmers: %s removed: %s (%.1f%%)",
           remain_nkmers_str, removed_nkmers_str,
           (100.0*removed_nkmers)/init_nkmers);
  }

  db_node_buf_dealloc(&nbuf);
  covg_buf_dealloc(&cbuf);

  return covg_threshold <= 1 ? 0 : covg_threshold;
}

void cleaning_dump_covg_histogram(const char *path, uint64_t *hist, size_t len)
{
  if(len == 0) return;
  ctx_assert(len >= 2);
  ctx_assert(hist[0] == 0);

  status("[cleaning] Writing covg distribution to: %s", path);
  size_t i, end;
  FILE *fout = fopen(path, "w");

  if(fout == NULL) {
    warn("Couldn't write size distribution! file: %s", path);
    return;
  }

  fprintf(fout, "Covg,Supernodes\n");
  fprintf(fout, "1,%zu\n", (size_t)hist[1]);
  for(end = len-1; end > 1 && hist[end] == 0; end--);
  for(i = 2; i <= end; i++) {
    if(hist[i] > 0)
      fprintf(fout, "%zu,%zu\n", i, (size_t)hist[i]);
  }
  fclose(fout);
}
