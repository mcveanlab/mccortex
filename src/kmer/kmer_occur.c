#include "global.h"
#include "kmer_occur.h"
#include "seq_reader.h"
#include "util.h"

//
// This file provides a datastore for loading sequences and recording where
// each kmer occurs in the sequences. Used in breakpoint_caller.c.
//

// Malloc's and returns chrom list - remember to free
// Returns NULL if num_reads == 0
static void generate_chrom_list(KOGraph *kograph,
                                const read_t *reads, size_t num_reads)
{
  kograph->chrom_name_buf = NULL;
  kograph->nchroms = num_reads;
  kograph->chroms = NULL;

  if(num_reads == 0) return;

  kograph->chroms = ctx_malloc(num_reads*sizeof(KOChrom));
  size_t i, name_len_sum = 0;

  // Concatenate names into one string
  for(i = 0; i < num_reads; i++) name_len_sum += reads[i].name.end+1;
  kograph->chrom_name_buf = ctx_malloc(name_len_sum);
  char *names_concat = kograph->chrom_name_buf;

  for(i = 0; i < num_reads; i++)
  {
    kograph->chroms[i] = (KOChrom){.id = i,
                                   .length = reads[i].seq.end,
                                   .name = names_concat};
  
    memcpy(names_concat, reads[i].name.b, reads[i].name.end);
    reads[i].name.b[reads[i].name.end] = '\0';
    names_concat += reads[i].name.end+1;
  }
}

// static void bkmer_update_counts_add(BinaryKmer bkmer, KONodeList *klists,
//                                     dBGraph *db_graph)
// {
//   bool found;
//   hkey_t hkey = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
//   klists[hkey].count++;
// }

// Add missing kmers and edges to the graph whilst keeping track of the count
// of how many times each kmer in the graph is seen in sequence (with klists)
// Threadsafe
static size_t add_read_to_graph_mt(const char *seq, size_t len,
                                   KONodeList *klists, dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  if(len < kmer_size) return 0;

  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  size_t i, num_novel_kmers = 0;
  bool found;

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
  __sync_fetch_and_add((volatile uint32_t*)&klists[prev.key].count, 1); // count++
  num_novel_kmers += !found;

  for(i = kmer_size; i < len; i++, prev = curr)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
    __sync_fetch_and_add((volatile uint32_t*)&klists[curr.key].count, 1); // count++
    db_graph_add_edge_mt(db_graph, 0, prev, curr);
    num_novel_kmers += !found;
  }

  return num_novel_kmers;
}

// Same as above but don't add missing kmers
// Threadsafe
static inline void bkmer_update_counts_find_mt(BinaryKmer bkmer,
                                               KONodeList *klists,
                                               const dBGraph *db_graph)
{
  hkey_t hkey = hash_table_find(&db_graph->ht, bkmer);
  if(hkey != HASH_NOT_FOUND)
    __sync_fetch_and_add((volatile uint32_t*)&klists[hkey].count, 1); // count++
}

struct ReadUpdateCounts {
  const read_t *r;
  KONodeList *klists;
  bool add_missing_kmers;
  dBGraph *db_graph;
};

// Multithreaded core function to store kmer occurances
// Also adds read to the graph if `add_missing_kmers` is true
static void read_update_counts(void *arg)
{
  struct ReadUpdateCounts data = *(struct ReadUpdateCounts*)arg;
  const read_t *r = data.r;

  if(data.add_missing_kmers) {
    add_read_to_graph_mt(r->seq.b, r->seq.end, data.klists, data.db_graph);
  }
  else {
    LoadingStats stats = LOAD_STATS_INIT_MACRO;
    READ_TO_BKMERS(r, data.db_graph->kmer_size, 0, 0, &stats,
                   bkmer_update_counts_find_mt, data.klists, data.db_graph);
  }
}

static void bkmer_store_kmer_pos(BinaryKmer bkmer, const dBGraph *db_graph,
                                 KONodeList *klists, KOccur *koccurs,
                                 size_t chrom_id, uint64_t offset)
{
  // bkmers were already added to graph -> don't need to find_or_insert
  // if missing kmers weren't added then kmer might be missing -> skip
  dBNode node = db_graph_find(db_graph, bkmer);

  if(node.key != HASH_NOT_FOUND)
  {
    KONodeList *kl = &klists[node.key];
    koccurs[kl->start+kl->count] = (KOccur){.chrom = chrom_id,
                                            .offset = offset,
                                            .orient = node.orient};
    kl->count++;
  }
}

static void read_store_kmer_pos(const read_t *r, size_t chrom_id,
                                KONodeList *klists, KOccur *koccurs,
                                const dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  LoadingStats stats = LOAD_STATS_INIT_MACRO;
  READ_TO_BKMERS(r, kmer_size, 0, 0, &stats, bkmer_store_kmer_pos,
                 db_graph, klists, koccurs, chrom_id, _offset);
}

// BEWARE: We add the reads to the graph
KOGraph kograph_create(const read_t *reads, size_t num_reads,
                       bool add_missing_kmers, size_t num_threads,
                       dBGraph *db_graph)
{
  if(add_missing_kmers) ctx_assert(db_graph->num_edge_cols == 1);
  ctx_assert(sizeof(KONodeList) == 12);

  KOGraph kograph;
  generate_chrom_list(&kograph, reads, num_reads);

  kograph.klists = ctx_calloc(db_graph->ht.capacity, sizeof(KONodeList));

  // 1. Loop through reads, add to graph and record kmer counts
  size_t i;
  struct ReadUpdateCounts *updates = ctx_malloc(num_reads * sizeof(struct ReadUpdateCounts));
  for(i = 0; i < num_reads; i++) {
    updates[i] = (struct ReadUpdateCounts){.r = &reads[i],
                                           .klists = kograph.klists,
                                           .add_missing_kmers = add_missing_kmers,
                                           .db_graph = db_graph};
  }

  util_run_threads(updates, num_reads, sizeof(struct ReadUpdateCounts),
                   num_threads, read_update_counts);

  ctx_free(updates);

  // 2. alloc lists
  uint64_t offset = 0, total_read_length = 0;

  for(i = 0; i < num_reads; i++)
  {
    kograph.klists[i].start = offset;
    offset += kograph.klists[i].count;
    kograph.klists[i].count = 0;
    total_read_length += reads[i].seq.end;
  }

  // offset should be less than sum of read lengths, which is why we wait to
  // parse all reads before mallocing to reduce memory footprint 
  ctx_assert(offset < total_read_length);

  if(offset == 0) kograph.koccurs = NULL;
  else {
    kograph.koccurs = ctx_malloc(offset * sizeof(KOccur));

    // 3. Loop through reads, record kmer pos
    for(i = 0; i < num_reads; i++) {
      read_store_kmer_pos(&reads[i], i, kograph.klists, kograph.koccurs, db_graph);
    }
  }

  return kograph;
}

void kograph_free(KOGraph kograph)
{
  if(kograph.chrom_name_buf != NULL) ctx_free(kograph.chrom_name_buf);
  if(kograph.chroms != NULL) ctx_free(kograph.chroms);
  if(kograph.klists != NULL) ctx_free(kograph.klists);
  if(kograph.koccurs != NULL) ctx_free(kograph.koccurs);
}

// Sort by chromosome then by position
static inline int kmer_occur_cmp(KOccur occur0, KOccur occur1, size_t dist)
{
  size_t offset0 = occur0.offset + dist, offset1 = occur1.offset;
  if(occur0.chrom != occur1.chrom) return (int)occur0.chrom - occur1.chrom;
  if(offset0 != offset1) return (offset0 < offset1 ? -1 : 1);
  return 0;
}

// occur0 should be before occur1
// kolist0 and kolist1 should already be sorted
void kograph_get_kmer_run(const KOccur *restrict kolist0, size_t num0,
                          const KOccur *restrict kolist1, size_t num1,
                          size_t dist, KOccurBuffer *restrict kobuf)
{
  size_t i, j;
  int cmp;
  for(i = j = 0; i < num0 && j < num1; ) {
    cmp = kmer_occur_cmp(kolist0[i], kolist1[j], dist);
    if(cmp < 0) i++;
    else if(cmp > 0) j++;
    else { kmer_occur_buf_add(kobuf, kolist0[i]); i++; j++; }
  }
  for(; i < num0; i++) kmer_occur_buf_add(kobuf, kolist0[i]);
  for(; j < num1; j++) kmer_occur_buf_add(kobuf, kolist1[j]);
}

/*
// Filter regions down to only those that stretch the whole distance
void kograph_filter_stretch(KOGraph kograph, const dBNode *nodes, size_t len,
                            KOccurBuffer *kobuf, KOccurBuffer *kobuftmp)
{
  kmer_occur_buf_reset(kobuf);
  kmer_occur_buf_reset(kobuftmp);
  if(len == 0) return;

  KOccur *kolist0, *kolist1, *kolist2;
  size_t n;

  kolist0 = kograph_get_check(kograph, nodes[0].key);
  if(kolist0 == NULL) return;
  if(len == 1) {
    kmer_occur_buf_append(kobuf, kolist0, kograph_num(nodes[0].key));
    return;
  }

  for(i = 1; i < len; i++) {
    kolist1 = kograph_get(kograph, nodes[1].key);
    n = kograph_num(kograph, nodes[1].key);
    if(n == 0) { kmer_occur_buf_reset(kobuf); return; }

  }
}
*/
