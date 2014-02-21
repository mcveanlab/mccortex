#include "global.h"
#include "build_graph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "async_read_io.h"
#include "loading_stats.h"

#include <pthread.h>
#include "seq_file.h"

typedef struct
{
  pthread_t thread;
  dBGraph *const db_graph;
  MsgPool *pool;
  LoadingStats *file_stats; // Array of stats for diff input files
} BuildGraphWorker;

//
// Check for PCR duplicates
//

// Read start (duplicate removal during read loading)
#define db_node_has_read_start_mt(graph,hkey,or) \
        bitset_get_mt((volatile uint8_t*)(graph)->readstrt, 2*(hkey)+(or))
#define db_node_set_read_start_mt(graph,hkey,or) \
        bitset_set_mt((volatile uint8_t*)(graph)->readstrt, 2*(hkey)+(or))

boolean seq_reads_are_novel(const read_t *r1, const read_t *r2,
                            uint8_t qual_cutoff1, uint8_t qual_cutoff2,
                            uint8_t hp_cutoff,
                            dBGraph *db_graph)
{
  hkey_t node1 = HASH_NOT_FOUND, node2 = HASH_NOT_FOUND;
  Orientation or1 = FORWARD, or2 = FORWARD;
  BinaryKmer curr_kmer, tmp_key;
  boolean found1 = false, found2 = false;
  size_t start1, start2, kmer_size = db_graph->kmer_size;

  start1 = seq_contig_start(r1, 0, kmer_size, qual_cutoff1, hp_cutoff);
  start2 = seq_contig_start(r2, 0, kmer_size, qual_cutoff2, hp_cutoff);

  boolean got_kmer1 = (start1 < r1->seq.end);
  boolean got_kmer2 = (start2 < r2->seq.end);

  if(!got_kmer1 && !got_kmer2) return 1;

  // Look up first kmer
  if(got_kmer1)
  {
    curr_kmer = binary_kmer_from_str(r1->seq.b+start1, kmer_size);
    tmp_key = bkmer_get_key(curr_kmer, kmer_size);
    node1 = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found1);
    or1 = bkmer_get_orientation(curr_kmer, tmp_key);
  }

  // Look up second kmer
  if(got_kmer2)
  {
    curr_kmer = binary_kmer_from_str(r2->seq.b+start2, kmer_size);
    tmp_key = bkmer_get_key(curr_kmer, kmer_size);
    node2 = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found2);
    or2 = bkmer_get_orientation(curr_kmer, tmp_key);
  }

  // Each read gives no kmer or a duplicate kmer
  // used find_or_insert so if we have a kmer we have a graph node
  if((!got_kmer1 || db_node_has_read_start_mt(db_graph, node1, or1)) &&
     (!got_kmer2 || db_node_has_read_start_mt(db_graph, node2, or2)))
  {
    return false;
  }

  // Read is novel
  if(got_kmer1) db_node_set_read_start_mt(db_graph, node1, or1);
  if(got_kmer2) db_node_set_read_start_mt(db_graph, node2, or2);
  return true;
}

boolean seq_read_is_novel(const read_t *r, dBGraph *db_graph,
                          uint8_t qual_cutoff, uint8_t hp_cutoff)
{
  hkey_t node = HASH_NOT_FOUND;
  Orientation or;
  BinaryKmer bkmer, tmp_key;
  boolean found = false;

  size_t start = seq_contig_start(r, 0, db_graph->kmer_size,
                                  qual_cutoff, hp_cutoff);

  if(start == r->seq.end) return 1;

  bkmer = binary_kmer_from_str(r->seq.b+start, db_graph->kmer_size);
  tmp_key = bkmer_get_key(bkmer, db_graph->kmer_size);
  node = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found);
  or = bkmer_get_orientation(bkmer, tmp_key);

  if(db_node_has_read_start_mt(db_graph, node, or)) return false;

  // Read is novel
  db_node_set_read_start_mt(db_graph, node, or);
  return true;
}


//
// Add to the de bruijn graph
//

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
void build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                             const char *seq, size_t len)
{
  ctx_assert(len >= db_graph->kmer_size);
  const size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  size_t i;

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = db_graph_find_or_add_node_mt(db_graph, bkmer, colour);

  for(i = kmer_size; i < len; i++)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = db_graph_find_or_add_node_mt(db_graph, bkmer, colour);
    db_graph_add_edge_mt(db_graph, colour, prev, curr);
    prev = curr;
  }
}

static void load_read(const read_t *r, dBGraph *db_graph,
                      uint8_t qual_cutoff, uint8_t hp_cutoff,
                      Colour colour, LoadingStats *stats)
{
  const size_t kmer_size = db_graph->kmer_size;
  if(r->seq.end < kmer_size) { stats->total_bad_reads++; return; }

  size_t search_start = 0;
  size_t contig_start, contig_end = 0;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    size_t contig_len = contig_end - contig_start;

    build_graph_from_str_mt(db_graph, colour, r->seq.b+contig_start, contig_len);

    size_t kmers_loaded = contig_len + 1 - kmer_size;
    stats->total_bases_loaded += contig_len;
    stats->kmers_loaded += kmers_loaded;
    stats->contigs_loaded++;
  }

  // contig_end == 0 if no contigs from this read
  if(contig_end == 0) stats->total_bad_reads++;
  else stats->total_good_reads++;
}

static boolean check_pcr_dupe(read_t *r1, read_t *r2,
                              uint8_t fq_cutoff1, uint8_t fq_cutoff2,
                              uint8_t hp_cutoff,
                              boolean remove_dups_se, boolean remove_dups_pe,
                              ReadMateDir matedir, dBGraph *db_graph)
{
  if(!remove_dups_se && !remove_dups_pe) return false;

  boolean samdupe1, samdupe2;
  samdupe1 = (r1->from_sam && r1->bam->core.flag & BAM_FDUP);
  samdupe2 = (r2 == NULL || (r2->from_sam && r2->bam->core.flag & BAM_FDUP));

  if(remove_dups_pe && samdupe1 && samdupe2) return true;
  if(remove_dups_se && samdupe1) return true;

  seq_reader_orient_mp_FF(r1, r2, matedir);

  // Use first kmer from each read to check if duplicate
  if(remove_dups_pe && r2 != NULL &&
     !seq_reads_are_novel(r1, r2, fq_cutoff1, fq_cutoff2, hp_cutoff, db_graph))
    return true;

  if(remove_dups_se && r2 == NULL &&
     !seq_read_is_novel(r1, db_graph, fq_cutoff1, hp_cutoff))
    return true;

  return false;
}

void build_graph_from_reads(read_t *r1, read_t *r2,
                            uint8_t fq_offset1, uint8_t fq_offset2,
                            uint8_t fq_cutoff, uint8_t hp_cutoff,
                            boolean remove_dups_se, boolean remove_dups_pe,
                            ReadMateDir matedir,
                            LoadingStats *stats, size_t colour,
                            dBGraph *db_graph)
{
  // status("r1: '%s' '%s'", r1->name.b, r1->seq.b);
  // if(r2) status("r2: '%s' '%s'", r2->name.b, r2->seq.b);

  stats->total_bases_read += r1->seq.end + (r2 ? r2->seq.end : 0);

  uint8_t fq_cutoff1, fq_cutoff2;
  fq_cutoff1 = fq_cutoff2 = fq_cutoff;

  if(fq_cutoff > 0)
  {
    fq_cutoff1 += fq_offset1;
    fq_cutoff2 += fq_offset2;
  }

  if(check_pcr_dupe(r1, r2, fq_cutoff1, fq_cutoff2, hp_cutoff,
                    remove_dups_se, remove_dups_pe, matedir, db_graph))
  {
    stats->total_dup_reads += (r2 == NULL ? 1 : 2);
    return;
  }

  load_read(r1, db_graph, fq_cutoff1, hp_cutoff, colour, stats);
  if(r2 != NULL) load_read(r2, db_graph, fq_cutoff2, hp_cutoff, colour, stats);
}

// pthread method, loop: reads from pool, add to graph
static void* grab_reads_from_pool(void *ptr)
{
  BuildGraphWorker *wrkr = (BuildGraphWorker*)ptr;
  MsgPool *pool = wrkr->pool;
  BuildGraphTask *task;
  AsyncIOData *data;
  int pos;

  while((pos = msgpool_claim_read(pool)) != -1)
  {
    memcpy(&data, msgpool_get_ptr(pool, pos), sizeof(AsyncIOData*));
    task = (BuildGraphTask*)data->ptr;

    build_graph_from_reads(&data->r1, &data->r2,
                           data->fq_offset1, data->fq_offset2,
                           task->fq_cutoff, task->hp_cutoff,
                           task->remove_dups_se, task->remove_dups_pe,
                           task->matedir,
                           &wrkr->file_stats[task->idx],
                           task->colour, wrkr->db_graph);

    msgpool_release(pool, pos, MPOOL_EMPTY);
  }

  pthread_exit(NULL);
}

// One thread used per input file, num_build_threads used to add reads to graph
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads)
{
  ctx_assert(db_graph->bktlocks != NULL);

  size_t i, f;

  AsyncIOData *data = malloc2(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  MsgPool pool;
  msgpool_alloc_yield(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*));
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  // Start async io reading
  AsyncIOReadTask *async_tasks = malloc2(num_files * sizeof(AsyncIOReadTask));

  for(f = 0; f < num_files; f++) {
    files[f].idx = f;
    AsyncIOReadTask aio_task = {.file1 = files[f].file1,
                                .file2 = files[f].file2,
                                .fq_offset = files[f].fq_offset,
                                .ptr = &files[f]};

    memcpy(&async_tasks[f], &aio_task, sizeof(AsyncIOReadTask));
  }

  BuildGraphWorker *workers = malloc2(num_build_threads * sizeof(BuildGraphWorker));

  for(i = 0; i < num_build_threads; i++)
  {
    BuildGraphWorker tmp_wrkr = {.db_graph = db_graph, .pool = &pool};
    tmp_wrkr.file_stats = calloc2(num_files, sizeof(LoadingStats));
    memcpy(&workers[i], &tmp_wrkr, sizeof(BuildGraphWorker));
  }

  // Create a lot of workers to build the graph
  asyncio_run_threads(&pool, async_tasks, num_files, grab_reads_from_pool,
                      workers, num_build_threads, sizeof(BuildGraphWorker));

  // start_build_graph_workers(&pool, db_graph, files, num_files, num_build_threads);

  free(async_tasks);
  msgpool_dealloc(&pool);

  // Clean up workers one by one...
  for(i = 0; i < num_build_threads; i++) {
    // Merge stats
    for(f = 0; f < num_files; f++)
      loading_stats_merge(&files[f].stats, &workers[i].file_stats[f]);

    // Free memory
    free(workers[i].file_stats);
  }

  free(workers);

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  free(data);
}
