#include "global.h"

#include <pthread.h>
#include <unistd.h> // usleep

#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "file_util.h"
#include "seq_reader.h"
#include "loading_stats.h"
#include "add_path_workers.h"
#include "add_read_paths.h"

// Temp data store for adding paths to graph
typedef struct
{
  const size_t wid; // Worker id
  dBGraph *const db_graph;
  GraphWalker wlk; // traversing covg gaps
  RepeatWalker rptwlk; // traversing covg gaps
  dBNodeBuffer nodebuf; // Contigs constructed here
  // Biggest gap between reads we'll try to traverse
  uint64_t *const insert_sizes, *const gap_sizes; // length gap_arr_cap
  size_t gap_arr_cap;
  volatile boolean got_job, completed;
} AddPathsWorker;

struct PathsWorkerPool
{
  // Threads
  size_t num_of_threads;
  pthread_t *threads;
  AddPathsWorker *workers;

  // Shared data
  dBGraph *db_graph;
  boolean seen_pe;

  // Data currently being read
  AddPathsJob rjob;
};

//
// Multithreading
//
pthread_attr_t thread_attr;
pthread_mutex_t data_written_mutex, data_read_mutex;
pthread_cond_t data_written_cond, data_read_cond;

// Incoming data
volatile boolean data_waiting = false, input_ended = false;
volatile AddPathsJob next_job; // job waiting

static void paths_worker_alloc(AddPathsWorker *worker, size_t wid,
                               size_t gap_arr_cap, dBGraph *db_graph)
{
  uint64_t *insert_sizes, *gap_sizes;
  insert_sizes = calloc2(2*gap_arr_cap, sizeof(uint64_t));
  gap_sizes = insert_sizes + gap_arr_cap;

  AddPathsWorker tmp = {.wid = wid, .db_graph = db_graph,
                        .insert_sizes = insert_sizes, .gap_sizes = gap_sizes,
                        .gap_arr_cap = gap_arr_cap,
                        .got_job = false, .completed = false};
  memcpy(worker, &tmp, sizeof(AddPathsWorker));

  db_node_buf_alloc(&worker->nodebuf, 4096);
  graph_walker_alloc(&worker->wlk);
  rpt_walker_alloc(&worker->rptwlk, db_graph->ht.capacity, 22); // 4MB

  // add_path_job_alloc(&worker->job);
}

static void paths_worker_dealloc(AddPathsWorker *worker)
{
  free(worker->insert_sizes);
  rpt_walker_dealloc(&worker->rptwlk);
  graph_walker_dealloc(&worker->wlk);
  db_node_buf_dealloc(&worker->nodebuf);
  // add_path_job_dealloc(&worker->job);
}

// Multithreaded notes
// pthread_cond_wait releases the lock

void* add_paths_thread(void *ptr)
{
  AddPathsWorker *worker = (AddPathsWorker*)ptr;
  dBGraph *db_graph = worker->db_graph;
  const size_t kmer_size = db_graph->kmer_size;
  worker->got_job = false;

  AddPathsJob job;
  add_path_job_alloc(&job);

  while(data_waiting || !input_ended)
  {
    pthread_mutex_lock(&data_written_mutex);

    // read in job => wait while [input not ended] and [data not ready]
    while(!data_waiting && !input_ended)
      pthread_cond_wait(&data_written_cond, &data_written_mutex);

    if(data_waiting) {
      AddPathsJob tmp_job;
      SWAP(job, next_job, tmp_job);
      worker->got_job = true;

      pthread_mutex_lock(&data_read_mutex);
      data_waiting = false;
      pthread_cond_signal(&data_read_cond);
      pthread_mutex_unlock(&data_read_mutex);
    }

    pthread_mutex_unlock(&data_written_mutex);

    // Do work
    if(worker->got_job)
    {
      if(job.r2.seq.end >= kmer_size) seq_read_reverse_complement(&job.r2);
      add_read_paths(&job, &worker->nodebuf,
                     &worker->wlk, &worker->rptwlk,
                     worker->insert_sizes, worker->gap_sizes, worker->db_graph);
      worker->got_job = false;
    }
  }

  add_path_job_dealloc(&job);

  worker->completed = true;
  pthread_exit(NULL);
}

// This function is passed to parse_filelist to load paths from sequence data
static void load_paths(read_t *r1, read_t *r2,
                       uint8_t fq_offset1, uint8_t fq_offset2,
                       const SeqReadingPrefs prefs, SeqLoadingStats *stats,
                       void *ptr)
{
  // Don't bother checking for duplicates
  (void)stats;

  PathsWorkerPool *pool = (PathsWorkerPool*)ptr;

  // printf("READ1: %s\n", r1->seq.b);
  // if(r2 != NULL) printf("READ2: %s\n", r2->seq.b);

  uint8_t qcutoff1 = prefs.quality_cutoff;
  uint8_t qcutoff2 = prefs.quality_cutoff;

  if(prefs.quality_cutoff > 0)
  {
    qcutoff1 += fq_offset1;
    qcutoff2 += fq_offset2;
  }

  // Wait until data has been read from next_job
  if(data_waiting) {
    pthread_mutex_lock(&data_read_mutex);
    while(data_waiting) pthread_cond_wait(&data_read_cond, &data_read_mutex);
    pthread_mutex_unlock(&data_read_mutex);
  }

  next_job.qcutoff1 = qcutoff1;
  next_job.qcutoff2 = qcutoff2;
  next_job.hp_cutoff = prefs.homopolymer_cutoff;
  next_job.ctp_col = pool->rjob.ctp_col;
  next_job.ctx_col = pool->rjob.ctx_col;
  next_job.ins_gap_min = pool->rjob.ins_gap_min;
  next_job.ins_gap_max = pool->rjob.ins_gap_max;

  read_t tmp_read;
  SWAP(*r1, next_job.r1, tmp_read);
  if(r2 != NULL) SWAP(*r2, next_job.r2, tmp_read);
  else {
    next_job.r2.seq.end = 0;
    next_job.r2.seq.b[0] = '\0';
  }

  // Pass to a thread
  pthread_mutex_lock(&data_written_mutex);
  data_waiting = true;
  pthread_cond_signal(&data_written_cond);
  pthread_mutex_unlock(&data_written_mutex);

  // Update stats
  stats->total_good_reads += 1 + (r2 != NULL);
}

static void _path_workers_wait_til_data_read()
{
  // Wait until all data consumed by workers
  while(data_waiting) {
    // Signal data waiting
    pthread_mutex_lock(&data_written_mutex);
    pthread_cond_signal(&data_written_cond);
    pthread_mutex_unlock(&data_written_mutex);
    // Wait for data to be taken
    pthread_mutex_lock(&data_read_mutex);
    if(data_waiting)
      pthread_cond_wait(&data_read_cond, &data_read_mutex);
    pthread_mutex_unlock(&data_read_mutex);
  }
  status("Jobs finished");
}

void path_workers_wait_til_finished(PathsWorkerPool *pool)
{
  status("Waiting for jobs to complete...");
  size_t i;
  _path_workers_wait_til_data_read();
  // Wait until all workers finish job
  while(1) {
    for(i = 0; i < pool->num_of_threads && !pool->workers[i].got_job; i++) {}
    if(i < pool->num_of_threads) usleep(500);
    else break;
  }
}

void path_workers_wait_til_completed(PathsWorkerPool *pool)
{
  status("Waiting for threads to complete...");
  size_t i;
  _path_workers_wait_til_data_read();
  // Release all threads waiting for data
  while(1) {
    for(i = 0; i < pool->num_of_threads && pool->workers[i].completed; i++) {}
    if(i < pool->num_of_threads) pthread_cond_signal(&data_written_cond);
    else break;
  }
}

PathsWorkerPool* path_workers_pool_new(size_t num_of_threads,
                                       dBGraph *db_graph, size_t gap_arr_cap)
{
  PathsWorkerPool *pool = calloc2(1, sizeof(PathsWorkerPool));

  size_t i;
  pool->workers = malloc2(sizeof(AddPathsWorker) * num_of_threads);
  pool->num_of_threads = num_of_threads;
  pool->db_graph = db_graph;
  pool->seen_pe = false;

  if(!seq_read_alloc(&pool->rjob.r1) || !seq_read_alloc(&pool->rjob.r2))
    die("Out of memory");

  for(i = 0; i < num_of_threads; i++)
    paths_worker_alloc(&pool->workers[i], i, gap_arr_cap, db_graph);

  // Start threads
  pool->threads = malloc2(num_of_threads * sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  if(pthread_mutex_init(&data_written_mutex, NULL) != 0) die("mutex init failed");
  if(pthread_mutex_init(&data_read_mutex, NULL) != 0) die("mutex init failed");

  if(pthread_cond_init(&data_written_cond, NULL) != 0) die("pthread_cond init failed");
  if(pthread_cond_init(&data_read_cond, NULL) != 0) die("pthread_cond init failed");

  if(seq_read_alloc((read_t*)&next_job.r1) == NULL ||
     seq_read_alloc((read_t*)&next_job.r2) == NULL)
    die("Out of memory");

  data_waiting = input_ended = false;

  int rc;
  for(i = 0; i < num_of_threads; i++)
  {
    rc = pthread_create(pool->threads+i, &thread_attr, add_paths_thread,
                        (void*)(pool->workers+i));
    if(rc != 0) die("Creating thread failed");
  }

  return pool;
}

// Read numbers are for printing out only
void path_workers_pool_dealloc(PathsWorkerPool *pool,
                               size_t num_se_reads, size_t num_pe_readpairs)
{
  AddPathsWorker *workers = pool->workers;
  size_t i, j, gap_arr_cap = workers[0].gap_arr_cap;

  input_ended = true;

  // Wait for workers to finish
  path_workers_wait_til_completed(pool);

  int rc;
  for(i = 0; i < pool->num_of_threads; i++) {
    rc = pthread_join(pool->threads[i], NULL);
    if(rc != 0) die("Joining thread failed");
  }

  // Sum insert sizes and gap sizes into worker 0
  uint64_t *insert_sizes = workers[0].insert_sizes;
  uint64_t *gap_sizes = workers[0].gap_sizes;

  for(i = 1; i < pool->num_of_threads; i++) {
    for(j = 0; j < gap_arr_cap; j++) {
      insert_sizes[j] += workers[i].insert_sizes[j];
      gap_sizes[j] += workers[i].gap_sizes[j];
    }
    paths_worker_dealloc(&workers[i]);
  }

  (void)num_se_reads; (void)num_pe_readpairs;
  // Print mp gap size / insert stats to a file
  // size_t kmer_size = pool->db_graph->kmer_size;
  // dump_gap_sizes("gap_sizes.%u.csv", gap_sizes, gap_arr_cap,
  //                kmer_size, false, num_se_reads + num_pe_readpairs*2);

  // if(pool->seen_pe) {
  //   dump_gap_sizes("mp_sizes.%u.csv", insert_sizes, gap_arr_cap,
  //                  kmer_size, true, num_pe_readpairs);
  // }

  paths_worker_dealloc(&workers[0]);

  free(pool->threads);
  free(pool->workers);

  seq_read_dealloc(&pool->rjob.r1);
  seq_read_dealloc(&pool->rjob.r2);

  seq_read_dealloc((read_t*)&next_job.r1);
  seq_read_dealloc((read_t*)&next_job.r2);

  pthread_cond_destroy(&data_written_cond);
  pthread_cond_destroy(&data_read_cond);

  pthread_mutex_destroy(&data_written_mutex);
  pthread_mutex_destroy(&data_read_mutex);
  pthread_attr_destroy(&thread_attr);

  free(pool);
}

void path_workers_add_paths_to_graph(PathsWorkerPool *pool,
                                     seq_file_t *sf1, seq_file_t *sf2,
                                     size_t ctx_col, size_t ctp_col,
                                     size_t ins_gap_min, size_t ins_gap_max,
                                     const SeqReadingPrefs prefs,
                                     SeqLoadingStats *stats)
{
  pool->rjob.ins_gap_min = ins_gap_min;
  pool->rjob.ins_gap_max = ins_gap_max;
  pool->rjob.ctx_col = ctx_col;
  pool->rjob.ctp_col = ctp_col;

  if(sf2 == NULL) {
    seq_parse_se_sf(sf1, &pool->rjob.r1, &pool->rjob.r2,
                    prefs, stats, &load_paths, pool);
  }
  else {
    pool->seen_pe = true;
    seq_parse_pe_sf(sf1, sf2, &pool->rjob.r1, &pool->rjob.r2,
                    prefs, stats, &load_paths, pool);
  }
}
