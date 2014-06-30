#include "global.h"
#include "async_read_io.h"
#include "seq_reader.h"
#include "file_util.h"
#include "util.h" // util_run_threads()

#include <pthread.h>

struct AsyncIOWorker
{
  pthread_t thread;
  MsgPool *const pool;
  AsyncIOInput task;
  size_t *const num_running;
};


// if out_base != NULL, we expect an output string as well:
//   -1, --seq <in>:<out>
//   -2, --seq2 <in1>:<in2>:<out>
//   -i, --seqi <in>:<out>
// if out_base == NULL, we expect:
//   -1, --seq <in>
//   -2, --seq2 <in1>:<in2>
//   -i, --seqi <in>
// If `out_base` is != NULL, it is set to point to the <out> string
void asyncio_task_parse(AsyncIOInput *task, char shortopt, char *path_arg,
                        uint8_t fq_offset, char **out_base)
{
  ctx_assert(shortopt == '1' || shortopt == '2' || shortopt == 'i');

  seq_file_t *sf1 = NULL, *sf2 = NULL;
  bool se = (shortopt == '1');
  bool pe = (shortopt == '2');
  bool il = (shortopt == 'i');

  if(!se && !pe && !il)
    die("Unknown command argument: -%c", shortopt);

  char *paths[3];
  size_t npaths, exp_npaths;

  npaths = string_split_str(path_arg, ':', paths, 3);
  if(npaths == 1)
    npaths = string_split_str(path_arg, ',', paths, 3);

  exp_npaths = (pe ? 2 : 1) + (out_base != NULL ? 1 : 0);
  if(npaths != exp_npaths) {
    const char *in_str = (se || il) ? "<in>" : "<in1>:<in2>";
    const char *out_str = out_base == NULL ? "" : ":<out>";
    die("Expected -%c %s%s", shortopt, in_str, out_str);
  }

  if(out_base != NULL)
    out_base[0] = paths[pe ? 2 : 1];

  if(pe) {
    if((sf1 = seq_open(paths[0])) == NULL)
      die("Cannot open %c file: %s", shortopt, paths[0]);
    if((sf2 = seq_open(paths[1])) == NULL)
      die("Cannot open %c file: %s", shortopt, paths[1]);
  }
  else {
    if((sf1 = seq_open(paths[0])) == NULL)
      die("Cannot open -%c file: %s", shortopt, paths[0]);
  }

  AsyncIOInput tmp = {.file1 = sf1, .file2 = sf2,
                         .fq_offset = fq_offset, .interleaved = il,
                         .ptr = NULL};
  memcpy(task, &tmp, sizeof(AsyncIOInput));
}

void asyncio_task_close(AsyncIOInput *task)
{
  if(task->file1 != NULL) seq_close(task->file1);
  if(task->file2 != NULL) seq_close(task->file2);
}

void asynciodata_alloc(AsyncIOData *iod)
{
  if(seq_read_alloc(&iod->r1) == NULL ||
     seq_read_alloc(&iod->r2) == NULL) die("Out of memory");
}

void asynciodata_dealloc(AsyncIOData *iod)
{
  seq_read_dealloc(&iod->r1);
  seq_read_dealloc(&iod->r2);
}

void asynciodata_pool_init(void *el, size_t idx, void *args)
{
  // status("alloc: %zu %p", idx, el);
  AsyncIOData *store = (AsyncIOData*)args, *data = store + idx;
  memcpy(el, &data, sizeof(AsyncIOData*));
}

// No memory allocated for io worker
static void async_io_worker_init(AsyncIOWorker *wrkr,
                                 const AsyncIOInput *task,
                                 MsgPool *pool, size_t *num_running)
{
  ctx_assert(pool->elsize == sizeof(AsyncIOData*));
  AsyncIOWorker tmp = {.pool = pool, .task = *task, .num_running = num_running};
  memcpy(wrkr, &tmp, sizeof(AsyncIOWorker));
}

static void add_to_pool(read_t *r1, read_t *r2,
                        uint8_t fq_offset1, uint8_t fq_offset2,
                        void *ptr)
{
  AsyncIOWorker *wrkr = (AsyncIOWorker*)ptr;
  MsgPool *pool = wrkr->pool;
  int pos;
  AsyncIOData *data;

  // Swap reads and parameters into the data obj
  pos = msgpool_claim_write(pool);
  memcpy(&data, msgpool_get_ptr(pool, pos), sizeof(AsyncIOData*));

  data->fq_offset1 = fq_offset1;
  data->fq_offset2 = fq_offset2;
  data->ptr = wrkr->task.ptr;

  SWAP(data->r1, *r1);

  if(r2) SWAP(data->r2, *r2);
  else seq_read_reset(&data->r2);

  msgpool_release(pool, pos, MPOOL_FULL);
}

static void* async_io_reader(void *ptr) __attribute__((noreturn));

static void* async_io_reader(void *ptr)
{
  AsyncIOWorker *wrkr = (AsyncIOWorker*)ptr;
  AsyncIOInput *task = &wrkr->task;

  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

  if(task->interleaved)
  {
    seq_parse_interleaved_sf(task->file1, task->fq_offset,
                             &r1, &r2, add_to_pool, wrkr);
  } else {
    seq_parse_pe_sf(task->file1, task->file2, task->fq_offset,
                    &r1, &r2, add_to_pool, wrkr);
  }

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  // Check if we are the last thread to finish, if so close the pool
  size_t n = __sync_sub_and_fetch((volatile size_t*)wrkr->num_running, 1);

  if(n == 0) {
    msgpool_close(wrkr->pool);
    ctx_free(wrkr->num_running);
  }

  pthread_exit(NULL);
}

// Start loading into a pool
// returns an array of AsyncIOWorker of length len_files, each is a running
// thread putting reading into the pool passed.
static AsyncIOWorker* asyncio_read_start(MsgPool *pool,
                                         const AsyncIOInput *inputs,
                                         size_t num_inputs)
{
  if(num_inputs == 0) return NULL;

  size_t i;
  int rc;

  // Initiate all reads in the pool
  ctx_assert(pool->elsize == sizeof(AsyncIOData*));

  // Create workers
  AsyncIOWorker *workers = ctx_malloc(num_inputs * sizeof(AsyncIOWorker));

  // Keep a counter of how many threads are still running
  // last thread to finish closes the pool
  size_t *num_running = ctx_malloc(sizeof(size_t));
  *num_running = num_inputs;

  for(i = 0; i < num_inputs; i++)
    async_io_worker_init(&workers[i], &inputs[i], pool, num_running);

  // Start threads
  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  for(i = 0; i < num_inputs; i++) {
    rc = pthread_create(&workers[i].thread, &thread_attr,
                        async_io_reader, (void*)&workers[i]);
    if(rc != 0) die("Creating thread failed: %s", strerror(rc));
  }

  pthread_attr_destroy(&thread_attr);
  return workers;
}

// Wait until the pool is empty
static void asyncio_read_finish(AsyncIOWorker *workers, size_t num_workers)
{
  if(num_workers == 0) return;

  size_t i;
  int rc;

  // Wait for threads to finish
  for(i = 0; i < num_workers; i++) {
    rc = pthread_join(workers[i].thread, NULL);
    if(rc != 0) die("Joining thread failed: %s", strerror(rc));
  }

  MsgPool *pool = workers[0].pool;
  msgpool_close(pool);
  msgpool_wait_til_empty(pool);
  ctx_assert(pool->num_full == 0);

  ctx_free(workers);
}

void asyncio_run_threads(MsgPool *pool,
                         AsyncIOInput *asyncio_inputs, size_t num_inputs,
                         void (*job)(void*),
                         void *args, size_t num_readers, size_t elsize)
{
  if(!num_inputs) return;
  ctx_assert(num_readers > 0);

  status("[asyncio] Inputs: %zu; Threads: %zu", num_inputs, num_readers);

  // Start async io reading
  AsyncIOWorker *asyncio_workers;
  asyncio_workers = asyncio_read_start(pool, asyncio_inputs, num_inputs);

  util_run_threads(args, num_readers, elsize, num_readers, job);

  // Finish with the async io (waits until queue is empty)
  asyncio_read_finish(asyncio_workers, num_inputs);
}

typedef struct {
  MsgPool *pool;
  void (*func)(AsyncIOData *_data, void *_arg);
  void *arg;
} PoolFuncPair;

// pthread method, loop: reads from pool, call function
static void grab_reads_from_pool(void *arg)
{
  PoolFuncPair wrkr = *(PoolFuncPair*)arg;
  int pos;
  AsyncIOData *data = NULL;

  while((pos = msgpool_claim_read(wrkr.pool)) != -1)
  {
    memcpy(&data, msgpool_get_ptr(wrkr.pool, pos), sizeof(AsyncIOData*));
    wrkr.func(data, wrkr.arg);
    msgpool_release(wrkr.pool, pos, MPOOL_EMPTY);
  }
}

// `num_inputs` number of threads pushing reads into the pool
// `num_readers` number of threads pulling reads from the pool
void asyncio_run_pool(AsyncIOInput *asyncio_inputs, size_t num_inputs,
                      void (*job)(AsyncIOData *_data, void *_arg),
                      void *args, size_t num_readers, size_t elsize)
{
  size_t i;
  AsyncIOData *data = ctx_malloc(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  MsgPool pool;
  msgpool_alloc(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*), USE_MSG_POOL);
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  PoolFuncPair poolfunc[num_readers];

  for(i = 0; i < num_readers; i++) {
    poolfunc[i] = (PoolFuncPair){.pool = &pool, .func = job,
                                 .arg = (char*)args+i*elsize};
  }

  asyncio_run_threads(&pool, asyncio_inputs, num_inputs, grab_reads_from_pool,
                      &poolfunc, num_readers, sizeof(PoolFuncPair));

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  ctx_free(data);
  msgpool_dealloc(&pool);
}

// Guess numer of kmers
size_t asyncio_input_nkmers(const AsyncIOInput *io)
{
  size_t i, est_num_bases = 0;
  for(i = 0; i < 2; i++) {
    seq_file_t *sf = i ? io->file1 : io->file2;
    if(sf) {
      off_t fsize = futil_get_file_size(sf->path);
      if(fsize < 0) {
        warn("Cannot get file size: %s", sf->path);
        return SIZE_MAX;
      }
      else {
        if(seq_is_fastq(sf) || seq_is_sam(sf))
          est_num_bases += fsize / 2;
        else
          est_num_bases += fsize;
      }
    }
  }
  // Heuristic: multiply by 5 to correct for compression
  return est_num_bases * 5;
}
