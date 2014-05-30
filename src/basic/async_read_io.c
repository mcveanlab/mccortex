#include "global.h"
#include "async_read_io.h"
#include "seq_reader.h"
#include "util.h" // util_run_threads()

#include <pthread.h>

struct AsyncIOWorker
{
  pthread_t thread;
  MsgPool *const pool;
  AsyncIOReadTask task;
  size_t *const num_running;
};

// path can be a single path or a pair of paths separated by a comma or colon
void asyncio_task_parse(AsyncIOReadTask *task, char shortopt, char *path,
                        uint8_t fq_offset)
{
  seq_file_t *sf1 = NULL, *sf2 = NULL;
  bool se = false, pe = false, il = false;
  se = (shortopt == '1');
  pe = (shortopt == '2');
  il = (shortopt == 'i');
  if(!se && !pe && !il) die("Unknown command argument: ");
  if(se || il) {
    if((sf1 = seq_open(path)) == NULL)
      die("Cannot open -%c file: %s", shortopt, path);
  }
  if(pe) {
    char *sep = strchr(path,':'), *path2;
    if(sep == NULL) sep = strchr(path,',');
    if(sep == NULL || sep == path+1 || sep[1] == '\0')
      die("%s <in1:in2> expects paths separated by a colon ':'", path);

    *sep = '\0';
    path2 = sep+1;

    if((sf1 = seq_open(path)) == NULL)
      die("Cannot open %c file: %s", shortopt, path);
    if((sf2 = seq_open(path2)) == NULL)
      die("Cannot open %c file: %s", shortopt, path2);
  }

  AsyncIOReadTask tmp = {.file1 = sf1, .file2 = sf2,
                         .fq_offset = fq_offset, .interleaved = il,
                         .ptr = NULL};
  memcpy(task, &tmp, sizeof(AsyncIOReadTask));
}

void asyncio_task_close(AsyncIOReadTask *task)
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
                                 const AsyncIOReadTask *task,
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
  AsyncIOReadTask *task = &wrkr->task;

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
                                         const AsyncIOReadTask *tasks,
                                         size_t num_tasks)
{
  if(num_tasks == 0) return NULL;

  size_t i;
  int rc;

  // Initiate all reads in the pool
  ctx_assert(pool->elsize == sizeof(AsyncIOData*));

  // Create workers
  AsyncIOWorker *workers = ctx_malloc(num_tasks * sizeof(AsyncIOWorker));

  // Keep a counter of how many threads are still running
  // last thread to finish closes the pool
  size_t *num_running = ctx_malloc(sizeof(size_t));
  *num_running = num_tasks;

  for(i = 0; i < num_tasks; i++)
    async_io_worker_init(&workers[i], &tasks[i], pool, num_running);

  // Start threads
  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  for(i = 0; i < num_tasks; i++) {
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
                         AsyncIOReadTask *asyncio_tasks, size_t num_inputs,
                         void (*job)(void*),
                         void *args, size_t num_readers, size_t elsize)
{
  if(!num_inputs) return;
  ctx_assert(num_readers > 0);

  status("[asyncio] Inputs: %zu; Threads: %zu", num_inputs, num_readers);

  // Start async io reading
  AsyncIOWorker *asyncio_workers;
  asyncio_workers = asyncio_read_start(pool, asyncio_tasks, num_inputs);

  util_run_threads(args, num_readers, elsize, num_readers, job);

  // Finish with the async io (waits until queue is empty)
  asyncio_read_finish(asyncio_workers, num_inputs);
}
