#include "global.h"
#include "async_read_io.h"
#include "seq_reader.h"

#include <pthread.h>

struct AsyncIOWorker
{
  pthread_t thread;
  MsgPool *const pool;
  AsyncIOReadTask task;
  AsyncIOData data;
  volatile size_t *const num_running;
};

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

static void asynciodata_pool_init(char *el, size_t idx, void *args)
{
  (void)idx; (void)args;
  AsyncIOData d;
  asynciodata_alloc(&d);
  memcpy(el, &d, sizeof(AsyncIOData));
}

static void asynciodata_pool_destroy(char *el, size_t idx, void *args)
{
  (void)idx; (void)args;
  AsyncIOData d;
  memcpy(&d, el, sizeof(AsyncIOData));
  asynciodata_dealloc(&d);
}

static void async_io_worker_alloc(AsyncIOWorker *wrkr,
                                  const AsyncIOReadTask *task,
                                  MsgPool *pool, volatile size_t *num_running)
{
  AsyncIOWorker tmp = {.pool = pool, .task = *task, .num_running = num_running};
  asynciodata_alloc(&tmp.data);
  *wrkr = tmp;
}

static void async_io_worker_dealloc(AsyncIOWorker *wrkr)
{
  asynciodata_dealloc(&wrkr->data);
}

static void add_to_pool(read_t *r1, read_t *r2,
                        uint8_t fq_offset1, uint8_t fq_offset2, void *ptr)
{
  (void)r1; (void)r2; // r1,r2 already point to reads in wrkr->data.r{1,2}

  AsyncIOWorker *wrkr = (AsyncIOWorker*)ptr;

  // Set up out temporary data struct
  wrkr->data.fq_offset1 = fq_offset1;
  wrkr->data.fq_offset2 = fq_offset2;
  wrkr->data.ptr = wrkr->task.ptr;

  assert(&wrkr->data.r1 == r1);
  assert(r2 == NULL || &wrkr->data.r2 == r2);

  // Swap tmp with empty in the pool
  AsyncIOData data;
  msgpool_write(wrkr->pool, &wrkr->data, &data);
  wrkr->data = data;
}

static void* async_io_reader(void *ptr)
{
  AsyncIOWorker *wrkr = (AsyncIOWorker*)ptr;
  AsyncIOReadTask *task = &wrkr->task;

  seq_parse_pe_sf(task->file1, task->file2, task->fq_offset,
                  &wrkr->data.r1, &wrkr->data.r2, add_to_pool, wrkr);

  // Check if we are the last thread to finish, if so close the pool
  if(__sync_sub_and_fetch(wrkr->num_running, 1) == 0)
    msgpool_close(wrkr->pool);

  pthread_exit(NULL);
}

// Start loading into a pool
// returns an array of AsyncIOWorker of length len_files, each is a running
// thread putting reading into the pool passed.
AsyncIOWorker* asyncio_read_start(MsgPool *pool,
                                  const AsyncIOReadTask *tasks,
                                  size_t num_tasks)
{
  assert(pool->elsize == sizeof(AsyncIOData));

  size_t i;
  int rc;

  // Initiate all reads in the pool
  msgpool_iterate(pool, asynciodata_pool_init, NULL);

  // Create workers
  AsyncIOWorker *workers = malloc(num_tasks * sizeof(AsyncIOWorker));

  // Keep a counter of how many threads are still running
  // last thread to finish closes the pool
  size_t *num_running = malloc2(sizeof(size_t));
  *num_running = num_tasks;

  for(i = 0; i < num_tasks; i++)
    async_io_worker_alloc(&workers[i], &tasks[i], pool, num_running);

  // Thread attribute joinable
  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  // Start threads
  for(i = 0; i < num_tasks; i++) {
    rc = pthread_create(&workers[i].thread, &thread_attr,
                        async_io_reader, (void*)&workers[i]);
    if(rc != 0) die("Creating thread failed");
  }

  // Finished with thread attribute
  if(pthread_attr_destroy(&thread_attr))
    warn("Bad return value when disposing of pthread_attr");

  return workers;
}

// Wait until the pool is empty
void asyncio_read_finish(AsyncIOWorker *workers, size_t num_workers)
{
  int rc;
  size_t i;

  // Wait for threads to finish
  for(i = 0; i < num_workers; i++) {
    rc = pthread_join(workers[i].thread, NULL);
    if(rc != 0) die("Joining thread failed");
  }

  MsgPool *pool = workers[0].pool;
  msgpool_close(pool);
  msgpool_wait_til_empty(pool);

  // Free memory in the pool
  msgpool_iterate(pool, asynciodata_pool_destroy, NULL);

  for(i = 0; i < num_workers; i++) async_io_worker_dealloc(&workers[i]);
  free(workers);
}

