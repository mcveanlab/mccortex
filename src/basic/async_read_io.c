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
  size_t *const num_running;
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

void asynciodata_pool_init(void *el, size_t idx, void *args)
{
  (void)idx; (void)args;
  // status("alloc: %zu %p", idx, el);
  AsyncIOData d;
  asynciodata_alloc(&d);
  memcpy(el, &d, sizeof(AsyncIOData));
}

void asynciodata_pool_destroy(void *el, size_t idx, void *args)
{
  (void)idx; (void)args;
  // status("destruct: %zu %p", idx, el);
  AsyncIOData d;
  memcpy(&d, el, sizeof(AsyncIOData));
  asynciodata_dealloc(&d);
}

static void async_io_worker_alloc(AsyncIOWorker *wrkr,
                                  const AsyncIOReadTask *task,
                                  MsgPool *pool, size_t *num_running)
{
  assert(pool->elsize == sizeof(AsyncIOData));
  AsyncIOWorker tmp = {.pool = pool, .task = *task, .num_running = num_running};
  asynciodata_alloc(&tmp.data);
  memcpy(wrkr, &tmp, sizeof(AsyncIOWorker));
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
  assert(r1 != r2);

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
  size_t n = __sync_sub_and_fetch((volatile size_t*)wrkr->num_running, 1);

  if(n == 0) {
    msgpool_close(wrkr->pool);
    free(wrkr->num_running);
  }

  pthread_exit(NULL);
}

// Start loading into a pool
// returns an array of AsyncIOWorker of length len_files, each is a running
// thread putting reading into the pool passed.
AsyncIOWorker* asyncio_read_start(MsgPool *pool,
                                  const AsyncIOReadTask *tasks,
                                  size_t num_tasks)
{
  size_t i;
  int rc;

  // Initiate all reads in the pool
  assert(pool->elsize == sizeof(AsyncIOData));

  // Create workers
  AsyncIOWorker *workers = malloc2(num_tasks * sizeof(AsyncIOWorker));

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

  for(i = 0; i < num_workers; i++) async_io_worker_dealloc(&workers[i]);
  free(workers);
}

void asyncio_run_threads(MsgPool *pool,
                         AsyncIOReadTask *asyncio_tasks, size_t num_inputs,
                         void* (*job)(void*),
                         void *args, size_t num_readers, size_t elsize)
{
  size_t i; int rc;

  if(!num_inputs) return;
  assert(num_readers > 0);

  status("[asyncio] Threads: %zu input %zu reading", num_inputs, num_readers);

  // Start async io reading
  AsyncIOWorker *asyncio_workers;
  asyncio_workers = asyncio_read_start(pool, asyncio_tasks, num_inputs);

  // Run the workers until the pool is closed
  // Thread attribute joinable
  pthread_attr_t thread_attr;
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  pthread_t *threads = malloc2(num_readers * sizeof(pthread_t));

  // Start threads
  for(i = 0; i < num_readers; i++) {
    rc = pthread_create(&threads[i], &thread_attr, job, args+i*elsize);
    if(rc != 0) die("Creating thread failed");
  }

  // Finished with thread attribute
  pthread_attr_destroy(&thread_attr);

  // Join threads
  for(i = 0; i < num_readers; i++) {
    rc = pthread_join(threads[i], NULL);
    if(rc != 0) die("Joining thread failed");
  }

  free(threads);

  // Finish with the async io (waits until queue is empty)
  asyncio_read_finish(asyncio_workers, num_readers);
}
