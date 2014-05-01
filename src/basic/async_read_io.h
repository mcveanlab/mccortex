#ifndef ASYNC_READ_IO_H_
#define ASYNC_READ_IO_H_

#include "seq_file.h"
#include "msg-pool/msgpool.h"

#include "loading_stats.h"

typedef struct
{
  seq_file_t *const file1, *const file2;
  void *const ptr; // general porpoise pointer is passes into AsyncIOData
  const uint8_t fq_offset;
  const bool interleaved; // if file1 is an interleaved PE file
} AsyncIOReadTask;

typedef struct
{
  read_t r1, r2;
  void *ptr;
  uint8_t fq_offset1, fq_offset2;
} AsyncIOData;

void asyncio_task_close(AsyncIOReadTask *task);

void asynciodata_alloc(AsyncIOData *iod);
void asynciodata_dealloc(AsyncIOData *iod);

typedef struct AsyncIOWorker AsyncIOWorker;

void asynciodata_pool_init(void *el, size_t idx, void *args);
void asynciodata_pool_destroy(void *el, size_t idx, void *args);

void asyncio_run_threads(MsgPool *pool,
                         AsyncIOReadTask *asyncio_tasks, size_t num_inputs,
                         void (*job)(void*),
                         void *args, size_t num_readers, size_t elsize);

#endif /* ASYNC_READ_IO_H_ */
