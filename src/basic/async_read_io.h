#ifndef ASYNC_READ_IO_H_
#define ASYNC_READ_IO_H_

#include "seq_file.h"
#include "msgpool.h"

#include "loading_stats.h"

typedef struct
{
  seq_file_t *const file1, *const file2;
  SeqLoadingStats *const stats; // stats are written to here
  void *const ptr;
  const uint8_t fq_offset;
} AsyncIOReadTask;

typedef struct
{
  read_t r1, r2;
  SeqLoadingStats *stats;
  void *ptr;
  uint8_t fq_offset1, fq_offset2;
} AsyncIOData;

void asynciodata_alloc(AsyncIOData *iod);
void asynciodata_dealloc(AsyncIOData *iod);

typedef struct AsyncIOWorker AsyncIOWorker;

// Pool is of type AsyncIOData
// Start loading into a pool
AsyncIOWorker* asyncio_read_start(MsgPool *pool,
                                  const AsyncIOReadTask *files,
                                  size_t num_files);

// Function blocks until all reads are loaded into the pool and the pool is empty
// Wait until the pool is empty
// frees workers memory
// Each worker is a thread reading from a file (or pair of files)
void asyncio_read_finish(AsyncIOWorker *workers, size_t num_workers);

#endif /* ASYNC_READ_IO_H_ */
