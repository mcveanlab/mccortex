#ifndef ASYNC_READ_IO_H_
#define ASYNC_READ_IO_H_

#include "seq_file.h"
#include "msg-pool/msgpool.h"

#include "loading_stats.h"

typedef struct
{
  seq_file_t *file1, *file2;
  void *ptr; // general porpoise pointer is passes into AsyncIOData
  const uint8_t fq_offset;
  const bool interleaved; // if file1 is an interleaved PE file
} AsyncIOReadInput;

typedef struct
{
  read_t r1, r2;
  void *ptr;
  uint8_t fq_offset1, fq_offset2;
} AsyncIOData;

#define async_task_pe_output(a) ((a)->file2 != NULL || (a)->interleaved)

// if out_base != NULL, we expect an output string as well:
//   -1, --seq <in>:<out>
//   -2, --seq2 <in1>:<in2>:<out>
//   -i, --seqi <in>:<out>
// if out_base == NULL, we expect:
//   -1, --seq <in>
//   -2, --seq2 <in1>:<in2>
//   -i, --seqi <in>
// If `out_base` is != NULL, it is set to point to the <out> string
void asyncio_task_parse(AsyncIOReadInput *task, char shortopt, char *path_arg,
                        uint8_t fq_offset, char **out_base);

void asyncio_task_close(AsyncIOReadInput *task);

void asynciodata_alloc(AsyncIOData *iod);
void asynciodata_dealloc(AsyncIOData *iod);

typedef struct AsyncIOWorker AsyncIOWorker;

void asynciodata_pool_init(void *el, size_t idx, void *args);
void asynciodata_pool_destroy(void *el, size_t idx, void *args);

void asyncio_run_threads(MsgPool *pool,
                         AsyncIOReadInput *asyncio_tasks, size_t num_inputs,
                         void (*job)(void*),
                         void *args, size_t num_readers, size_t elsize);

#endif /* ASYNC_READ_IO_H_ */
