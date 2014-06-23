#ifndef READ_THREAD_CMD_H_
#define READ_THREAD_CMD_H_

#include "cmd_mem.h"
#include "graph_file_reader.h"
#include "gpath_reader.h"
#include "correct_aln_input.h"

//
// ctx_thread.c and ctx_correct.c use many of the same command line arguments
// Therefore this file provides the common functionality to both
//

struct ReadThreadCmdArgs
{
  size_t num_of_threads;
  struct MemArgs memargs;
  char *graph_path, *out_ctp_path;
  bool use_new_paths, clean_paths;
  char *dump_seq_sizes, *dump_mp_sizes;
  int clean_threshold; // 0 => no cleaning, -1 => auto
  size_t colour; // ctx_correct only

  GraphFileReader gfile;
  GPathFileBuffer gpfiles;
  CorrectAlnInputBuffer inputs;

  size_t max_gap_limit; // max of inputs[].crt_params.ins_gap_max
  size_t path_max_usedcols; // max colours in path files
};

#define READ_THREAD_CMD_ARGS_INIT {.num_of_threads = 0,                \
                                   .memargs = MEM_ARGS_INIT,           \
                                   .graph_path = NULL,                 \
                                   .out_ctp_path = NULL,               \
                                   .use_new_paths = false,             \
                                   .clean_paths = false,               \
                                   .dump_seq_sizes = NULL,             \
                                   .dump_mp_sizes = NULL,              \
                                   .clean_threshold = 0,               \
                                   .colour = 0,                        \
                                   .gfile = INIT_GRAPH_READER_MACRO,   \
                                   .gpfiles = OBJBUF_INIT,             \
                                   .inputs = OBJBUF_INIT,              \
                                   .max_gap_limit = 0,                 \
                                   .path_max_usedcols = 0}

void read_thread_args_alloc(struct ReadThreadCmdArgs *args);
void read_thread_args_dealloc(struct ReadThreadCmdArgs *args);

// If `correct_cmd` is true:
//  - require --seq <in>:<out> instead of just --seq <in>
//  - allow multiple colours to be loaded from one graph
// If `correct_cmd` is false:
//  - do not take <out> argument with sequence files (e.g. --seq <in>)
//  - do not allow a multicolour graph to be loaded
void read_thread_args_parse(struct ReadThreadCmdArgs *args,
                            int argc, char **argv,
                            const struct option *longopts, bool correct_cmd);

#endif /* READ_THREAD_CMD_H_ */
