#include "global.h"

#include "msgpool.h"
#include <pthread.h>

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "correct_alignment.h"
#include "correct_reads_input.h"
#include "async_read_io.h"
#include "loading_stats.h"
#include "seq_reader.h"

// DEV: add back lost bases from edges of contigs

const char correct_usage[] =
"usage: "CMD" correct [options] <input.ctx> [...]\n"
"  Pull out contigs from the graph, print statistics\n"
"\n"
"  Options:\n"
"    --memory <M>             How much memory to use\n"
"    --paths <in.ctp>         Path file to load [can specify multiple times]\n"
"    --col <colour>           Colour to thread through\n"
"    --seq <in> <out>         Correct reads from file (supports sam,bam,fq,*.gz)\n"
"    --seq2 <in1> <in2> <out> Correct PE reads\n"
"    --twoway                 Use two-way gap filling\n"
"    --oneway                 Use one-way gap filling\n"
"    --fq_threshold <fq>      FASTQ quality threshold\n"
"    --fq_offset <qual>       FASTQ quality score offset\n"
"    --cut_hp <N>             Cut reads afer <N> consecutive bases\n"
"    --FR --FF --RF           Mate pair orientation [default: FR]\n"
"\n"
" --seq outputs <out>.fa.gz, --seq2 outputs <out>.1.fa.gz, <out>.2.fa.gz\n"
" --seq must come AFTER two/oneway options. Output may be slightly shuffled.\n";


typedef struct {
  StrBuf out1, out2;
  gzFile gz1, gz2;
} CorrectedOutput;

typedef struct {
  dBGraph *db_graph;
  MsgPool *pool;
  dBAlignment aln;
  CorrectAlnWorker corrector;
  CorrectAlnReadsTask *input;
  LoadingStats stats;
  pthread_mutex_t lock;
} CorrectReadsWorker;

// Returns true on success, false on error
static boolean corrected_output_alloc(CorrectedOutput *out,
                                      CorrectAlnReadsTask *in)
{
  out->gz1 = out->gz2 = NULL;
  const char *base = (const char*)in->ptr;
  StrBuf *out1 = &out->out1, *out2 = &out->out2;
  strbuf_alloc(out1, 512);
  strbuf_alloc(out2, 512);

  if(in->file2 == NULL) {
    strbuf_append_str(out1, base);
    strbuf_append_str(out1, ".fa.gz");
    if(futil_file_exists(out1->buff)) {
      warn("File already exists: %s", out1->buff);
      return false;
    }
    out->gz1 = gzopen(out1->buff, "w");
    if(out->gz1 == NULL) warn("Cannot write to: %s", out1->buff);
    return (out->gz1 != NULL);
  } else {
    strbuf_append_str(out1, base);
    strbuf_append_str(out2, base);
    strbuf_append_str(out1, ".1.fa.gz");
    strbuf_append_str(out2, ".2.fa.gz");

    boolean fexists = false;
    if(futil_file_exists(out1->buff)) {
      warn("File already exists: %s", out1->buff);
      fexists = true;
    }
    if(futil_file_exists(out2->buff)) {
      warn("File already exists: %s", out2->buff);
      fexists = true;
    }
    if(fexists) return true;

    out->gz1 = gzopen(out1->buff, "w");
    out->gz2 = gzopen(out2->buff, "w");
    if(out->gz1 == NULL) warn("Cannot write to: %s", out1->buff);
    if(out->gz2 == NULL) warn("Cannot write to: %s", out2->buff);
    return (out->gz1 != NULL && out->gz2 != NULL);
  }
}

static void corrected_output_clean_up(CorrectedOutput *out)
{
  if(out->gz1 != NULL) { gzclose(out->gz1); unlink(out->out1.buff); }
  if(out->gz2 != NULL) { gzclose(out->gz2); unlink(out->out2.buff); }
}

static void correct_reads_worker_alloc(CorrectReadsWorker *wrkr,
                                       MsgPool *pool,  dBGraph *db_graph)
{
  wrkr->db_graph = db_graph;
  wrkr->pool = pool;
  wrkr->input = NULL;
  db_alignment_alloc(&wrkr->aln);
  correct_aln_worker_alloc(&wrkr->corrector, db_graph);
  loading_stats_init(&wrkr->stats);
  if(pthread_mutex_init(&wrkr->lock, NULL) != 0) die("Mutex init failed");
}

static void correct_reads_worker_dealloc(CorrectReadsWorker *wrkr)
{
  pthread_mutex_destroy(&wrkr->lock);
  db_alignment_dealloc(&wrkr->aln);
  correct_aln_worker_dealloc(&wrkr->corrector);
}

static boolean handle_read(CorrectReadsWorker *wrkr, gzFile gz,
                           read_t *r, uint8_t fq_cutoff, uint8_t hp_cutoff)
{
  dBNodeBuffer *nbuf;
  boolean output = false;

  db_alignment_from_reads(&wrkr->aln, r, NULL,
                          fq_cutoff, 0, hp_cutoff,
                          wrkr->db_graph);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, wrkr->input->crt_params);

  while((nbuf = correct_alignment_nxt(&wrkr->corrector)) != NULL) {
    if(!output) gzprintf(gz, ">%s\n", r->name.b);
    else gzprintf(gz, "N");
    db_nodes_gzprint(nbuf->data, nbuf->len, wrkr->db_graph, gz);
    output = true;
  }
  if(output) gzprintf(gz, "\n");

  return output;
}

static void correct_read(CorrectReadsWorker *wrkr, AsyncIOData *data)
{
  uint8_t fq_cutoff1, fq_cutoff2, hp_cutoff;

  CorrectAlnReadsTask *input = wrkr->input;
  CorrectedOutput *output = (CorrectedOutput*)input->ptr;

  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  fq_cutoff1 = fq_cutoff2 = input->fq_cutoff;

  if(fq_cutoff1 > 0) {
    fq_cutoff1 += data->fq_offset1;
    fq_cutoff2 += data->fq_offset2;
  }

  hp_cutoff = input->hp_cutoff;

  // Update stats
  if(r2 == NULL) wrkr->stats.num_se_reads++;
  else wrkr->stats.num_pe_reads += 2;

  boolean print1, print2 = false;
  gzFile gz1 = output->gz1, gz2 = (output->gz2 ? output->gz2 : output->gz1);

  // Don't want two threads writing to the same file at the same time: get mutex
  pthread_mutex_lock(&wrkr->lock);

  print1 = handle_read(wrkr, gz1, r1, fq_cutoff1, hp_cutoff);
  print2 = (r2 != NULL && handle_read(wrkr, gz2, r2, fq_cutoff2, hp_cutoff));

  if(r2 != NULL && print1 != print2) {
    if(print1) gzprintf(gz2, ">%s\n\n", r2->name.b);
    else gzprintf(gz1, ">%s\n\n", r1->name.b);
  }

  // release mutex
  pthread_mutex_unlock(&wrkr->lock);
}

// pthread method, loop: grabs job, does processing
static void* correct_reads_thread(void *ptr)
{
  CorrectReadsWorker *wrkr = (CorrectReadsWorker*)ptr;

  // incoming is empty, outgoing is allocated
  AsyncIOData data0, data1, *incoming = &data0, *outgoing = &data1, *tmpdataptr;
  asynciodata_alloc(&data1);

  // Read into incoming, write out outgoing
  while(msgpool_read(wrkr->pool, incoming, outgoing))
  {
    correct_read(wrkr, incoming);
    // Now incoming is allocated, outgoing is not, so switch
    SWAP(incoming, outgoing, tmpdataptr);
  }

  asynciodata_dealloc(outgoing);

  pthread_exit(NULL);
}


int ctx_correct(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 2 arguments

  size_t max_ninputs = argc / 2, num_inputs = 0;
  CorrectAlnReadsTask *inputs = malloc2(max_ninputs * sizeof(CorrectAlnReadsTask));
  size_t i, j, num_work_threads = args->max_work_threads;
  int argi; // arg index to continue from

  // Load args
  // argi = load_args(argc, argv, inputs, &num_inputs);
  argi = correct_reads_parse(argc, argv, inputs, num_inputs, false, true);

  if(argi == argc) cmd_print_usage("Expected at least one graph file");
  size_t num_gfiles = argc - argi;
  char **graph_paths = &argv[argi];

  //
  // Open graph gfiles
  //
  GraphFileReader gfiles[num_gfiles];
  size_t ctx_max_kmers = 0, ctx_total_cols = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], graph_paths[i], true);

    // stack graphs
    file_filter_update_intocol(&gfiles[i].fltr,
                               gfiles[i].fltr.intocol + ctx_total_cols);

    ctx_total_cols = graph_file_usedcols(&gfiles[i]);
    ctx_max_kmers = MAX2(ctx_max_kmers, gfiles[i].hdr.num_of_kmers);
  }

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t ctp_max_mem = 0, ctp_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    ctp_max_mem = MAX2(ctp_max_mem, pfiles[i].hdr.num_path_bytes);
    ctp_max_usedcols = MAX2(ctp_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(gfiles, num_gfiles, pfiles, num_pfiles);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;
  char path_mem_str[100];

  // 1 bit needed per kmer if we need to keep track of noreseed
  bits_per_kmer = sizeof(Edges)*8 + ctx_max_kmers + sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, false, &graph_mem);

  // Paths memory
  size_t tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem = ctp_max_mem + tmp_path_mem;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s", path_mem_str);

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Check we can read all output files
  //
  // Loop over inputs, change input->ptr from char* to CorrectedOutput
  CorrectedOutput outputs[num_inputs];

  for(i = 0; i < num_inputs && corrected_output_alloc(&outputs[i], &inputs[i]); i++)
    inputs[i].ptr = &outputs[i];

  // Check if something went wrong
  if(i < num_inputs) {
    for(j = 0; j < i; j++)
      corrected_output_clean_up(&outputs[i]);
    die("Couldn't open output files");
  }

  // futil_is_file_writable() creates the output file - we have to delete these
  // if we then fail
  for(i = 0; i < num_inputs; i++) {
    if(futil_file_exists((char*)inputs[i].ptr)) {
      for(j = 0; j < i; j++) unlink((char*)inputs[j].ptr);
      die("File already exists: %s", (char*)inputs[i].ptr);
    }
    if(!futil_is_file_writable((char*)inputs[i].ptr)) {
      for(j = 0; j < i; j++) unlink((char*)inputs[j].ptr);
      die("Cannot write to output file: %s", (char*)inputs[i].ptr);
    }
  }

  // Allocate
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ctx_total_cols, 1, kmers_in_hash);

  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);

  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = calloc2(bytes_per_col * ctx_total_cols, sizeof(uint8_t));
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(PathIndex));
  memset(db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(PathIndex));

  path_store_alloc(&db_graph.pdata, ctp_max_mem, tmp_path_mem, ctp_max_usedcols);

  // Run alignment
  MsgPool pool;
  msgpool_alloc_spinlock(&pool, MSGPOOLRSIZE, sizeof(AsyncIOData));
  msgpool_iterate(&pool, asynciodata_pool_init, NULL);

  CorrectReadsWorker *wrkrs = malloc2(num_work_threads * sizeof(CorrectReadsWorker));

  for(i = 0; i < num_work_threads; i++)
    correct_reads_worker_alloc(&wrkrs[i], &pool, &db_graph);

  AsyncIOReadTask *asyncio_tasks = malloc2(num_inputs * sizeof(AsyncIOReadTask));
  correct_reads_input_to_asycio(asyncio_tasks, inputs, num_inputs);

  // Load input files num_io_threads at a time
  size_t num_io_threads = args->max_io_threads, n;

  for(i = 0; i < num_inputs; i += num_io_threads) {
    n = MIN2(num_inputs - i, num_io_threads);
    asyncio_run_threads(&pool, asyncio_tasks+i, n,
                        correct_reads_thread, wrkrs, num_work_threads);
  }

  for(i = 0; i < num_work_threads; i++)
    correct_reads_worker_dealloc(&wrkrs[i]);

  free(wrkrs);
  free(asyncio_tasks);
  msgpool_iterate(&pool, asynciodata_pool_destroy, NULL);
  msgpool_dealloc(&pool);

  for(i = 0; i < num_inputs; i++)
    corrected_output_alloc(&outputs[i], &inputs[i]);

  free(inputs);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);
  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
