#include "global.h"

#include "msgpool.h"
#include <pthread.h>

#include "commands.h"
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

// DEV: print read sequence in lower case instead of N

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

#define MAX_CONTEXT 50

typedef struct {
  StrBuf out1, out2;
  gzFile gz1, gz2;
  pthread_mutex_t lock;
} CorrectedOutput;

typedef struct {
  dBGraph *db_graph;
  MsgPool *pool;
  dBAlignment aln;
  CorrectAlnWorker corrector;
  CorrectAlnReadsTask *input;
  LoadingStats stats;
  // For filling in gaps
  GraphWalker wlk;
  RepeatWalker rptwlk;
  StrBuf buf1, buf2;
  dBNodeBuffer tmpnbuf;
} CorrectReadsWorker;

// Returns true on success, false on error
static bool corrected_output_open(CorrectedOutput *out,
                                     CorrectAlnReadsTask *in)
{
  out->gz1 = out->gz2 = NULL;
  const char *base = (const char*)in->ptr;
  StrBuf *out1 = &out->out1, *out2 = &out->out2;
  strbuf_alloc(out1, 512);
  strbuf_alloc(out2, 512);

  if(pthread_mutex_init(&out->lock, NULL) != 0) die("Mutex init failed");

  if(in->file2 == NULL) {
    strbuf_append_str(out1, base);
    strbuf_append_str(out1, ".fa.gz");
    if(futil_file_exists(out1->buff)) {
      warn("File already exists: %s", out1->buff);
      return false;
    }
    out->gz1 = gzopen(out1->buff, "w");
    if(out->gz1 == NULL) {
      warn("Cannot write to: %s", out1->buff);
      return false;
    }
  } else {
    strbuf_append_str(out1, base);
    strbuf_append_str(out2, base);
    strbuf_append_str(out1, ".1.fa.gz");
    strbuf_append_str(out2, ".2.fa.gz");

    bool fexists = false;
    if(futil_file_exists(out1->buff)) {
      warn("File already exists: %s", out1->buff);
      fexists = true;
    }
    if(futil_file_exists(out2->buff)) {
      warn("File already exists: %s", out2->buff);
      fexists = true;
    }
    if(fexists) return false;

    out->gz1 = gzopen(out1->buff, "w");
    out->gz2 = gzopen(out2->buff, "w");
    if(out->gz1 == NULL) warn("Cannot write to: %s", out1->buff);
    if(out->gz2 == NULL) warn("Cannot write to: %s", out2->buff);
    if(out->gz1 == NULL || out->gz2 == NULL) return false;
  }

  return true;
}

static void corrected_output_close(CorrectedOutput *out)
{
  if(out->gz1 != NULL) gzclose(out->gz1);
  if(out->gz2 != NULL) gzclose(out->gz2);
  out->gz1 = out->gz2 = NULL;
  strbuf_dealloc(&out->out1);
  strbuf_dealloc(&out->out2);
  pthread_mutex_destroy(&out->lock);
}

static void corrected_output_delete(CorrectedOutput *out)
{
  if(out->gz1 != NULL) unlink(out->out1.buff);
  if(out->gz2 != NULL) unlink(out->out2.buff);
  corrected_output_close(out);
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
  graph_walker_alloc(&wrkr->wlk);
  rpt_walker_alloc(&wrkr->rptwlk, db_graph->ht.capacity, 22); // 4MB bloom
  strbuf_alloc(&wrkr->buf1, 1024);
  strbuf_alloc(&wrkr->buf2, 1024);
  db_node_buf_alloc(&wrkr->tmpnbuf, 512);
}

static void correct_reads_worker_dealloc(CorrectReadsWorker *wrkr)
{
  db_alignment_dealloc(&wrkr->aln);
  correct_aln_worker_dealloc(&wrkr->corrector);
  graph_walker_dealloc(&wrkr->wlk);
  rpt_walker_dealloc(&wrkr->rptwlk);
  strbuf_dealloc(&wrkr->buf1);
  strbuf_dealloc(&wrkr->buf2);
  db_node_buf_dealloc(&wrkr->tmpnbuf);
}

static void handle_read(CorrectReadsWorker *wrkr, const read_t *r, StrBuf *buf,
                        uint8_t fq_cutoff, uint8_t hp_cutoff)
{
  dBNodeBuffer *nbuf, *tmpnbuf = &wrkr->tmpnbuf;
  dBAlignment *aln = &wrkr->aln;
  GraphWalker *wlk = &wrkr->wlk;
  RepeatWalker *rptwlk = &wrkr->rptwlk;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = db_graph->kmer_size;

  size_t i, idx, gap, num_n, nbases;
  size_t init_len, end_len;
  BinaryKmer bkmer;
  Nucleotide nuc;
  char bkmerstr[MAX_KMER_SIZE+1];

  db_alignment_from_reads(&wrkr->aln, r, NULL,
                          fq_cutoff, 0, hp_cutoff, db_graph);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, wrkr->input->crt_params);

  // Get first alignment
  nbuf = correct_alignment_nxt(&wrkr->corrector);

  // Extend left
  strbuf_reset(buf);
  strbuf_append_char(buf, '>');
  strbuf_append_strn(buf, r->name.b, r->name.end);
  strbuf_append_char(buf, '\n');

  if(nbuf == NULL) {
    for(i = 0; i < r->seq.end; i++) strbuf_append_char(buf, 'N');
    strbuf_append_char(buf, '\n');
    return;
  }

  // extend left
  size_t left_gap = aln->gaps.data[0], right_gap = aln->r1enderr;

  if(left_gap > 0)
  {
    // Walk left
    graph_walker_prime(wlk, nbuf->data, nbuf->len, MAX_CONTEXT, false,
                       0, 0, db_graph);

    db_node_buf_reset(tmpnbuf);
    db_node_buf_ensure_capacity(nbuf, left_gap);

    while(tmpnbuf->len < left_gap && graph_traverse(wlk) &&
          rpt_walker_attempt_traverse(rptwlk, wlk))
    {
      tmpnbuf->data[tmpnbuf->len++] = wlk->node;
    }

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, tmpnbuf->data, tmpnbuf->len);

    // Add Ns for bases we couldn't resolve
    for(i = tmpnbuf->len; i < left_gap; i++) strbuf_append_char(buf, 'N');

    // Append bases
    for(i = tmpnbuf->len-1; i != SIZE_MAX; i--) {
      nuc = db_node_get_first_nuc(db_node_reverse(tmpnbuf->data[i]), db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }
  }

  // Append first contig
  bkmer = db_graph_oriented_bkmer(db_graph, nbuf->data[0]);
  binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
  strbuf_append_strn(buf, bkmerstr, kmer_size);
  for(i = 1; i < nbuf->len; i++) {
    nuc = db_node_get_last_nuc(nbuf->data[i], db_graph);
    strbuf_append_char(buf, dna_nuc_to_char(nuc));
  }

  while(correct_alignment_get_endidx(&wrkr->corrector) < aln->nodes.len)
  {
    nbuf = correct_alignment_nxt(&wrkr->corrector);
    ctx_assert(nbuf != NULL);
    idx = correct_alignment_get_strtidx(&wrkr->corrector);
    gap = aln->gaps.data[idx];
    num_n = gap < kmer_size ? 0 : gap - kmer_size + 1;
    for(i = 0; i < num_n; i++) strbuf_append_char(buf, 'N');

    nbases = MIN2(gap+1, kmer_size);
    binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
    strbuf_append_strn(buf, bkmerstr+kmer_size-nbases, nbases);

    for(i = 1; i < nbuf->len; i++) {
      nuc = db_node_get_last_nuc(nbuf->data[i], db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }
  }

  // extend right
  if(right_gap > 0)
  {
    // walk right
    graph_walker_prime(wlk, nbuf->data, nbuf->len, MAX_CONTEXT, true,
                       0, 0, db_graph);

    init_len = nbuf->len;
    end_len = init_len + right_gap;
    db_node_buf_ensure_capacity(nbuf, end_len);

    while(nbuf->len < end_len && graph_traverse(wlk) &&
          rpt_walker_attempt_traverse(rptwlk, wlk))
    {
      nbuf->data[nbuf->len++] = wlk->node;
    }

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, nbuf->data+init_len, nbuf->len-init_len);

    // Copy added bases into buffer
    for(i = init_len; i < nbuf->len; i++) {
      nuc = db_node_get_first_nuc(nbuf->data[i], db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }

    // Add Ns for bases we couldn't resolve
    for(i = nbuf->len; i < end_len; i++) strbuf_append_char(buf, 'N');
  }

  strbuf_append_char(buf, '\n');
}

static void correct_read(CorrectReadsWorker *wrkr, AsyncIOData *data)
{
  uint8_t fq_cutoff1, fq_cutoff2, hp_cutoff;

  CorrectAlnReadsTask *input = wrkr->input;
  CorrectedOutput *output = (CorrectedOutput*)input->ptr;
  StrBuf *buf1 = &wrkr->buf1, *buf2 = &wrkr->buf2;

  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  fq_cutoff1 = fq_cutoff2 = input->fq_cutoff;

  if(fq_cutoff1 > 0) {
    fq_cutoff1 += data->fq_offset1;
    fq_cutoff2 += data->fq_offset2;
  }

  hp_cutoff = input->hp_cutoff;

  strbuf_reset(buf1);
  strbuf_reset(buf2);

  // Update stats
  if(r2 == NULL) wrkr->stats.num_se_reads++;
  else wrkr->stats.num_pe_reads += 2;

  gzFile gz1 = output->gz1, gz2 = (output->gz2 ? output->gz2 : output->gz1);

  handle_read(wrkr, r1, buf1, fq_cutoff1, hp_cutoff);
  if(r2 != NULL) handle_read(wrkr, r2, buf2, fq_cutoff2, hp_cutoff);

  // Don't want two threads writing to the same file at the same time: get mutex
  pthread_mutex_lock(&output->lock);

  gzputs(gz1, buf1->buff);
  if(r2 != NULL) gzputs(gz2 ? gz2 : gz1, buf2->buff);

  // release mutex
  pthread_mutex_unlock(&output->lock);
}

// pthread method, loop: grabs job, does processing
static void* correct_reads_thread(void *ptr)
{
  CorrectReadsWorker *wrkr = (CorrectReadsWorker*)ptr;
  MsgPool *pool = wrkr->pool;
  AsyncIOData *data;
  int pos;

  while((pos = msgpool_claim_read(pool)) != -1)
  {
    memcpy(&data, msgpool_get_ptr(pool, pos), sizeof(AsyncIOData*));
    correct_read(wrkr, data);
    msgpool_release(pool, pos, MPOOL_EMPTY);
  }

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
  argi = correct_reads_parse(argc, argv, inputs, &num_inputs,
                             false, true, NULL, NULL);

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
    ctx_max_kmers = MAX2(ctx_max_kmers, gfiles[i].num_of_kmers);
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

  for(i = 0; i < num_inputs && corrected_output_open(&outputs[i], &inputs[i]); i++)
    inputs[i].ptr = &outputs[i];

  // Check if something went wrong
  if(i < num_inputs) {
    for(j = 0; j < i; j++)
      corrected_output_delete(&outputs[i]);
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
  AsyncIOData *data = malloc2(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  MsgPool pool;
  msgpool_alloc_yield(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*));
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  CorrectReadsWorker *wrkrs = malloc2(num_work_threads * sizeof(CorrectReadsWorker));

  for(i = 0; i < num_work_threads; i++)
    correct_reads_worker_alloc(&wrkrs[i], &pool, &db_graph);

  AsyncIOReadTask *asyncio_tasks = malloc2(num_inputs * sizeof(AsyncIOReadTask));
  correct_reads_input_to_asycio(asyncio_tasks, inputs, num_inputs);

  // Load input files num_io_threads at a time
  size_t num_io_threads = args->max_io_threads, n;

  for(i = 0; i < num_inputs; i += num_io_threads) {
    n = MIN2(num_inputs - i, num_io_threads);
    asyncio_run_threads(&pool, asyncio_tasks+i, n, correct_reads_thread,
                        wrkrs, num_work_threads, sizeof(CorrectReadsWorker));
  }

  for(i = 0; i < num_work_threads; i++)
    correct_reads_worker_dealloc(&wrkrs[i]);

  for(i = 0; i < num_inputs; i++) {
    if(inputs[i].file1 != NULL) seq_close(inputs[i].file1);
    if(inputs[i].file2 != NULL) seq_close(inputs[i].file2);
    corrected_output_close(&outputs[i]);
  }

  free(wrkrs);
  free(asyncio_tasks);
  msgpool_dealloc(&pool);

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  free(data);

  free(inputs);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);
  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
