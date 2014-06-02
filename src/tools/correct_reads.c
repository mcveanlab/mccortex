#include "global.h"
#include "correct_reads.h"
#include "correct_alignment.h"
#include "async_read_io.h"
#include "loading_stats.h"
#include "seq_reader.h"
#include "file_util.h"
#include "msg-pool/msgpool.h"

#define MAX_CONTEXT 50

typedef struct
{
  const dBGraph *db_graph;
  MsgPool *pool;
  dBAlignment aln;
  CorrectAlnWorker corrector;
  LoadingStats stats;
  // For filling in gaps
  GraphWalker wlk;
  RepeatWalker rptwlk;
  StrBuf buf1, buf2;
  dBNodeBuffer tmpnbuf;
} CorrectReadsWorker;

// Returns true on success, false on error
bool corrected_output_open(CorrectedOutput *out, const CorrectAlnTask *in)
{
  out->gzse = out->gzpe1 = out->gzpe2 = NULL;
  const char *base = (const char*)in->ptr;
  StrBuf *outse = &out->outse, *outpe1 = &out->outpe1, *outpe2 = &out->outpe2;
  strbuf_alloc(outse, 512);
  strbuf_alloc(outpe1, 512);
  strbuf_alloc(outpe2, 512);

  if(pthread_mutex_init(&out->lockse, NULL) != 0) die("Mutex init failed");
  if(pthread_mutex_init(&out->lockpe, NULL) != 0) die("Mutex init failed");

  if(in->files.interleaved || in->files.file2 == NULL)
  {
    strbuf_append_str(outse, base);
    strbuf_append_str(outse, ".fa.gz");
    if(futil_file_exists(outse->buff)) {
      warn("SE File already exists: %s", outse->buff);
      return false;
    }
    out->gzse = gzopen(outse->buff, "w");
    if(out->gzse == NULL) {
      warn("Cannot write to: %s", outse->buff);
      return false;
    }
  }

  if(in->files.interleaved || in->files.file2 != NULL)
  {
    strbuf_append_str(outpe1, base);
    strbuf_append_str(outpe2, base);
    strbuf_append_str(outpe1, ".1.fa.gz");
    strbuf_append_str(outpe2, ".2.fa.gz");

    bool fexists = false;
    if(futil_file_exists(outpe1->buff)) {
      warn("File already exists: %s", outpe1->buff);
      fexists = true;
    }
    if(futil_file_exists(outpe2->buff)) {
      warn("File already exists: %s", outpe2->buff);
      fexists = true;
    }
    if(fexists) return false;

    out->gzpe1 = gzopen(outpe1->buff, "w");
    out->gzpe2 = gzopen(outpe2->buff, "w");
    if(out->gzpe1 == NULL) warn("Cannot write to: %s", outpe1->buff);
    if(out->gzpe2 == NULL) warn("Cannot write to: %s", outpe2->buff);
    if(out->gzpe1 == NULL || out->gzpe2 == NULL) return false;
  }

  return true;
}

void corrected_output_close(CorrectedOutput *out)
{
  if(out->gzse != NULL) gzclose(out->gzse);
  if(out->gzpe1 != NULL) gzclose(out->gzpe1);
  if(out->gzpe2 != NULL) gzclose(out->gzpe2);
  out->gzse = out->gzpe1 = out->gzpe2 = NULL;
  strbuf_dealloc(&out->outse);
  strbuf_dealloc(&out->outpe1);
  strbuf_dealloc(&out->outpe2);
  pthread_mutex_destroy(&out->lockse);
  pthread_mutex_destroy(&out->lockpe);
}

void corrected_output_delete(CorrectedOutput *out)
{
  if(out->gzse != NULL) unlink(out->outse.buff);
  if(out->gzpe1 != NULL) unlink(out->outpe1.buff);
  if(out->gzpe2 != NULL) unlink(out->outpe2.buff);
  corrected_output_close(out);
}

static void correct_reads_worker_alloc(CorrectReadsWorker *wrkr,
                                       MsgPool *pool,
                                       const dBGraph *db_graph)
{
  wrkr->db_graph = db_graph;
  wrkr->pool = pool;
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

static void handle_read(CorrectReadsWorker *wrkr,
                        const CorrectAlnTask *input,
                        const read_t *r, StrBuf *buf,
                        uint8_t fq_cutoff, uint8_t hp_cutoff)
{
  dBNodeBuffer *nbuf, *tmpnbuf = &wrkr->tmpnbuf;
  dBAlignment *aln = &wrkr->aln;
  GraphWalker *wlk = &wrkr->wlk;
  RepeatWalker *rptwlk = &wrkr->rptwlk;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = db_graph->kmer_size;
  const size_t ctxcol = input->crt_params.ctxcol;
  const size_t ctpcol = input->crt_params.ctpcol;

  size_t i, idx, gap, num_n, nbases;
  size_t init_len, end_len;
  BinaryKmer bkmer;
  Nucleotide nuc;
  char bkmerstr[MAX_KMER_SIZE+1];

  db_alignment_from_reads(&wrkr->aln, r, NULL,
                          fq_cutoff, 0, hp_cutoff, db_graph, -1);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, input->crt_params);

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
                       ctxcol, ctpcol, db_graph);

    db_node_buf_reset(tmpnbuf);
    db_node_buf_ensure_capacity(nbuf, left_gap);

    while(tmpnbuf->len < left_gap && graph_walker_next(wlk) &&
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
  bkmer = db_node_oriented_bkmer(db_graph, nbuf->data[0]);
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

    while(nbuf->len < end_len && graph_walker_next(wlk) &&
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

  CorrectAlnTask *input = (CorrectAlnTask*)data->ptr;
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

  if(r2 == NULL)
  {
    // Single ended read
    handle_read(wrkr, input, r1, buf1, fq_cutoff1, hp_cutoff);
    pthread_mutex_lock(&output->lockse);
    gzputs(output->gzse, buf1->buff);
    pthread_mutex_unlock(&output->lockse);

    // Update stats
    wrkr->stats.num_se_reads++;
  }
  else
  {
    // Paired-end reads
    handle_read(wrkr, input, r1, buf1, fq_cutoff1, hp_cutoff);
    handle_read(wrkr, input, r2, buf2, fq_cutoff2, hp_cutoff);
    pthread_mutex_lock(&output->lockpe);
    gzputs(output->gzpe1, buf1->buff);
    gzputs(output->gzpe2, buf2->buff);
    pthread_mutex_unlock(&output->lockpe);

    // Update stats
    wrkr->stats.num_pe_reads += 2;
  }
}

// pthread method, loop: grabs job, does processing
static void correct_reads_thread(void *ptr)
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
}

// `input` and `outputs` should both be of length `num_inputs`
void correct_reads(size_t num_threads, size_t max_io_threads,
                   CorrectAlnTask *inputs, CorrectedOutput *outputs,
                   size_t num_inputs, const dBGraph *db_graph)
{
  size_t i, n;
  AsyncIOData *data = ctx_malloc(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  // Create pool of AsyncIOData* that point to elements in the above array
  // -> swapping of pointers faster than whole AsyncIOData elements
  MsgPool pool;
  msgpool_alloc(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*), USE_MSG_POOL);
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  CorrectReadsWorker *wrkrs = ctx_calloc(num_threads, sizeof(CorrectReadsWorker));

  for(i = 0; i < num_threads; i++)
    correct_reads_worker_alloc(&wrkrs[i], &pool, db_graph);

  AsyncIOReadTask *asyncio_tasks = ctx_calloc(num_inputs, sizeof(AsyncIOReadTask));
  correct_reads_input_to_asycio(asyncio_tasks, inputs, num_inputs);

  // Load input files max_io_threads at a time
  for(i = 0; i < num_inputs; i += max_io_threads) {
    n = MIN2(num_inputs - i, max_io_threads);
    asyncio_run_threads(&pool, asyncio_tasks+i, n, correct_reads_thread,
                        wrkrs, num_threads, sizeof(CorrectReadsWorker));
  }

  for(i = 0; i < num_threads; i++)
    correct_reads_worker_dealloc(&wrkrs[i]);

  for(i = 0; i < num_inputs; i++) {
    asyncio_task_close(&inputs[i].files);
    corrected_output_close(&outputs[i]);
  }

  ctx_free(wrkrs);
  ctx_free(asyncio_tasks);
  msgpool_dealloc(&pool);

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  ctx_free(data);
}
