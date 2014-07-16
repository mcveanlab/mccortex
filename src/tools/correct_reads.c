#include "global.h"
#include "correct_reads.h"
#include "correct_alignment.h"
#include "async_read_io.h"
#include "loading_stats.h"
#include "seq_reader.h"
#include "file_util.h"
#include "msg-pool/msgpool.h"

typedef struct
{
  const dBGraph *db_graph;
  volatile size_t *rcounter;
  dBAlignment aln;
  CorrectAlnWorker corrector;
  // For filling in gaps
  GraphWalker wlk;
  RepeatWalker rptwlk;
  StrBuf buf1, buf2;
  dBNodeBuffer tmpnbuf;
} CorrectReadsWorker;

static void correct_reads_worker_alloc(CorrectReadsWorker *wrkr,
                                       const dBGraph *db_graph)
{
  wrkr->db_graph = db_graph;
  db_alignment_alloc(&wrkr->aln);
  correct_aln_worker_alloc(&wrkr->corrector, db_graph);
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

// Prints read sequence in lower case instead of N
static void handle_read(CorrectReadsWorker *wrkr,
                        const CorrectAlnInput *input,
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

  // Print read in FASTA format
  strbuf_reset(buf);
  strbuf_append_char(buf, '>');
  strbuf_append_strn(buf, r->name.b, r->name.end);
  strbuf_append_char(buf, '\n');

  // Get de Bruijn graph alignment
  db_alignment_from_reads(&wrkr->aln, r, NULL,
                          fq_cutoff, 0, hp_cutoff, db_graph, -1);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, input->crt_params);

  // Get first alignment
  nbuf = correct_alignment_nxt(&wrkr->corrector);

  if(nbuf == NULL) {
    // Alignment failed - copy read in lower case
    size_t offset = buf->end;
    strbuf_append_strn(buf, r->seq.b, r->seq.end);
    for(i = offset; i < buf->end; i++) buf->b[i] = tolower(buf->b[i]);
    strbuf_append_char(buf, '\n');
    return;
  }

  // extend left
  size_t left_gap = aln->gaps.data[0], right_gap = aln->r1enderr;
  size_t bases_printed = 0;

  if(left_gap > 0)
  {
    // Walk left
    graph_walker_prime(wlk, nbuf->data, nbuf->len,
                       input->crt_params.max_context, false,
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
    size_t unresolved = left_gap - tmpnbuf->len;
    for(; bases_printed < unresolved; bases_printed++) {
      strbuf_append_char(buf, tolower(r->seq.b[bases_printed])); // N
    }

    // Append bases
    for(i = tmpnbuf->len-1; i != SIZE_MAX; i--) {
      nuc = db_node_get_first_nuc(db_node_reverse(tmpnbuf->data[i]), db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }

    bases_printed += tmpnbuf->len;
  }

  // Append first contig
  bkmer = db_node_oriented_bkmer(db_graph, nbuf->data[0]);
  binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
  strbuf_append_strn(buf, bkmerstr, kmer_size);
  for(i = 1; i < nbuf->len; i++) {
    nuc = db_node_get_last_nuc(nbuf->data[i], db_graph);
    strbuf_append_char(buf, dna_nuc_to_char(nuc));
  }
  bases_printed += kmer_size + nbuf->len - 1;

  while(correct_alignment_get_endidx(&wrkr->corrector) < aln->nodes.len)
  {
    nbuf = correct_alignment_nxt(&wrkr->corrector);
    ctx_assert(nbuf != NULL);
    idx = correct_alignment_get_strtidx(&wrkr->corrector);
    gap = aln->gaps.data[idx];
    num_n = gap < kmer_size ? 0 : gap - kmer_size + 1;

    for(i = 0; i < num_n; i++) {
      strbuf_append_char(buf, tolower(r->seq.b[bases_printed++])); // N
    }

    nbases = MIN2(gap+1, kmer_size);
    bkmer = db_node_oriented_bkmer(db_graph, nbuf->data[0]);
    binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
    strbuf_append_strn(buf, bkmerstr+kmer_size-nbases, nbases);

    for(i = 1; i < nbuf->len; i++) {
      nuc = db_node_get_last_nuc(nbuf->data[i], db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }

    bases_printed += nbases + nbuf->len - 1;
  }

  // extend right
  if(right_gap > 0)
  {
    // walk right
    graph_walker_prime(wlk, nbuf->data, nbuf->len,
                       input->crt_params.max_context, true,
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
      nuc = db_node_get_last_nuc(nbuf->data[i], db_graph);
      strbuf_append_char(buf, dna_nuc_to_char(nuc));
    }
    bases_printed += nbuf->len - init_len;

    // Add Ns for bases we couldn't resolve
    for(; bases_printed < r->seq.end; bases_printed++) {
      strbuf_append_char(buf, tolower(r->seq.b[bases_printed])); // N
    }
  }

  strbuf_append_char(buf, '\n');
}

static void correct_read(CorrectReadsWorker *wrkr, AsyncIOData *data)
{
  uint8_t fq_cutoff1, fq_cutoff2, hp_cutoff;

  CorrectAlnInput *input = (CorrectAlnInput*)data->ptr;
  SeqOutput *output = input->output;
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
    pthread_mutex_lock(&output->lock_se);
    gzwrite(output->gzout_se, buf1->b, buf1->end);
    pthread_mutex_unlock(&output->lock_se);
  }
  else
  {
    // Paired-end reads
    handle_read(wrkr, input, r1, buf1, fq_cutoff1, hp_cutoff);
    handle_read(wrkr, input, r2, buf2, fq_cutoff2, hp_cutoff);
    pthread_mutex_lock(&output->lock_pe);
    gzwrite(output->gzout_pe[0], buf1->b, buf1->end);
    gzwrite(output->gzout_pe[1], buf2->b, buf2->end);
    pthread_mutex_unlock(&output->lock_pe);
  }
}

// pthread method, loop: grabs job, does processing
static void correct_reads_thread(AsyncIOData *data, void *ptr)
{
  CorrectReadsWorker *wrkr = (CorrectReadsWorker*)ptr;
  correct_read(wrkr, data);

  // Print progress
  size_t n = __sync_add_and_fetch(wrkr->rcounter, 1);
  ctx_update("CorrectReads", n);
}

// Correct reads against the graph, and print out
void correct_reads(CorrectAlnInput *inputs, size_t num_inputs,
                   const char *dump_seqgap_hist_path,
                   const char *dump_fraglen_hist_path,
                   size_t num_threads, const dBGraph *db_graph)
{
  size_t i, n, read_counter = 0;

  CorrectReadsWorker *wrkrs = ctx_calloc(num_threads, sizeof(CorrectReadsWorker));

  for(i = 0; i < num_threads; i++) {
    correct_reads_worker_alloc(&wrkrs[i], db_graph);
    wrkrs[i].rcounter = &read_counter;
  }

  AsyncIOInput *asyncio_tasks = ctx_calloc(num_inputs, sizeof(AsyncIOInput));
  correct_aln_input_to_asycio(asyncio_tasks, inputs, num_inputs);

  // Load input files MAX_IO_THREADS at a time
  for(i = 0; i < num_inputs; i += MAX_IO_THREADS) {
    n = MIN2(num_inputs - i, MAX_IO_THREADS);
    asyncio_run_pool(asyncio_tasks+i, n, correct_reads_thread,
                     wrkrs, num_threads, sizeof(CorrectReadsWorker));
  }

  // Merge stats into workers[0]
  for(i = 1; i < num_threads; i++)
    correct_aln_merge_stats(&wrkrs[0].corrector, &wrkrs[i].corrector);

  LoadingStats *stats = &wrkrs[0].corrector.stats;
  CorrectAlnStats *gapstats = &wrkrs[0].corrector.gapstats;

  correct_aln_dump_stats(stats, gapstats,
                         dump_seqgap_hist_path,
                         dump_fraglen_hist_path,
                         db_graph->ht.num_kmers);

  for(i = 0; i < num_threads; i++)
    correct_reads_worker_dealloc(&wrkrs[i]);

  ctx_free(wrkrs);
  ctx_free(asyncio_tasks);
}
