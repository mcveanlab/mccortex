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
  CorrectAlnWorker corrector;
  // For filling in gaps
  GraphWalker wlk;
  RepeatWalker rptwlk;
  StrBuf rbuf1, rbuf2, qbuf;

  // Corrected alignment
  dBNodeBuffer nodebuf; Int32Buffer posbuf;
} CorrectReadsWorker;

static void correct_reads_worker_alloc(CorrectReadsWorker *wrkr,
                                       const dBGraph *db_graph)
{
  wrkr->db_graph = db_graph;
  correct_aln_worker_alloc(&wrkr->corrector, db_graph);
  graph_walker_alloc(&wrkr->wlk);
  rpt_walker_alloc(&wrkr->rptwlk, db_graph->ht.capacity, 22); // 4MB bloom
  strbuf_alloc(&wrkr->rbuf1, 1024); // read1
  strbuf_alloc(&wrkr->rbuf2, 1024); // read2
  strbuf_alloc(&wrkr->qbuf, 1024); // quality scores
  db_node_buf_alloc(&wrkr->nodebuf, 512);
  int32_buf_alloc(&wrkr->posbuf, 512);
}

static void correct_reads_worker_dealloc(CorrectReadsWorker *wrkr)
{
  correct_aln_worker_dealloc(&wrkr->corrector);
  graph_walker_dealloc(&wrkr->wlk);
  rpt_walker_dealloc(&wrkr->rptwlk);
  strbuf_dealloc(&wrkr->rbuf1);
  strbuf_dealloc(&wrkr->rbuf2);
  strbuf_dealloc(&wrkr->qbuf);
  db_node_buf_dealloc(&wrkr->nodebuf);
  int32_buf_dealloc(&wrkr->posbuf);
}

/*
// Don't need this now
static inline void _buf_append_qual_scores(StrBuf *buf, const read_t *r)
{
  // Print quality scores
  strbuf_append_str(buf, "+\n");
  if(r->qual.end) {
    strbuf_append_strn(buf, r->qual.b, r->qual.end);
  } else {
    strbuf_ensure_capacity(buf, buf->end + r->qual.end + 2);
    memset(buf->b+buf->end, '.', r->seq.end);
    buf->end += r->seq.end;
    buf->b[buf->end] = '\0';
  }
  strbuf_append_char(buf, '\n');
}
*/

static inline void strbuf_append_strn_lc(StrBuf *buf, const char *str, size_t len)
{
  strbuf_ensure_capacity(buf, buf->end + len);
  char *to = buf->b + buf->end;
  const char *end = str + len;
  for(; str < end; str++, to++) *to = tolower(*str);
  buf->end += len;
  buf->b[buf->end] = '\0';
}

static inline void strbuf_append_strn_uc(StrBuf *buf, const char *str, size_t len)
{
  strbuf_ensure_capacity(buf, buf->end + len);
  char *to = buf->b + buf->end;
  const char *end = str + len;
  for(; str < end; str++, to++) *to = toupper(*str);
  buf->end += len;
  buf->b[buf->end] = '\0';
}

// Append char `c` `n` times
static inline void strbuf_append_charn(StrBuf *buf, char c, size_t n)
{
  strbuf_ensure_capacity(buf, buf->end + n);
  memset(buf->b+buf->end, c, n);
  buf->end += n;
  buf->b[buf->end] = '\0';
}

/*
// Quality score may be too short, so this function is useful
static inline void _append_quals(StrBuf *buf, const read_t *r,
                                 size_t start, size_t len)
{
  size_t end = start + len, limit = MIN2(end, r->qual.end), n;
  strbuf_ensure_capacity(buf, buf->end + len);
  if(start < limit)
  {
    // Copy quality scores
    n = limit - start;
    memcpy(buf->b+buf->end, r->qual.b + start, n);
    start += n;
    buf->end += n;
  }
  if(start < end)
  {
    // Fill in missing quality scores with '.'
    n = end - start;
    memset(buf->b+buf->end, '.', n);
    buf->end += n;
  }
  buf->b[buf->end] = '\0';
}
*/

// Returns the new number of bases printed
static inline
size_t _print_read_kmer(const read_t *r, StrBuf *rbuf, StrBuf *qbuf,
                        size_t pos, size_t num_bases_printed,
                        const dBGraph *db_graph)
{
  size_t n, kmer_size = db_graph->kmer_size;

  // printf("first pos: %zu num_bases_printed: %zu\n", pos, num_bases_printed);

  if(pos > num_bases_printed) {
    // Fill in missing seq
    n = pos - num_bases_printed;
    // printf("  fill: %zu\n", n);
    strbuf_append_strn_lc(rbuf, r->seq.b+num_bases_printed,  n);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b+num_bases_printed, n);
    num_bases_printed = pos;
  }

  // Append bases that match a kmer
  if(pos + kmer_size > num_bases_printed) {
    n = pos + kmer_size - num_bases_printed;
    // printf("  mtch: %zu\n", n);
    strbuf_append_strn_uc(rbuf, r->seq.b  + num_bases_printed, n);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b + num_bases_printed, n);
  }

  return pos + kmer_size;
}

/*
// Get run of -1 values
static size_t _get_run_neg(const int32_t *arr, size_t len)
{
  size_t i = 0;
  while(i < len && arr[i] < 0) i++;
  return i;
}
*/

// Prints read sequence in lower case instead of N
static void handle_read3(CorrectReadsWorker *wrkr,
                         const CorrectAlnParam *params,
                         const read_t *r, StrBuf *rbuf, StrBuf *qbuf,
                         uint8_t fq_cutoff, uint8_t hp_cutoff,
                         dBNodeBuffer *nodebuf, Int32Buffer *posbuf)
{
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = db_graph->kmer_size;
  CorrectAlnWorker *corrector = &wrkr->corrector;

  correct_aln_read(corrector, params, r, fq_cutoff, hp_cutoff,
                   nodebuf, posbuf);

  ctx_assert(nodebuf->len == posbuf->len);

  if(nodebuf->len == 0)
  {
    // Didn't get any kmers
    // Copy whole sequence as lowercase
    strbuf_append_strn_lc(rbuf, r->seq.b,  r->seq.end);
    // Copy whole quality score
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b, r->qual.end);
    return;
  }

  // There must be at least one kmer where rpos != -1

  size_t i = 0, j, num_neg;
  size_t bases_printed = 0;
  Nucleotide nuc;

  const size_t num_nodes = nodebuf->len;
  const dBNode *node_arr = nodebuf->data;
  const int32_t *pos_arr = posbuf->data;

  // Deal with neg kmers at start here
  for(num_neg = 0; num_neg < num_nodes && pos_arr[num_neg] == -1; num_neg++) {}
  ctx_assert(num_neg < num_nodes);

  if(num_neg > 0)
  {
    // Append existing up until first "filled in" kmer
    strbuf_append_strn_lc(rbuf, r->seq.b, pos_arr[num_neg] - num_neg);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b, pos_arr[num_neg] - num_neg);

    strbuf_ensure_capacity(rbuf, rbuf->end + num_neg);
    strbuf_ensure_capacity(qbuf, qbuf->end + num_neg);

    // Copy first base from each kmer
    for(j = 0; j < num_neg; j++) {
      nuc = db_node_get_first_nuc(node_arr[j], db_graph);
      rbuf->b[rbuf->end++] = dna_nuc_to_char(nuc);
      qbuf->b[qbuf->end++] = '.';
    }
  }


  while(i < num_nodes)
  {
    j = i;
    while(j < num_nodes && pos_arr[j] < 0) j++;
    num_neg = j - i;

    if(num_neg == 0) {
      // Print non-neg
      bases_printed = _print_read_kmer(r, rbuf, qbuf, pos_arr[i],
                                       bases_printed, db_graph);
      i++; // Go to next node
    }
    else if(i + num_neg == num_nodes) break;
    else {
      // size_t end = (size_t)pos_arr[i + num_neg];
      // int32_t exp_kmers = pos_arr[i+num_neg] + kmer_size - bases_printed - 2;

      if(num_neg >= kmer_size) {
        // Deal with run of missing kmers
        size_t end = i + num_neg - kmer_size;

        for(; i < end; i++) {
          // Grab last base
          nuc = db_node_get_last_nuc(node_arr[i+j], db_graph);
          strbuf_append_char(rbuf, dna_nuc_to_char(nuc)); // prints uppercase
          strbuf_append_char(qbuf, '.');
        }

        // bases_printed = pos_arr[i + num_neg];
        // bases_printed = pos_arr[i + num_neg] + exp_kmers - num_neg;
      }
      else {
        // bases_printed = pos_arr[i + num_neg] + exp_kmers - num_neg;
      }

      size_t nextpos = pos_arr[i + num_neg];
      if(num_neg < kmer_size) nextpos += kmer_size - num_neg;

      bases_printed = MAX2(bases_printed, nextpos);

      // bases_printed = pos_arr[i + num_neg] + kmer_size -
      //                 (num_neg >= kmer_size ? kmer_size : num_neg) - 1;
        // bases_printed = MAX2(bases_printed + num_neg - kmer_size - 1,
        //                      (size_t)pos_arr[i + num_neg] + exp_kmers - num_neg);

      i += num_neg;
    }
  }

  // Deal with filled in right hand side
  if(i < num_nodes)
  {
    ctx_assert(num_neg > 0);

    strbuf_ensure_capacity(rbuf, rbuf->end + num_neg);
    strbuf_ensure_capacity(qbuf, qbuf->end + num_neg);

    // Copy first base from each kmer
    for(j = 0; j < num_neg; j++) {
      nuc = db_node_get_last_nuc(node_arr[j], db_graph);
      rbuf->b[rbuf->end++] = dna_nuc_to_char(nuc);
      qbuf->b[qbuf->end++] = '.';
    }

    bases_printed += num_neg;
    size_t rem_bases = r->seq.end - bases_printed;
    strbuf_append_strn_lc(rbuf, r->seq.b  + bases_printed, rem_bases);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b+bases_printed, rem_bases);
  }

  if(r->qual.end == 0) strbuf_reset(qbuf);
}

// Print read in FASTA, FASTQ or PLAIN format
static void handle_read2(CorrectReadsWorker *wrkr,
                         const CorrectAlnParam *params,
                         const read_t *r, StrBuf *rbuf, StrBuf *qbuf,
                         uint8_t fq_cutoff, uint8_t hp_cutoff,
                         dBNodeBuffer *nodebuf, Int32Buffer *posbuf,
                         seq_format format)
{
  strbuf_reset(rbuf);
  strbuf_reset(qbuf); // quality scores go here

  if((format & SEQ_FMT_FASTQ) && r->qual.end && r->seq.end != r->qual.end) {
    die("Read qual scores don't match seq length (%zu vs %zu): %s",
        r->qual.end, r->seq.end, r->name.b);
  }

  if(format & (SEQ_FMT_FASTA | SEQ_FMT_FASTQ)) {
    strbuf_append_char(rbuf, format == SEQ_FMT_FASTA ? '>' : '@');
    strbuf_append_strn(rbuf, r->name.b, r->name.end);
    strbuf_append_char(rbuf, '\n');
  }

  // Write read sequence string to rbuf, quality scores string to qbuf
  size_t orig_len = rbuf->end;
  handle_read3(wrkr, params, r, rbuf, qbuf, fq_cutoff, hp_cutoff, nodebuf, posbuf);
  size_t readlen = rbuf->end - orig_len;

  // Copy quality scores to read buffer ready to print
  if(format & SEQ_FMT_FASTQ) {
    strbuf_append_str(rbuf, "\n+\n");
    if(r->qual.end == 0) strbuf_append_charn(rbuf, '.', readlen);
    else                 strbuf_append_strn(rbuf, qbuf->b, qbuf->end);
    strbuf_append_char(rbuf, '\n');
  }

  strbuf_reset(qbuf);
}

/*
// Prints read sequence in lower case instead of N
static void handle_read(CorrectReadsWorker *wrkr,
                        const CorrectAlnParam *params,
                        const read_t *r, StrBuf *buf,
                        uint8_t fq_cutoff, uint8_t hp_cutoff,
                        seq_format format)
{
  dBNodeBuffer *nbuf, *tmpnbuf = &wrkr->tmpnbuf;
  GraphWalker *wlk = &wrkr->wlk;
  RepeatWalker *rptwlk = &wrkr->rptwlk;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = db_graph->kmer_size;
  const size_t ctxcol = params->ctxcol;
  const size_t ctpcol = params->ctpcol;

  size_t i, idx, gap, num_n, nbases;
  size_t init_len, end_len;
  BinaryKmer bkmer;
  Nucleotide nuc;
  char bkmerstr[MAX_KMER_SIZE+1];

  // Print read in FASTA format
  strbuf_reset(buf);

  if(format & (SEQ_FMT_FASTA | SEQ_FMT_FASTQ)) {
    strbuf_append_char(buf, format == SEQ_FMT_FASTA ? '>' : '@');
    strbuf_append_strn(buf, r->name.b, r->name.end);
    strbuf_append_char(buf, '\n');
  }

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, params,
                         r, NULL, fq_cutoff, 0, hp_cutoff);

  // Get first alignment
  nbuf = correct_alignment_nxt(&wrkr->corrector);

  if(nbuf == NULL) {
    // Alignment failed - copy read in lower case
    size_t offset = buf->end;
    strbuf_append_strn(buf, r->seq.b, r->seq.end);
    for(i = offset; i < buf->end; i++) buf->b[i] = tolower(buf->b[i]);
    strbuf_append_char(buf, '\n');
    if(format == SEQ_FMT_FASTQ) _buf_append_qual_scores(buf, r);
    return;
  }

  // extend left
  dBAlignment *aln = &wrkr->corrector.aln;
  size_t left_gap = aln->rpos.data[0], right_gap = aln->r1enderr;
  size_t bases_printed = 0;

  if(left_gap > 0)
  {
    // Walk left
    graph_walker_prime(wlk, nbuf->data, nbuf->len,
                       params->max_context, false,
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
    ctx_assert(idx > 0);
    gap = aln->rpos.data[idx] - aln->rpos.data[idx-1] - 1;
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
                       params->max_context, true,
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
  if(format == SEQ_FMT_FASTQ) _buf_append_qual_scores(buf, r);
}
*/

static void correct_read(CorrectReadsWorker *wrkr, AsyncIOData *data)
{
  uint8_t fq_cutoff1, fq_cutoff2, hp_cutoff;

  CorrectAlnInput *input = (CorrectAlnInput*)data->ptr;
  const CorrectAlnParam *params = &input->crt_params;
  SeqOutput *output = input->output;
  StrBuf *rbuf1 = &wrkr->rbuf1, *rbuf2 = &wrkr->rbuf2, *qbuf = &wrkr->qbuf;
  dBNodeBuffer *nodebuf = &wrkr->nodebuf;
  Int32Buffer *posbuf = &wrkr->posbuf;
  seq_format format = output->fmt;

  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  fq_cutoff1 = fq_cutoff2 = input->fq_cutoff;

  if(fq_cutoff1 > 0) {
    fq_cutoff1 += data->fq_offset1;
    fq_cutoff2 += data->fq_offset2;
  }

  hp_cutoff = input->hp_cutoff;

  strbuf_reset(rbuf1);
  strbuf_reset(rbuf2);

  if(r2 == NULL)
  {
    // Single ended read
    handle_read2(wrkr, params, r1, rbuf1, qbuf, fq_cutoff1, hp_cutoff,
                 nodebuf, posbuf, format);
    pthread_mutex_lock(&output->lock_se);
    gzwrite(output->gzout_se, rbuf1->b, rbuf1->end);
    pthread_mutex_unlock(&output->lock_se);
  }
  else
  {
    // Paired-end reads
    handle_read2(wrkr, params, r1, rbuf1, qbuf, fq_cutoff1, hp_cutoff,
                 nodebuf, posbuf, format);
    handle_read2(wrkr, params, r2, rbuf2, qbuf, fq_cutoff2, hp_cutoff,
                 nodebuf, posbuf, format);
    pthread_mutex_lock(&output->lock_pe);
    gzwrite(output->gzout_pe[0], rbuf1->b, rbuf1->end);
    gzwrite(output->gzout_pe[1], rbuf2->b, rbuf2->end);
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
