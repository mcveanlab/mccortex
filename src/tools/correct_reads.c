#include "global.h"
#include "correct_reads.h"
#include "correct_alignment.h"
#include "async_read_io.h"
#include "seq_loading_stats.h"
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
  char fq_zero; // character to use to zero fastq [default: '.']
  bool append_orig_seq; // append sequence to name ">name prev=OLDSEQ"

  // Corrected alignment
  dBNodeBuffer nodebuf; Int32Buffer posbuf;
} CorrectReadsWorker;

static void correct_reads_worker_alloc(CorrectReadsWorker *wrkr,
                                       size_t *read_cntr_ptr,
                                       bool append_orig_seq, char fq_zero,
                                       const dBGraph *db_graph)
{
  wrkr->rcounter = read_cntr_ptr;
  wrkr->fq_zero = fq_zero;
  wrkr->append_orig_seq = append_orig_seq;
  wrkr->db_graph = db_graph;
  correct_aln_worker_alloc(&wrkr->corrector, false, db_graph);
  graph_walker_alloc(&wrkr->wlk, db_graph);
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
    // printf("  fill: %zu '%.*s'\n", n, (int)n, r->seq.b + num_bases_printed);
    strbuf_append_strn_lc(rbuf, r->seq.b+num_bases_printed,  n);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b+num_bases_printed, n);
    num_bases_printed = pos;
  }

  // Append bases that match a kmer
  if(pos + kmer_size > num_bases_printed) {
    n = pos + kmer_size - num_bases_printed;
    // printf("  mtch: %zu '%.*s'\n", n, (int)n, r->seq.b + num_bases_printed);
    strbuf_append_strn_uc(rbuf, r->seq.b + num_bases_printed, n);
    if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b + num_bases_printed, n);
  }

  return pos + kmer_size;
}

// Prints read sequence in lower case instead of N
static void handle_read2(CorrectReadsWorker *wrkr,
                         const CorrectAlnParam *params,
                         const read_t *r, StrBuf *rbuf, StrBuf *qbuf,
                         uint8_t fq_cutoff, uint8_t hp_cutoff,
                         dBNodeBuffer *nodebuf, Int32Buffer *posbuf)
{
  const char fq_zero = wrkr->fq_zero;
  const dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = db_graph->kmer_size;
  CorrectAlnWorker *corrector = &wrkr->corrector;

  correct_aln_read(corrector, params, r, fq_cutoff, hp_cutoff,
                   nodebuf, posbuf);

  // db_alignment_print(&corrector->aln);

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
  const dBNode *node_arr = nodebuf->b;
  const int32_t *pos_arr = posbuf->b;

  // Deal with neg kmers at start here
  for(num_neg = 0; num_neg < num_nodes && pos_arr[num_neg] == -1; num_neg++) {}
  ctx_assert(num_neg < num_nodes);

  // Append existing up until first "filled in" kmer
  strbuf_append_strn_lc(rbuf, r->seq.b, pos_arr[num_neg] - num_neg);
  if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b, pos_arr[num_neg] - num_neg);

  if(num_neg > 0)
  {
    strbuf_ensure_capacity(rbuf, rbuf->end + num_neg);
    strbuf_ensure_capacity(qbuf, qbuf->end + num_neg);

    // Copy first base from each kmer
    for(j = 0; j < num_neg; j++) {
      // printf("%zu: %zu\n", j, node_arr[j].key);
      ctx_assert(HASH_ENTRY_ASSIGNED(db_graph->ht.table[node_arr[j].key]));
      nuc = db_node_get_first_nuc(node_arr[j], db_graph);
      rbuf->b[rbuf->end++] = dna_nuc_to_char(nuc);
      qbuf->b[qbuf->end++] = fq_zero;
    }

    rbuf->b[rbuf->end] = qbuf->b[qbuf->end] = '\0';
  }

  i = num_neg;
  bases_printed = pos_arr[num_neg];

  while(i < num_nodes)
  {
    // get num_neg -- the number of kmers filling a gap
    for(j = i; j < num_nodes && pos_arr[j] < 0; j++) {}
    num_neg = j - i;

    if(num_neg == 0) {
      // Print non-neg
      // printf(" pos_arr[i:%zu]: %i\n", i, pos_arr[i]);
      bases_printed = _print_read_kmer(r, rbuf, qbuf, pos_arr[i],
                                       bases_printed, db_graph);
      i++; // Go to next node
    }
    else if(i + num_neg == num_nodes) break;
    else
    {
      size_t exp_kmers = pos_arr[i+num_neg] - pos_arr[i-1] - 1;
      ctx_assert(pos_arr[i+num_neg] >= 0);

      // status("num_neg: %zu exp_kmers: %zu kmer_size: %zu [%s]",
      //        num_neg, exp_kmers, kmer_size, r->name.b);

      size_t nprint = 0;
      if(num_neg >= kmer_size) nprint = num_neg - kmer_size + 1;
      if(num_neg > exp_kmers) nprint = num_neg - exp_kmers;

      if(nprint) {
        // Deal with run of missing kmers
        size_t end = i + nprint;
        // status("print: %zu", nprint);
        for(j = i; j < end; j++) {
          // Grab last base
          nuc = db_node_get_last_nuc(node_arr[j], db_graph);
          strbuf_append_char(rbuf, dna_nuc_to_char(nuc)); // prints uppercase
          strbuf_append_char(qbuf, fq_zero);
        }
      }

      size_t nextpos = pos_arr[i + num_neg];
      if(num_neg < kmer_size) nextpos += kmer_size - num_neg - 1;
      bases_printed = MAX2(bases_printed, nextpos);

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
    for(j = i; j < num_nodes; j++) {
      nuc = db_node_get_last_nuc(node_arr[j], db_graph);
      rbuf->b[rbuf->end++] = dna_nuc_to_char(nuc);
      qbuf->b[qbuf->end++] = fq_zero;
    }

    ctx_assert(bases_printed == pos_arr[i-1]+kmer_size);
    bases_printed += num_neg;
  }

  size_t rem_bases = r->seq.end - bases_printed;
  strbuf_append_strn_lc(rbuf, r->seq.b + bases_printed, rem_bases);
  if(r->qual.end > 0) strbuf_append_strn(qbuf, r->qual.b+bases_printed, rem_bases);
  bases_printed += rem_bases;

  ctx_assert2(bases_printed == r->seq.end, "%zu vs %zu", bases_printed, r->seq.end);

  if(r->qual.end == 0) strbuf_reset(qbuf);
}

// Print read in FASTA, FASTQ or PLAIN format
static void handle_read(CorrectReadsWorker *wrkr,
                        const CorrectAlnParam *params,
                        const read_t *r, StrBuf *rbuf, StrBuf *qbuf,
                        uint8_t fq_cutoff, uint8_t hp_cutoff,
                        dBNodeBuffer *nodebuf, Int32Buffer *posbuf,
                        seq_format format, bool append_orig_seq)
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
    if(append_orig_seq) {
      strbuf_append_str(rbuf, " orig=");
      strbuf_append_strn(rbuf, r->seq.b, r->seq.end);
    }
    strbuf_append_char(rbuf, '\n');
  }

  // Write read sequence string to rbuf, quality scores string to qbuf
  size_t orig_len = rbuf->end;
  handle_read2(wrkr, params, r, rbuf, qbuf, fq_cutoff, hp_cutoff, nodebuf, posbuf);
  size_t readlen = rbuf->end - orig_len;

  // Copy quality scores to read buffer ready to print
  if(format & SEQ_FMT_FASTQ) {
    strbuf_append_str(rbuf, "\n+\n");
    if(r->qual.end == 0) strbuf_append_charn(rbuf, wrkr->fq_zero, readlen);
    else                 strbuf_append_strn(rbuf, qbuf->b, qbuf->end);
  }

  strbuf_append_char(rbuf, '\n');
  strbuf_reset(qbuf);
}


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
    handle_read(wrkr, params, r1, rbuf1, qbuf, fq_cutoff1, hp_cutoff,
                nodebuf, posbuf, format, wrkr->append_orig_seq);
    pthread_mutex_lock(&output->lock_se);
    gzwrite(output->gzout_se, rbuf1->b, rbuf1->end);
    pthread_mutex_unlock(&output->lock_se);
  }
  else
  {
    // Paired-end reads
    handle_read(wrkr, params, r1, rbuf1, qbuf, fq_cutoff1, hp_cutoff,
                nodebuf, posbuf, format, wrkr->append_orig_seq);
    handle_read(wrkr, params, r2, rbuf2, qbuf, fq_cutoff2, hp_cutoff,
                nodebuf, posbuf, format, wrkr->append_orig_seq);
    pthread_mutex_lock(&output->lock_pe);
    gzwrite(output->gzout_pe[0], rbuf1->b, rbuf1->end);
    gzwrite(output->gzout_pe[1], rbuf2->b, rbuf2->end);
    pthread_mutex_unlock(&output->lock_pe);
  }
}

// pthread method, loop: grabs job, does processing
static void correct_reads_thread(AsyncIOData *data, size_t threadid, void *ptr)
{
  (void)threadid;
  CorrectReadsWorker *wrkr = (CorrectReadsWorker*)ptr;
  correct_read(wrkr, data);

  // Print progress
  size_t n = __sync_add_and_fetch(wrkr->rcounter, 1);
  ctx_update("CorrectReads", n);
}

// Correct reads against the graph, and print out
// @param fq_zero use to fill quality scores; defaults to '.' if zero
void correct_reads(CorrectAlnInput *inputs, size_t num_inputs,
                   const char *dump_seqgap_hist_path,
                   const char *dump_fraglen_hist_path,
                   char fq_zero, bool append_orig_seq,
                   size_t num_threads, const dBGraph *db_graph)
{
  size_t i, n, read_counter = 0;

  if(!fq_zero) fq_zero = '.';

  CorrectReadsWorker *wrkrs = ctx_calloc(num_threads, sizeof(CorrectReadsWorker));

  for(i = 0; i < num_threads; i++) {
    correct_reads_worker_alloc(&wrkrs[i], &read_counter,
                               fq_zero, append_orig_seq,
                               db_graph);
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

  SeqLoadingStats *load_stats = &wrkrs[0].corrector.load_stats;
  CorrectAlnStats *aln_stats = &wrkrs[0].corrector.aln_stats;

  correct_aln_dump_stats(aln_stats, load_stats,
                         dump_seqgap_hist_path,
                         dump_fraglen_hist_path,
                         db_graph->ht.num_kmers);

  for(i = 0; i < num_threads; i++)
    correct_reads_worker_dealloc(&wrkrs[i]);

  ctx_free(wrkrs);
  ctx_free(asyncio_tasks);
}
