#include "global.h"
#include "db_alignment.h"
#include "seq_reader.h"


#define INIT_BUFLEN 1024

// Estimate initial memory required
size_t db_alignment_est_mem()
{
  return (sizeof(size_t)+sizeof(dBNode))*INIT_BUFLEN;
}

void db_alignment_alloc(dBAlignment *aln)
{
  db_node_buf_alloc(&aln->nodes, INIT_BUFLEN);
  int32_buf_alloc(&aln->rpos, INIT_BUFLEN);
}

void db_alignment_dealloc(dBAlignment *aln)
{
  db_node_buf_dealloc(&aln->nodes);
  int32_buf_dealloc(&aln->rpos);
  memset(aln, 0, sizeof(dBAlignment));
}

// if colour is -1 aligns to all colours, otherwise aligns to given colour only
// Returns number of kmers lost from the end
static size_t db_alignment_from_read(dBAlignment *aln, const read_t *r,
                                     uint8_t qcutoff, uint8_t hp_cutoff,
                                     const dBGraph *db_graph, int colour)
{
  size_t contig_start, contig_end = 0, search_start = 0;
  const size_t kmer_size = db_graph->kmer_size;

  BinaryKmer bkmer, tmp_key;
  Nucleotide nuc;
  hkey_t node;
  size_t i, offset, nxtbse;

  dBNodeBuffer *nodes = &aln->nodes;
  Int32Buffer *rpos = &aln->rpos;

  ctx_assert(nodes->len == rpos->len);
  size_t n = nodes->len, init_len = n;

  db_node_buf_capacity(nodes, n + r->seq.end);
  int32_buf_capacity(rpos, n + r->seq.end);

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qcutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qcutoff, hp_cutoff, &search_start);

    const char *contig = r->seq.b + contig_start;
    size_t contig_len = contig_end - contig_start;

    bkmer = binary_kmer_from_str(contig, kmer_size);
    bkmer = binary_kmer_right_shift_one_base(bkmer);

    for(offset=contig_start, nxtbse=kmer_size-1; nxtbse < contig_len; nxtbse++,offset++)
    {
      nuc = dna_char_to_nuc(contig[nxtbse]);
      bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
      tmp_key = binary_kmer_get_key(bkmer, kmer_size);
      node = hash_table_find(&db_graph->ht, tmp_key);

      if(node != HASH_NOT_FOUND &&
         (colour == -1 || db_node_has_col(db_graph, node, colour)))
      {
        nodes->data[n].key = node;
        nodes->data[n].orient = bkmer_get_orientation(bkmer, tmp_key);
        rpos->data[n] = offset;
        n++;
      }
    }
  }

  // Return number of bases from the last kmer found until read end
  size_t ret = (n == init_len ? r->seq.end /* No kmers found */
                              : r->seq.end - (rpos->data[n-1] + kmer_size));

  nodes->len = rpos->len = n;

  // Check for sequence gaps
  for(i = init_len; i+1 < nodes->len; i++) {
    if(rpos->data[i]+1 < rpos->data[i+1]) {
      aln->seq_gaps = true;
      break;
    }
  }

  return ret;
}


// if colour is -1 aligns to all colours, otherwise aligns to given colour only
// Assumes both reads are in FF orientation
void db_alignment_from_reads(dBAlignment *alignment,
                             const read_t *r1, const read_t *r2,
                             uint8_t qcutoff1, uint8_t qcutoff2,
                             uint8_t hp_cutoff,
                             const dBGraph *db_graph, int colour)
{
  ctx_assert(colour == -1 || db_graph->node_in_cols != NULL);

  db_node_buf_reset(&alignment->nodes);
  int32_buf_reset(&alignment->rpos);
  alignment->seq_gaps = false;
  alignment->r2enderr = 0;
  alignment->passed_r2 = (r2 != NULL);
  alignment->colour = colour;
  alignment->r1bases = r1->seq.end;
  alignment->r2bases = r2 ? r2->seq.end : 0;

  alignment->r1enderr = db_alignment_from_read(alignment, r1,
                                               qcutoff1, hp_cutoff,
                                               db_graph, colour);

  alignment->r2strtidx = alignment->nodes.len;

  if(r2 != NULL) {
    alignment->r2enderr = db_alignment_from_read(alignment, r2,
                                                 qcutoff2, hp_cutoff,
                                                 db_graph, colour);
  }

  alignment->used_r1 = (alignment->r1enderr < r1->seq.end);
  alignment->used_r2 = (r2 != NULL && alignment->r2enderr < r2->seq.end);

  #ifdef CTXVERBOSE
    db_alignment_print(alignment);
  #endif
}

// Returns index of node just after next gap,
// or aln->nodes.len if no more gaps
size_t db_alignment_next_gap(const dBAlignment *aln, size_t start)
{
  size_t i, end = aln->rpos.len;
  int32_t *rpos = aln->rpos.data;

  if(end == 0) return 0;

  if(aln->used_r1 && aln->used_r2 && start < aln->r2strtidx)
    end = aln->r2strtidx;

  for(i = start+1; i < end && rpos[i-1]+1 == rpos[i]; i++) {}

  // Return position after gap
  return i;
}

// @return true iff alignment has no gaps, all kmers found
bool db_alignment_is_perfect(const dBAlignment *aln)
{
  return (!aln->passed_r2 && aln->used_r1 && // we only given one read and used it
          !aln->seq_gaps && // there are no gaps
          aln->rpos.len > 0 && aln->rpos.data[0] == 0 && // no missing at start
          aln->r1enderr == 0); // no missing kmers at the end
}

//
// Debugging
//

void db_alignment_print(const dBAlignment *aln)
{
  pthread_mutex_lock(&ctx_biglock);

  printf("dBAlignment:\n");
  printf("  r1bases: %zu r2bases: %zu\n",   aln->r1bases,  aln->r2bases);
  printf("  r1enderr: %zu r2enderr: %zu\n", aln->r1enderr, aln->r2enderr);
  printf("  passed_r2: %s r2strtidx: %zu\n", aln->passed_r2 ? "yes" : "no", aln->r2strtidx);
  printf("  used_r1: %s used_r2: %s\n", aln->used_r1 ? "yes" : "no", aln->used_r2 ? "yes" : "no");
  printf("  seq_gaps: %s\n", aln->seq_gaps ? "yes" : "no");
  printf("  colour: %i\n", aln->colour);

  size_t i;

  for(i = 0; i < aln->nodes.len; i++) {
    printf("  [%zu] %i: %zu:%i\n", i, aln->rpos.data[i],
           (size_t)aln->nodes.data[i].key, aln->nodes.data[i].orient);
  }

  pthread_mutex_unlock(&ctx_biglock);
}

// Check all edges between ungapped adjacent nodes
bool db_alignment_check_edges(const dBAlignment *aln, const dBGraph *graph)
{
  size_t start, end;

  for(start = 0; start < aln->nodes.len; start = end)
  {
    end = db_alignment_next_gap(aln, start);
    if(db_graph_check_all_edges(graph, aln->nodes.data+start, end-start))
      return false;
  }

  return true;
}
