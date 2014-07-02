#include "global.h"
#include "db_alignment.h"
#include "seq_reader.h"


#define INIT_BUFLEN 1024

// Estimate initial memory required
size_t db_alignment_est_mem()
{
  return (sizeof(size_t)+sizeof(dBNode))*INIT_BUFLEN;
}

void db_alignment_alloc(dBAlignment *alignment)
{
  db_node_buf_alloc(&alignment->nodes, INIT_BUFLEN);
  uint32_buf_alloc(&alignment->gaps, INIT_BUFLEN);
}

void db_alignment_dealloc(dBAlignment *alignment)
{
  db_node_buf_dealloc(&alignment->nodes);
  uint32_buf_dealloc(&alignment->gaps);
}

// if colour is -1 aligns to all colours, otherwise aligns to given colour only
// Returns number of kmers lost from the end
static size_t db_alignment_from_read(dBAlignment *alignment, const read_t *r,
                                     uint8_t qcutoff, uint8_t hp_cutoff,
                                     const dBGraph *db_graph, int colour)
{
  size_t contig_start, contig_end = 0, search_start = 0, nxt_exp_kmer_offset = 0;
  const size_t kmer_size = db_graph->kmer_size;

  BinaryKmer bkmer, tmp_key;
  Nucleotide nuc;
  hkey_t node;
  size_t offset, nxtbse;

  dBNodeBuffer *nodes = &alignment->nodes;
  Uint32Buffer *gaps = &alignment->gaps;

  ctx_assert(nodes->len == gaps->len);
  size_t n = gaps->len;

  db_node_buf_ensure_capacity(nodes, n + r->seq.end);
  uint32_buf_ensure_capacity(gaps, n + r->seq.end);

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qcutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qcutoff, hp_cutoff, &search_start);

    const char *contig = r->seq.b + contig_start;
    size_t contig_len = contig_end - contig_start;

    bkmer = binary_kmer_from_str(contig, kmer_size);
    bkmer = binary_kmer_right_shift_one_base(bkmer);

    for(offset=0, nxtbse=kmer_size-1; nxtbse < contig_len; nxtbse++,offset++)
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
        gaps->data[n] = (uint32_t)(contig_start+offset - nxt_exp_kmer_offset);
        nxt_exp_kmer_offset = contig_start+offset+1;
        n++;
      }
      else alignment->seq_gaps = true;
    }
  }

  nodes->len = gaps->len = n;

  // If contig_end == 0 we never found any kmers
  size_t end_of_last_kmer = nxt_exp_kmer_offset+kmer_size-1;
  return contig_end == 0 ? r->seq.end : r->seq.end - end_of_last_kmer;
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
  uint32_buf_reset(&alignment->gaps);
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
    db_alignment_print(alignment, db_graph);
  #endif
}

// Returns index of node just after next gap,
// or aln->nodes.len if no more gaps
size_t db_alignment_next_gap(const dBAlignment *aln, size_t start)
{
  size_t i, end = aln->nodes.len;
  if(end == 0) return 0;

  if(aln->used_r1 && aln->used_r2 && start < aln->r2strtidx)
    end = aln->r2strtidx;

  for(i = start+1; i < end && aln->gaps.data[i] == 0; i++) {}
  return i;
}

//
// Debugging
//

void db_alignment_print(const dBAlignment *aln, const dBGraph *db_graph)
{
  pthread_mutex_lock(&ctx_biglock);

  printf("dBAlignment:\n");
  size_t i, start = 0, end = db_alignment_next_gap(aln, 0);
  while(start < aln->nodes.len)
  {
    if(start == aln->r2strtidx)
      printf("    gap: %zu -[ins]- %u\n", aln->r1enderr, aln->gaps.data[start]);
    else if(start == 0)
      printf("    start gap: %u\n", aln->gaps.data[start]);
    else
      printf("    gap: %u\n", aln->gaps.data[start]);

    printf("  %zu nodes\n", end-start);

    for(i = start; i < end; i++)
      printf(" %zu:%i", (size_t)aln->nodes.data[i].key, (int)aln->nodes.data[i].orient);
    printf(":\n");
    db_nodes_print_verbose(aln->nodes.data+start,end-start,db_graph,stdout);
    db_nodes_print_edges(aln->nodes.data+start,end-start,db_graph,stdout);
    printf("\n");

    start = end;
    end = db_alignment_next_gap(aln, end);
  }
  if(aln->passed_r2) {
    if(!aln->used_r2) printf("    [ins] unused r2: %u\n", aln->gaps.data[end]);
    else printf(" end gap [r2]: %zu\n", aln->r2enderr);
  }
  else printf(" end gap [r1]: %zu\n", aln->r1enderr);

  pthread_mutex_unlock(&ctx_biglock);
}

// Check all edges between ungapped adjacent nodes
bool db_alignment_check_edges(const dBAlignment *aln, const dBGraph *graph)
{
  size_t start, end, i;

  for(start = 0; start < aln->nodes.len; start = end) {
    end = db_alignment_next_gap(aln, start);

    // check edge between i and i+1
    for(i = start; i+1 < end; i++) {
      if(!db_graph_check_edges(graph, aln->nodes.data[i], aln->nodes.data[i+1]))
        return false;
    }
  }

  return true;
}
