#include "global.h"
#include "breakpoint_caller.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "multicol_walker.h"

typedef struct {
  size_t i, start, length;
} Chrom;

typedef struct {
  uint64_t start, count;
} KOccurList;

typedef struct
{
  uint64_t orient:1, pos:63;
} KOccur;

// Malloc's and returns chrom list - remember to free
static Chrom* generate_chrom_list(const read_t *reads, size_t num_reads)
{
  Chrom *chroms = malloc2(num_reads*sizeof(Chrom));
  size_t i, offset = 0;

  for(i = 0; i < num_reads; i++)
  {
    chroms[i] = (Chrom){.i = i, .start = offset, .length = reads[i].seq.end};
    offset += reads[i].seq.end;
  }

  return chroms;
}

static void bkmer_update_counts(BinaryKmer bkmer, dBGraph *db_graph,
                                KOccurList *klists)
{
  bool found;
  hkey_t hkey = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
  klists[hkey].count++;
}

static void read_update_counts(const read_t *r, dBGraph *db_graph,
                               KOccurList *klists)
{
  const size_t kmer_size = db_graph->kmer_size;
  LoadingStats stats;
  READ_TO_BKMERS(r, kmer_size, 0, 0, &stats, bkmer_update_counts,
                 db_graph, klists);
}

static void bkmer_store_kmer_pos(BinaryKmer bkmer, const dBGraph *db_graph,
                                 KOccurList *klists, KOccur *koccurs,
                                 uint64_t pos)
{
  KOccurList *kl;
  // bkmers were already added to graph -> don't need to find_or_insert
  dBNode node = db_graph_find(db_graph, bkmer);
  ctx_assert(node.key != HASH_NOT_FOUND);

  kl = &klists[node.key];
  koccurs[kl->start + kl->count] = (KOccur){.pos = pos, .orient = node.orient};
  kl->count++;
}

static void read_store_kmer_pos(const read_t *r, const dBGraph *db_graph,
                                KOccurList *klists, KOccur *koccurs,
                                uint64_t roffset)
{
  const size_t kmer_size = db_graph->kmer_size;
  LoadingStats stats;
  READ_TO_BKMERS(r, kmer_size, 0, 0, &stats, bkmer_store_kmer_pos,
                 db_graph, klists, koccurs, roffset+_offset);
}

static void process_contig(const dBGraph *db_graph,
                           const KOccurList *klists, const KOccur *koccurs,
                           MulticolWalker *walker,
                           size_t *cols_used, size_t *col_lengths,
                           dBNodeBuffer *nbuf, gzFile gzout)
{
  // Work backwards to find last place we met the ref
  // nbuf[0] is ref node
  // nbuf[1] is first node not in ref
  size_t i, end;
  for(i = nbuf->len-1; i > 1; i--)
    if(klists[nbuf->data[i].key].count > 0)
      break;

  if(i == 1) return;

  // Found second ref meeting - find first node in block that is in ref
  end = i;
  while(klists[nbuf->data[end-1].key].count > 0) end--;

  gzprintf(gzout, ">call.5pflank\n");
  // DEV
  gzprintf(gzout, "\n");

  gzprintf(gzout, ">call.3pflank\n");
  db_nodes_gzprint_cont(nbuf->data+end, nbuf->len - end, db_graph, gzout);

  gzprintf(gzout, ">call.path\n");
  db_nodes_gzprint_cont(nbuf->data, end, db_graph, gzout);
  // DEV
}

// DEV: for speed could remember where colours split off
static void assemble_contig(const dBGraph *db_graph,
                            const KOccurList *klists, const KOccur *koccurs,
                            MulticolWalker *walker,
                            size_t *cols, size_t *cols_used, size_t *col_lengths,
                            dBNodeBuffer *nbuf, gzFile gzout)
{
  // Follow path using multicol_walker_assemble_contig
  size_t i, ncols = db_graph->num_of_cols, ncols_used;
  for(i = 0; i < ncols; i++) cols[i] = i;

  while(ncols > 0)
  {
    ncols_used = multicol_walker_assemble_contig(walker, cols, ncols,
                                                 cols_used, col_lengths, nbuf);

    // Print with path
    process_contig(db_graph, klists, koccurs, walker, cols_used, col_lengths,
                   nbuf, gzout);

    // Remove the used colours from list of colours to use
    ncols = multicol_walker_rem_cols(cols, ncols, cols_used, ncols_used);
  }
}

// Walk the graph remembering the last time we met the ref
// When traversal fails, dump sequence up to last meeting with the ref
static void follow_break(dBNode node, const dBGraph *db_graph,
                         const KOccurList *klists, const KOccur *koccurs,
                         MulticolWalker *walker,
                         size_t *cols, size_t *cols2, size_t *cols3,
                         dBNodeBuffer *nbuf, gzFile gzout)
{
  size_t i, j, num_next;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);
  Edges edges = db_node_get_edges(db_graph, node.key, 0);

  num_next = db_graph_next_nodes(db_graph, bkey, node.orient, edges,
                                 next_nodes, next_nucs);

  // Filter out next nodes in the reference
  for(i = 0, j = 0; i < num_next; i++) {
    if(klists[next_nodes[i].key].count > 0) {
      next_nodes[j] = next_nodes[i];
      next_nucs[j] = next_nucs[i];
      j++;
    }
  }

  // Abandon if all options are in ref or none are
  if(j == num_next || j == 0) return;

  num_next = j;

  // 2. Follow all paths not in ref, in all colours
  for(i = 0; i < num_next; i++)
  {
    db_node_buf_reset(nbuf);
    nbuf->data[0] = node;
    nbuf->data[1] = next_nodes[i];
    nbuf->len = 2;
    assemble_contig(db_graph, klists, koccurs, walker, cols, cols2, cols3,
                    nbuf, gzout);
  }
}

struct BreakpointCaller {
  size_t idx;
  hkey_t start, end;
  const dBGraph *db_graph;
  const KOccurList *klists;
  const KOccur *koccurs;
  gzFile gzout;
  size_t *cols, *cols2, *cols3;
};

static void* call_breakpoints(void *ptr)
{
  const struct BreakpointCaller *caller = (struct BreakpointCaller*)ptr;
  const dBGraph *db_graph = caller->db_graph;
  const KOccurList *klists = caller->klists;
  const KOccur *koccurs = caller->koccurs;
  gzFile gzout = caller->gzout;
  size_t *cols = caller->cols;
  size_t *cols2 = caller->cols2;
  size_t *cols3 = caller->cols3;

  const BinaryKmer *table = db_graph->ht.table, *bkptr;
  const BinaryKmer *start = table + caller->start;
  const BinaryKmer *end = table + caller->end;

  ctx_assert(db_graph->num_edge_cols == 1);

  hkey_t hkey;
  Edges edges;

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 1024);

  MulticolWalker walker;
  multicol_walker_alloc(&walker, db_graph);

  for(bkptr = start; bkptr < end; bkptr++) {
    if(HASH_ENTRY_ASSIGNED(*bkptr)) {
      hkey = (hkey_t)(bkptr - table);
      // check node is in the ref
      if(klists[hkey].count > 0) {
        edges = db_node_get_edges(db_graph, hkey, 0);
        if(edges_get_outdegree(edges, FORWARD) > 1) {
          follow_break((dBNode){.key = hkey, .orient = FORWARD}, db_graph,
                       klists, koccurs, &walker, cols, cols2, cols3, &nbuf,
                       gzout);
        }
        if(edges_get_outdegree(edges, REVERSE) > 1) {
          follow_break((dBNode){.key = hkey, .orient = REVERSE}, db_graph,
                       klists, koccurs, &walker, cols, cols2, cols3, &nbuf,
                       gzout);
        }
      }
    }
  }

  multicol_walker_dealloc(&walker);
  db_node_buf_dealloc(&nbuf);

  return NULL;
}

static void breakpoints_print_header(gzFile gzout, const CmdArgs *args)
{
  char datestr[9], cwd[PATH_MAX + 1];

  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  gzprintf(gzout, "##fileFormat=CtxBreakpointsv0.1\n##fileDate=%s\n", datestr);
  gzprintf(gzout, "##cmd=\"%s\"\n", args->cmdline);
  if(futil_get_current_dir(cwd) != NULL) gzprintf(gzout, "##wkdir=%s\n", cwd);

  gzprintf(gzout, "##reference=TODO\n");
  gzprintf(gzout, "##ctxVersion=\""VERSION_STATUS_STR"\">\n");
  gzprintf(gzout, "##ctxKmerSize=%i", MAX_KMER_SIZE);
}

static void breakpoints_multithreaded(gzFile gzout, const dBGraph *db_graph,
                                      Chrom *chroms, size_t num_chroms,
                                      KOccurList *klists, KOccur *koccurs,
                                      gzFile *tmp_files, size_t num_of_threads,
                                      const CmdArgs *args)
{
  ctx_assert((num_of_threads == 1) == (tmp_files == NULL));

  // Print header to gzout
  breakpoints_print_header(gzout, args);

  const size_t ncols = db_graph->num_of_cols;
  struct BreakpointCaller *callers;
  callers = malloc2(num_of_threads * sizeof(struct BreakpointCaller));
  size_t *colour_lists = malloc2(num_of_threads * ncols * 3 * sizeof(size_t));

  hkey_t start, end, step = db_graph->ht.capacity / num_of_threads;
  gzFile gztmp;
  int rc;
  size_t i, *colptr;

  for(i = 0; i < num_of_threads; i++)
  {
    start = i * step;
    end = (i == num_of_threads-1) ? db_graph->ht.capacity : start + step;
    gztmp = (tmp_files ? tmp_files[i] : gzout);
    colptr = colour_lists+i*ncols*3;

    callers[i]
      = (struct BreakpointCaller){.idx = i, .start = start, .end = end,
                                  .db_graph = db_graph, .klists = klists,
                                  .koccurs = koccurs, .gzout = gztmp,
                                  .cols = colptr, .cols2 = colptr+ncols,
                                  .cols3 = colptr+ncols*2};
  }

  if(num_of_threads <= 1) {
    call_breakpoints(&callers[0]);
  }
  else
  {
    // Multithreaded version
    /* Initialize and set thread detached attribute */
    pthread_t *threads = malloc2(num_of_threads * sizeof(pthread_t));
    pthread_attr_t thread_attr;
    pthread_attr_init(&thread_attr);
    pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

    for(i = 0; i < num_of_threads; i++) {
      rc = pthread_create(&threads[i], &thread_attr,
                          call_breakpoints, (void*)&callers[i]);
      if(rc != 0) die("Creating thread failed");
    }

    /* wait for all threads to complete */
    for(i = 0; i < num_of_threads; i++) {
      rc = pthread_join(threads[i], NULL);
      if(rc != 0) die("Joining thread failed");
    }

    // Merge files and clean up
    futil_merge_tmp_gzfiles(tmp_files, num_of_threads, gzout);
    pthread_attr_destroy(&thread_attr);
    free(threads);
  }

  free(colour_lists);
  free(callers);
}

void breakpoints_call(gzFile gzout, dBGraph *db_graph,
                      const read_t *reads, size_t num_reads,
                      size_t num_of_threads, const CmdArgs *args)
{
  // Set up temporary files
  gzFile *tmp_files = NULL;

  if(num_of_threads > 1) {
    tmp_files = futil_create_tmp_gzfiles(num_of_threads);
  }

  Chrom *chroms = generate_chrom_list(reads, num_reads);
  KOccurList *klists = calloc2(db_graph->ht.capacity, sizeof(KOccurList));

  // 1. Loop through reads, add to graph and record kmer counts
  size_t i;
  for(i = 0; i < num_reads; i++)
    read_update_counts(&reads[i], db_graph, klists);

  // 2. alloc lists
  uint64_t offset = 0;
  for(i = 0; i < num_reads; i++) {
    klists[i].start = offset;
    offset += klists[i].count;
    klists[i].count = 0;
  }

  KOccur *koccurs = malloc2(offset * sizeof(KOccur));

  // 3. Loop through reads, record kmer pos
  for(i = 0; i < num_reads; i++)
    read_store_kmer_pos(&reads[i], db_graph, klists, koccurs, klists[i].start);

  // 4. Call
  breakpoints_multithreaded(gzout, db_graph, chroms, num_reads,
                            klists, koccurs, tmp_files, num_of_threads, args);

  if(tmp_files) free(tmp_files);

  free(klists);
  free(koccurs);
  free(chroms);
}
