#include "global.h"

#include <pthread.h>
#include "msg-pool/msgpool.h" // pool for getting jobs

#include "generate_paths.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "path_store.h"
#include "path_format.h"
#include "graph_paths.h"
#include "db_alignment.h"
#include "correct_alignment.h"
#include "async_read_io.h"
#include "seq_reader.h"
#include "binary_seq.h"

//
// Multithreaded code to add paths to the graph from sequence data
// Uses async_read_io to allow multiple readers, multiple workers
//

// #define CTXVERBOSE 1

struct GenPathWorker
{
  pthread_t thread;
  dBGraph *const db_graph;

  // We take jobs from the pool
  MsgPool *pool;
  volatile size_t *rcounter;

  // Current job
  AsyncIOData *data; // current data
  CorrectAlnReadsTask task; // current task
  LoadingStats stats;

  dBAlignment aln;
  CorrectAlnWorker corrector;

  // Nucleotides and positions of junctions
  // only one array allocated for each type, rev points to half way through
  uint8_t *pck_fw, *pck_rv;
  size_t *pos_fw, *pos_rv;
  size_t num_fw, num_rv, junc_arrsize;

  FILE *tmp_fh;
};

#define INIT_BUFLEN 1024

// Printing variables defined in correct_reads_input.h
// Used for printint output
volatile size_t print_contig_id = 0, print_path_id = 0;


size_t gen_paths_worker_est_mem(const dBGraph *db_graph)
{
  size_t job_mem, corrector_mem, junc_mem, packed_mem;
  job_mem = 1024*4; // Assume 1024 bytes per read, 2 reads, seq+qual
  corrector_mem = correct_aln_worker_est_mem(db_graph);
  junc_mem = 2 * INIT_BUFLEN * (sizeof(Nucleotide)+sizeof(size_t));
  packed_mem = INIT_BUFLEN;

  return job_mem + corrector_mem + junc_mem + packed_mem + sizeof(GenPathWorker);
}

#define binary_seq_mem(n) ((((n)+3)/4 + sizeof(PathLen))*4)

static void _gen_paths_worker_alloc(GenPathWorker *wrkr, dBGraph *db_graph,
                                    FILE *tmp_fh)
{
  GenPathWorker tmp = {.db_graph = db_graph, .pool = NULL,
                       .tmp_fh = tmp_fh};

  db_alignment_alloc(&tmp.aln);
  correct_aln_worker_alloc(&tmp.corrector, db_graph);
  loading_stats_init(&tmp.stats);

  // Junction data
  // only fw arrays are malloc'd, rv point to fw
  tmp.junc_arrsize = INIT_BUFLEN;
  tmp.pck_fw = ctx_calloc(2, binary_seq_mem(tmp.junc_arrsize));
  tmp.pck_rv = tmp.pck_fw + binary_seq_mem(tmp.junc_arrsize);

  tmp.pos_fw = ctx_malloc(tmp.junc_arrsize * sizeof(*tmp.pos_fw) * 2);
  tmp.pos_rv = tmp.pos_fw + tmp.junc_arrsize;
  tmp.num_fw = tmp.num_rv = 0;

  memcpy(wrkr, &tmp, sizeof(GenPathWorker));
}

static void _gen_paths_worker_dealloc(GenPathWorker *wrkr)
{
  db_alignment_dealloc(&wrkr->aln);
  correct_aln_worker_dealloc(&wrkr->corrector);
  ctx_free(wrkr->pck_fw);
  ctx_free(wrkr->pos_fw);
}


GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph, FILE *fout)
{
  size_t i;
  GenPathWorker *workers = ctx_malloc(n * sizeof(GenPathWorker));
  for(i = 0; i < n; i++) _gen_paths_worker_alloc(&workers[i], graph, fout);
  return workers;
}

void gen_paths_workers_dealloc(GenPathWorker *workers, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) _gen_paths_worker_dealloc(&workers[i]);
  ctx_free(workers);
}

static inline void worker_nuc_cap(GenPathWorker *wrkr, size_t req_cap)
{
  size_t old_cap = wrkr->junc_arrsize, old_pck_mem, new_pck_mem;

  if(req_cap > old_cap)
  {
    wrkr->junc_arrsize = roundup2pow(req_cap);

    // use recalloc which zeros new memory
    // Zero new memory to keep valgrind happy :(
    old_pck_mem = binary_seq_mem(old_cap);
    new_pck_mem = binary_seq_mem(wrkr->junc_arrsize);
    wrkr->pck_fw = ctx_recalloc(wrkr->pck_fw, old_pck_mem*2, new_pck_mem*2);
    wrkr->pck_rv = wrkr->pck_fw + new_pck_mem;

    wrkr->pos_fw = ctx_realloc(wrkr->pos_fw, 2 * wrkr->junc_arrsize * sizeof(size_t));
    wrkr->pos_rv = wrkr->pos_fw + wrkr->junc_arrsize;
  }
}

// Merge stats into workers[0]
static void generate_paths_merge_stats(GenPathWorker *wrkrs, size_t num_workers)
{
  size_t i;
  for(i = 1; i < num_workers; i++) {
    correct_aln_stats_merge(&wrkrs[0].corrector.gapstats, &wrkrs[i].corrector.gapstats);
    correct_aln_stats_zero(&wrkrs[i].corrector.gapstats);
    loading_stats_merge(&wrkrs[0].stats, &wrkrs[i].stats);
    loading_stats_init(&wrkrs[i].stats);
  }
}

/* Get stats */

CorrectAlnStats gen_paths_get_gapstats(GenPathWorker *wrkr)
{
  return wrkr->corrector.gapstats;
}

LoadingStats gen_paths_get_stats(const GenPathWorker *wrkr)
{
  return wrkr->stats;
}

// Returns number of paths added
// `pos_pl` is an array of positions in the nodes array of nodes to add paths to
// `packed_ptr` is <plen><seq> and is the nucleotides denoting this path
//   nuc at [3] is the nucleotide to take leaving nodes[pos_pl[3]]
// `pos_mn` is an array of junctions running in the opposite direction
// both pos_pl and pos_mn are sorted in their DIRECTION
// if pl_is_fw,  pos_pl is 1,3,5,6,12, pos_mn is 15,11,5,2
// if !pl_is_fw, pos_pl is 15,11,5,2,  pos_mn is 1,3,5,6,12
static inline size_t _juncs_to_paths(const size_t *restrict pos_pl,
                                     const size_t *restrict pos_mn,
                                     const size_t num_pl, const size_t num_mn,
                                     uint8_t *packed_ptr,
                                     const bool pl_is_fw,
                                     const dBNode *nodes,
                                     GenPathWorker *wrkr)
{
  size_t i, num_added = 0;
  const size_t ctpcol = wrkr->task.crt_params.ctpcol;

  dBGraph *db_graph = wrkr->db_graph;
  dBNode node;
  size_t start_mn, start_pl, pos;
  PathLen plen, plen_orient;
  bool added;
  PathIndex pindex = 0; // address of path once added
  bool printed = false;

  #ifdef CTXVERBOSE
    // char str[num_pl+1];
    // for(i = 0; i < num_pl; i++) str[i] = dna_nuc_to_char(nuc_pl[i]);
    // str[num_pl] = '\0';
    // status("[addpath] %s %s", pl_is_fw ? "fw" : "rv", str);
  #endif

  // <plen><seq> is in packed_ptr
  // create packed path with remaining 3 diff offsets (1..3)
  size_t pckd_memsize = sizeof(PathLen) + (num_pl+3)/4;
  uint8_t *packed_ptrs[4], *pckd = packed_ptr+sizeof(PathLen);

  for(i = 0; i < 4; i++) packed_ptrs[i] = packed_ptr + i*pckd_memsize;

  binary_seq_cpy(packed_ptrs[1]+sizeof(PathLen), pckd, 1, num_pl);
  binary_seq_cpy(packed_ptrs[2]+sizeof(PathLen), pckd, 2, num_pl);
  binary_seq_cpy(packed_ptrs[3]+sizeof(PathLen), pckd, 3, num_pl);

  // pl => plus in direction
  // mn => minus against direction

  // pl_is_fw:
  //  pos_pl: 0,1,2,3
  //  pos_mn: 3,2,1
  //!pl_is_fw:
  //  pos_pl: 3,2,1
  //  pos_mn: 0,1,2,3

  // start_pl if the first base in the path
  // nodes[pos_mn[start_mn]] is the node before which we add the path

  // Add paths longest -> shortest
  for(start_pl = 0, start_mn = num_mn-1; start_mn != SIZE_MAX; start_mn--)
  {
    if(pl_is_fw)
      while(start_pl < num_pl && pos_pl[start_pl] < pos_mn[start_mn]) start_pl++;
    else
      while(start_pl < num_pl && pos_pl[start_pl] > pos_mn[start_mn]) start_pl++;

    if(start_pl == num_pl) break;

    // nodes[pos] is the kmer we are adding this path to
    pos = (pl_is_fw ? pos_mn[start_mn] - 1 : pos_mn[start_mn] + 1);
    node = (pl_is_fw ? nodes[pos] : db_node_reverse(nodes[pos]));

    start_pl -= (start_pl > 0 && pos_pl[start_pl-1] == pos);

    // ^ deals with the case where >junction on kmer JUST before <junction:
    // e.g. backtrack to add the 'F'
    //
    // bCD ----+      +---> FGh
    //         |      |
    //         v      |           FW [BCD->FGH]: FH (added F)
    // BCD -> DEF -> EFG -> FGH   RV [FGH->BCD]: complement(B)
    //  |
    //  +---> DEf
    //

    // Check path is not too long (MAX_PATHLEN is the limit)
    plen = (PathLen)MIN2(num_pl - start_pl, MAX_PATHLEN);

    #ifdef CTXVERBOSE
      char kmerstr[MAX_KMER_SIZE+1];
      BinaryKmer tmpkmer = db_node_get_bkmer(db_graph, node.key);
      binary_kmer_to_str(tmpkmer, db_graph->kmer_size, kmerstr);
      printf(" %s:%i) start_pl: %zu start_mn: %zu {%zu}\n",
             kmerstr, node.orient, start_pl, start_mn, pos_mn[start_mn]);
    #endif

    // Write orient and length to packed representation
    plen_orient = packedpath_combine_lenorient(plen,node.orient);
    packed_ptr = packed_ptrs[start_pl&3] + start_pl/4;
    memcpy(packed_ptr, &plen_orient, sizeof(PathLen));

    // mask top byte!
    size_t top_idx = sizeof(PathLen) + (plen+3)/4 - 1;
    uint8_t top_byte = packed_ptr[top_idx];
    packed_ptr[top_idx] &= 0xff >> (8 - bits_in_top_byte(plen));

    added = graph_paths_find_or_add_mt(node, ctpcol, packed_ptr, plen,
                                       &db_graph->pstore, &pindex);

    packed_ptr[top_idx] = top_byte; // restore top byte

    #ifdef CTXVERBOSE
      printf("We %s\n", added ? "added" : "abandoned");
    #endif

    // If the path already exists, all of its subpaths also already exist
    if(!added && plen < MAX_PATHLEN) break;
    num_added++;

    if(gen_paths_print_paths && !printed)
    {
      // print path
      size_t start, end;
      if(pl_is_fw) { start = pos, end = pos_pl[num_pl-1]+1; }
      else { start = pos_pl[num_pl-1]-1, end = pos; }
      ctx_assert2(start < end, "start: %zu, end: %zu", start, end);

      pthread_mutex_lock(&biglock);
      fprintf(stdout, ">path%zu.%s\n", print_path_id++, pl_is_fw ? "fw" : "rv");
      db_nodes_print(nodes+start, end-start+1, db_graph, stdout);
      fputc('\n', stdout);
      db_nodes_print_edges(nodes+start, end-start+1, db_graph, stdout);
      fputc('\n', stdout);
      pthread_mutex_unlock(&biglock);

      printed = true;
    }

    #ifdef CTXCHECKS
      const size_t ctxcol = wrkr->task.crt_params.ctxcol;
      Colour cols[2] = {ctxcol, ctpcol};
      GraphPathPairing gp = {.ctxcols = cols, .ctpcols = cols+1, .n = 1};

      // Check path before we wrote it
      ctx_check2(graph_paths_check_valid(node, ctxcol, packed_ptr+sizeof(PathLen),
                                         plen, db_graph),
                 "read: %s %s", wrkr->data->r1.name.b, wrkr->data->r1.seq.b);

      // Check path after we wrote it
      ctx_check2(graph_paths_check_path(node.key, pindex, &gp, db_graph),
                 "read: %s %s", wrkr->data->r1.name.b, wrkr->data->r1.seq.b);
    #endif

    // status("Path is:...");
    // print_path(node.key, packed_ptr, &db_graph->pstore);
  }

  return num_added;
}

static void worker_junctions_to_paths(GenPathWorker *wrkr,
                                      const dBNodeBuffer *contig)
{
  size_t num_fw = wrkr->num_fw, num_rv = wrkr->num_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  uint8_t *pck_fw = wrkr->pck_fw, *pck_rv = wrkr->pck_rv;

  ctx_assert2(num_fw && num_rv, "num_fw: %zu num_rv: %zu", num_fw, num_rv);

  #ifdef CTXVERBOSE
    status("num_fw: %zu num_rv: %zu", num_fw, num_rv);
  #endif

  // Reverse pos_rv array; reverse nucleotides if needed later
  size_t i, j;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j]);
  }

  // _juncs_to_paths returns the number of paths added
  size_t n;
  const dBNode *nodes = contig->data;

  n = _juncs_to_paths(pos_fw, pos_rv, num_fw, num_rv, pck_fw, true, nodes, wrkr);
  if(n) {
    binary_seq_reverse_complement(pck_rv+sizeof(PathLen), num_rv);
    _juncs_to_paths(pos_rv, pos_fw, num_rv, num_fw, pck_rv, false, nodes, wrkr);
  }
}

static void worker_contig_to_junctions(GenPathWorker *wrkr,
                                       const dBNodeBuffer *contig)
{
  // status("nodebuf: %zu", wrkr->contig.len+MAX_KMER_SIZE+1);
  worker_nuc_cap(wrkr, contig->len);

  if(gen_paths_print_contigs) {
    pthread_mutex_lock(&biglock);
    fprintf(stdout, ">contig%zu\n", print_contig_id++);
    db_nodes_print(contig->data, contig->len, wrkr->db_graph, stdout);
    fputc('\n', stdout);
    pthread_mutex_unlock(&biglock);
  }

  // Find forks in this colour
  Edges edges;
  size_t  i, num_fw = 0, num_rv = 0;
  int indegree, outdegree;
  uint8_t *pck_fw, *pck_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide nuc;

  pck_fw = wrkr->pck_fw+sizeof(PathLen);
  pck_rv = wrkr->pck_rv+sizeof(PathLen);

  dBGraph *db_graph = wrkr->db_graph;

  const dBNode *nodes = contig->data;
  const size_t contig_len = contig->len;
  const size_t ctxcol = wrkr->task.crt_params.ctxcol;

  for(i = 0; i < contig_len; i++)
  {
    edges = db_node_get_edges(db_graph, nodes[i].key, ctxcol);
    outdegree = edges_get_outdegree(edges, nodes[i].orient);
    indegree = edges_get_indegree(edges, nodes[i].orient);

    if(indegree > 1 && i > 0)
    {
      nuc = db_node_get_first_nuc(nodes[i-1], db_graph);
      binary_seq_set(pck_rv, num_rv, nuc);
      pos_rv[num_rv++] = i;
    }

    // Only adding forward junctions after num_rv > 0 seems like a nice
    // optimisation but is actually a bad idea, in the case of outdegree > 1
    // just before first indegree >1, we need that junction
    if(outdegree > 1 && i+1 < contig_len)
    {
      nuc = db_node_get_last_nuc(nodes[i+1], db_graph);
      binary_seq_set(pck_fw, num_fw, nuc);
      pos_fw[num_fw++] = i;
    }
  }

  wrkr->num_fw = num_fw;
  wrkr->num_rv = num_rv;

  if(num_fw > 0 && num_rv > 0)
    worker_junctions_to_paths(wrkr, contig);
}

// wrkr->data and wrkr->task must be set before calling this functions
static void reads_to_paths(GenPathWorker *wrkr)
{
  AsyncIOData *data = wrkr->data;
  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  if(gen_paths_print_reads) {
    pthread_mutex_lock(&biglock);
    printf(">read %s %s\n%s %s\n", data->r1.name.b, data->r2.name.b,
                                   data->r1.seq.b, data->r2.seq.b);
    pthread_mutex_unlock(&biglock);
  }

  uint8_t fq_cutoff1, fq_cutoff2;
  fq_cutoff1 = fq_cutoff2 = wrkr->task.fq_cutoff;

  if(fq_cutoff1 > 0) {
    fq_cutoff1 += data->fq_offset1;
    fq_cutoff2 += data->fq_offset2;
  }

  uint8_t hp_cutoff = wrkr->task.hp_cutoff;

  // Second read is in reverse orientation - need in forward
  if(r2 != NULL)
    seq_reader_orient_mp_FF_or_RR(r1, r2, wrkr->task.matedir);

  // Update stats
  if(r2 == NULL) wrkr->stats.num_se_reads++;
  else wrkr->stats.num_pe_reads += 2;

  db_alignment_from_reads(&wrkr->aln, r1, r2,
                          fq_cutoff1, fq_cutoff2, hp_cutoff,
                          wrkr->db_graph, wrkr->task.crt_params.ctxcol);

  ctx_check2(db_alignment_check_edges(&wrkr->aln, wrkr->db_graph),
             "Edges missing: was read %s%s%s used to build the graph?",
             r1->name.b, r2 ? ", " : "", r2 ? r2->name.b : "");

  // For debugging
  // db_alignment_print(&wrkr->aln, wrkr->db_graph);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, wrkr->task.crt_params);

  dBNodeBuffer *nbuf;
  while((nbuf = correct_alignment_nxt(&wrkr->corrector)) != NULL)
    worker_contig_to_junctions(wrkr, nbuf);
}

// Print progress every 5M reads
#define REPORT_RATE 5000000
// Defragment collection every 10M reads
#define DEFRAG_RATE 10000000

static void gen_paths_print_progress(size_t n)
{
  if(n % REPORT_RATE == 0)
  {
    char num_str[100];
    long_to_str(n, num_str);
    status("[GenPaths] Read %s entries (reads / read pairs)", num_str);
  }
}

// pthread method, loop: grabs job, does processing
static void generate_paths_worker(void *ptr)
{
  GenPathWorker *wrkr = (GenPathWorker*)ptr;
  MsgPool *pool = wrkr->pool;
  AsyncIOData *data;
  int pos;

  while((pos = msgpool_claim_read(pool)) != -1)
  {
    memcpy(&data, msgpool_get_ptr(pool, pos), sizeof(AsyncIOData*));
    wrkr->data = data;
    memcpy(&wrkr->task, data->ptr, sizeof(CorrectAlnReadsTask));
    reads_to_paths(wrkr);
    msgpool_release(pool, pos, MPOOL_EMPTY);

    // Print progress
    size_t n = __sync_add_and_fetch(wrkr->rcounter, 1);
    gen_paths_print_progress(n);
  }
}

void gen_paths_worker_seq(GenPathWorker *wrkr, AsyncIOData *data,
                          const CorrectAlnReadsTask *task)
{
  // Copy task to worker
  wrkr->data = data;
  memcpy(&wrkr->task, task, sizeof(CorrectAlnReadsTask));

  reads_to_paths(wrkr);
}

// Function used in tests
void gen_paths_from_str_mt(GenPathWorker *gen_path_wrkr, char *seq,
                           CorrectAlnParam params)
{
  // Fake reads
  size_t seqlen = strlen(seq);
  char empty[10] = "", rname[20] = "Example";
  read_t r1 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = seq, .end = seqlen, .size = seqlen+1},
               .qual = {.b = empty, .end = 0, .size = 1}};

  read_t r2 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = empty, .end = 0, .size = 1},
               .qual = {.b = empty, .end = 0, .size = 1}};

  AsyncIOData iodata = {.r1 = r1, .r2 = r2, .ptr = NULL,
                        .fq_offset1 = 0, .fq_offset2 = 2};

  AsyncIOReadTask iotask = {.file1 = NULL, .file2 = NULL,
                            .fq_offset = 0, .interleaved = false};

  CorrectAlnReadsTask task = {.files = iotask, .fq_cutoff = 0, .hp_cutoff = 0,
                              .matedir = READPAIR_FF, .crt_params = params,
                              .ptr = NULL};

  gen_paths_worker_seq(gen_path_wrkr, &iodata, &task);
}

void generate_paths(CorrectAlnReadsTask *tasks, size_t num_inputs,
                    GenPathWorker *workers, size_t num_workers)
{
  size_t i;

  AsyncIOData *data = ctx_malloc(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  size_t read_counter = 0;
  MsgPool pool;
  msgpool_alloc(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*), USE_MSG_POOL);
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  for(i = 0; i < num_workers; i++) {
    workers[i].pool = &pool;
    workers[i].rcounter = &read_counter;
  }

  AsyncIOReadTask *asyncio_tasks = ctx_malloc(num_inputs * sizeof(AsyncIOReadTask));
  correct_reads_input_to_asycio(asyncio_tasks, tasks, num_inputs);

  asyncio_run_threads(&pool, asyncio_tasks, num_inputs, generate_paths_worker,
                      workers, num_workers, sizeof(GenPathWorker));

  ctx_free(asyncio_tasks);
  msgpool_dealloc(&pool);

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  ctx_free(data);

  // Merge gap counts into worker[0]
  generate_paths_merge_stats(workers, num_workers);
}
