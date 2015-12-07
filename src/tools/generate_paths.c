#include "global.h"

#include <pthread.h>
#include "msg-pool/msgpool.h" // pool for getting jobs

#include "generate_paths.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "db_alignment.h"
#include "correct_alignment.h"
#include "async_read_io.h"
#include "seq_reader.h"
#include "binary_seq.h"
#include "gpath_checks.h"

//
// Multithreaded code to add paths to the graph from sequence data
// Uses async_read_io to allow multiple readers, multiple workers
//

// #define CTXVERBOSE 1

#define GEN_PATHS_COUNTER_STEP 100
#define INIT_BUFLEN 1024

struct GenPathWorker
{
  pthread_t thread;
  dBGraph *const db_graph;

  // We take jobs from the pool
  // MsgPool *pool;
  size_t nreads;
  volatile size_t *shared_nreads;

  // Current job
  AsyncIOData *data; // current data
  CorrectAlnInput task; // current task

  CorrectAlnWorker corrector;

  // Nucleotides and positions of junctions
  // only one array allocated for each type, rev points to half way through
  uint8_t *pck_fw, *pck_rv;
  size_t *pos_fw, *pos_rv;
  size_t num_fw, num_rv, junc_arrsize;
};

// Printing variables defined in correct_aln_input.h
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

// GenPathWorker stores four buffers of size n
#define gworker_seq_buf(n) (binary_seq_mem(n)*4)

static void _gen_paths_worker_alloc(GenPathWorker *wrkr, dBGraph *db_graph)
{
  GenPathWorker tmp = {.db_graph = db_graph};

  correct_aln_worker_alloc(&tmp.corrector, true, db_graph);

  // Junction data
  // only fw arrays are malloc'd, rv point to fw
  tmp.junc_arrsize = INIT_BUFLEN;
  tmp.pck_fw = ctx_calloc(2, gworker_seq_buf(tmp.junc_arrsize));
  tmp.pck_rv = tmp.pck_fw + gworker_seq_buf(tmp.junc_arrsize);

  tmp.pos_fw = ctx_malloc(tmp.junc_arrsize * sizeof(*tmp.pos_fw) * 2);
  tmp.pos_rv = tmp.pos_fw + tmp.junc_arrsize;
  tmp.num_fw = tmp.num_rv = 0;

  memcpy(wrkr, &tmp, sizeof(GenPathWorker));
}

static void _gen_paths_worker_dealloc(GenPathWorker *wrkr)
{
  correct_aln_worker_dealloc(&wrkr->corrector);
  ctx_free(wrkr->pck_fw);
  ctx_free(wrkr->pos_fw);
}


GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph)
{
  size_t i;
  GenPathWorker *workers = ctx_malloc(n * sizeof(GenPathWorker));
  for(i = 0; i < n; i++)
    _gen_paths_worker_alloc(&workers[i], graph);
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
    old_pck_mem = gworker_seq_buf(old_cap);
    new_pck_mem = gworker_seq_buf(wrkr->junc_arrsize);
    wrkr->pck_fw = ctx_recallocarray(wrkr->pck_fw, old_pck_mem*2, new_pck_mem*2, 1);
    wrkr->pck_rv = wrkr->pck_fw + new_pck_mem;

    wrkr->pos_fw = ctx_realloc(wrkr->pos_fw, 2 * wrkr->junc_arrsize * sizeof(size_t));
    wrkr->pos_rv = wrkr->pos_fw + wrkr->junc_arrsize;
  }
}

/* Get stats */

CorrectAlnStats* gen_paths_get_aln_stats(GenPathWorker *wrkr)
{
  return &wrkr->corrector.aln_stats;
}

SeqLoadingStats* gen_paths_get_stats(GenPathWorker *wrkr)
{
  return &wrkr->corrector.load_stats;
}

// Returns number of paths added
// `pos_pl` is an array of positions in the nodes array of nodes to add paths to
// `packed_ptr` is <seq> and is the nucleotides denoting this path
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
  dBGraph *db_graph = wrkr->db_graph;
  const size_t ctpcol = wrkr->task.crt_params.ctpcol;
  GPathSet *gpset = &db_graph->gphash.gpstore->gpset;

  size_t i, num_added = 0;
  dBNode node;
  size_t start_mn, start_pl, pos;
  size_t plen;
  // bool added;
  // pkey_t pindex = 0; // address of path once added
  bool printed = false;

  #ifdef CTXVERBOSE
    // char str[num_pl+1];
    // for(i = 0; i < num_pl; i++) str[i] = dna_nuc_to_char(nuc_pl[i]);
    // str[num_pl] = '\0';
    // status("[addpath] %s %s", pl_is_fw ? "fw" : "rv", str);
  #endif

  // create packed path with remaining 3 diff offsets (1..3)
  size_t pckd_memsize = binary_seq_mem(num_pl);
  uint8_t *packed_ptrs[4], *pckd = packed_ptr;

  for(i = 0; i < 4; i++) packed_ptrs[i] = packed_ptr + i*pckd_memsize;

  binary_seq_cpy(packed_ptrs[1], pckd, 1, num_pl);
  binary_seq_cpy(packed_ptrs[2], pckd, 2, num_pl);
  binary_seq_cpy(packed_ptrs[3], pckd, 3, num_pl);

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

    // Check path is not too long (GPATH_MAX_JUNCS is the limit)
    plen = MIN2(num_pl - start_pl, GPATH_MAX_JUNCS);

    // Start and end nodes relative to `nodes` array
    size_t start, end; // index of first, last node in `nodes`
    if(pl_is_fw) { start = pos, end = pos_pl[start_pl+plen-1]+1; }
    else { start = pos_pl[start_pl+plen-1]-1, end = pos; }
    ctx_assert2(start < end, "start: %zu, end: %zu", start, end);

    // Write orient and length to packed representation
    packed_ptr = packed_ptrs[start_pl&3] + start_pl/4;

    // mask top byte!
    size_t top_idx = (plen+3)/4 - 1;
    uint8_t top_byte = packed_ptr[top_idx];
    packed_ptr[top_idx] &= 0xff >> (8 - bits_in_top_byte(plen));

    bool found = false;
    GPathNew newgpath = {.seq = packed_ptr,
                         .orient = node.orient, .num_juncs = plen,
                         .colset = NULL, .nseen = NULL};

    // #ifdef CTXVERBOSE
    //   char kmerstr[MAX_KMER_SIZE+1];
    //   BinaryKmer tmpkmer = db_node_get_bkmer(db_graph, node.key);
    //   binary_kmer_to_str(tmpkmer, db_graph->kmer_size, kmerstr);
    //   fputs(kmerstr, stdout);
    //   binary_seq_print(newgpath.seq, newgpath.num_juncs, stdout);
    //   fputc('\n', stdout);
    //   printf(" %s:%i) start_pl: %zu start_mn: %zu {%zu}\n",
    //          kmerstr, node.orient, start_pl, start_mn, pos_mn[start_mn]);
    // #endif

    GPath *gpath = gpath_hash_find_or_insert_mt(&db_graph->gphash, node.key,
                                                newgpath, &found);

    // Add colour
    bitset_set(gpath_get_colset(gpath, gpset->ncols), ctpcol);
    uint8_t *nseen = gpath_set_get_nseen(gpset, gpath);
    if(nseen != NULL) safe_add_uint8_mt(&nseen[ctpcol], 1);

    packed_ptr[top_idx] = top_byte; // restore top byte

    #ifdef CTXVERBOSE
      printf("We %s\n", found ? "abandoned" : "added");
    #endif

    // If the path already exists, all of its subpaths also already exist
    // if(found && plen < GPATH_MAX_JUNCS) break;
    num_added++;

    // Debugging
    if(gen_paths_print_paths && !printed)
    {
      // print path
      pthread_mutex_lock(&ctx_biglock);
      fprintf(stdout, ">path%zu.%s\n", print_path_id++, pl_is_fw ? "fw" : "rv");
      db_nodes_print(nodes+start, end-start+1, db_graph, stdout);
      fputc('\n', stdout);
      db_nodes_print_edges(nodes+start, end-start+1, db_graph, stdout);
      fputc('\n', stdout);
      pthread_mutex_unlock(&ctx_biglock);

      printed = true;
    }
  }

  return num_added;
}

static void worker_junctions_to_paths(GenPathWorker *wrkr,
                                      const dBNode *nodes, size_t num_nodes)
{
  (void)num_nodes;

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

  _juncs_to_paths(pos_fw, pos_rv, num_fw, num_rv, pck_fw, true, nodes, wrkr);
  binary_seq_reverse_complement(pck_rv, num_rv);
  _juncs_to_paths(pos_rv, pos_fw, num_rv, num_fw, pck_rv, false, nodes, wrkr);
}

static void worker_contig_to_junctions(GenPathWorker *wrkr,
                                       const dBNode *nodes, size_t num_nodes)
{
  // status("nodebuf: %zu", wrkr->contig.len+MAX_KMER_SIZE+1);
  worker_nuc_cap(wrkr, num_nodes);

  if(gen_paths_print_contigs) {
    pthread_mutex_lock(&ctx_biglock);
    fprintf(stdout, ">contig%zu\n", print_contig_id++);
    db_nodes_print(nodes, num_nodes, wrkr->db_graph, stdout);
    fputc('\n', stdout);
    pthread_mutex_unlock(&ctx_biglock);
  }

  // Find forks in this colour
  Edges edges;
  size_t  i, num_fw = 0, num_rv = 0;
  int indegree, outdegree;
  uint8_t *pck_fw = wrkr->pck_fw, *pck_rv = wrkr->pck_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide nuc;

  dBGraph *db_graph = wrkr->db_graph;
  const size_t ctxcol = wrkr->task.crt_params.ctxcol;

  for(i = 0; i < num_nodes; i++)
  {
    edges = db_node_get_edges(db_graph, nodes[i].key, ctxcol);
    outdegree = edges_get_outdegree(edges, nodes[i].orient);
    indegree = edges_get_indegree(edges, nodes[i].orient);

    if(i+1 < num_nodes)
    {
      nuc = db_node_get_last_nuc(nodes[i+1], db_graph);
      ctx_assert(edges_has_edge(edges, nuc, nodes[i].orient));

      // Only adding forward junctions after num_rv > 0 seems like a nice
      // optimisation but is actually a bad idea, in the case of outdegree > 1
      // just before first indegree >1, we need that junction
      if(outdegree > 1)
      {
        binary_seq_set(pck_fw, num_fw, nuc);
        pos_fw[num_fw++] = i;
      }
    }

    if(indegree > 1 && i > 0)
    {
      nuc = db_node_get_first_nuc(nodes[i-1], db_graph);
      binary_seq_set(pck_rv, num_rv, nuc);
      pos_rv[num_rv++] = i;
    }
  }

  wrkr->num_fw = num_fw;
  wrkr->num_rv = num_rv;

  if(num_fw > 0 && num_rv > 0)
    worker_junctions_to_paths(wrkr, nodes, num_nodes);
}

static bool error_printed_reads_overlap = false;

// wrkr->data and wrkr->task must be set before calling this functions
static void reads_to_paths(GenPathWorker *wrkr)
{
  AsyncIOData *data = wrkr->data;
  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  if(gen_paths_print_reads) {
    pthread_mutex_lock(&ctx_biglock);
    printf(">read %s %s\n%s %s\n", data->r1.name.b, data->r2.name.b,
                                   data->r1.seq.b, data->r2.seq.b);
    pthread_mutex_unlock(&ctx_biglock);
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

  // Check reads don't overlap in the middle of the fragment
  if(r1 && r2 && !error_printed_reads_overlap &&
     r1->seq.end + r2->seq.end > wrkr->task.crt_params.frag_len_min)
  {
    warn("Reads may overlap in fragment: %zu + %zu > frag len min: %u; max: %u",
         r1->seq.end, r2->seq.end,
         wrkr->task.crt_params.frag_len_min, wrkr->task.crt_params.frag_len_max);
    error_printed_reads_overlap = true;
  }

  correct_alignment_init(&wrkr->corrector, &wrkr->task.crt_params,
                         r1, r2, fq_cutoff1, fq_cutoff2, hp_cutoff);

  // ctx_check2(db_alignment_check_edges(&wrkr->corrector.aln, wrkr->db_graph),
  //            "Edges missing: was read %s%s%s used to build the graph?",
  //            r1->name.b, r2 ? ", " : "", r2 ? r2->name.b : "");

  dBNodeBuffer *nbuf;
  while((nbuf = correct_alignment_nxt(&wrkr->corrector)) != NULL)
  {
    worker_contig_to_junctions(wrkr, nbuf->b, nbuf->len);
  }
}

// pthread method, loop: grabs job, does processing
static void generate_paths_worker(AsyncIOData *data, size_t threadid, void *ptr)
{
  (void)threadid;
  GenPathWorker *wrkr = (GenPathWorker*)ptr;
  wrkr->data = data;
  memcpy(&wrkr->task, data->ptr, sizeof(CorrectAlnInput));
  reads_to_paths(wrkr);

  // Print progress
  wrkr->nreads++;
  if(wrkr->nreads >= GEN_PATHS_COUNTER_STEP) {
    // Update shared counter
    size_t n = __sync_fetch_and_add(wrkr->shared_nreads, wrkr->nreads);
    // if n .. n+wrkr->nreads
    ctx_update2("GenPaths", n, n+wrkr->nreads, CTX_UPDATE_REPORT_RATE);
    wrkr->nreads = 0;
  }
}

void gen_paths_worker_seq(GenPathWorker *wrkr, AsyncIOData *data,
                          const CorrectAlnInput *task)
{
  // Copy task to worker
  wrkr->data = data;
  memcpy(&wrkr->task, task, sizeof(CorrectAlnInput));

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
               .seq  = {.b = seq,   .end = seqlen, .size = seqlen+1},
               .qual = {.b = empty, .end = 0,      .size = 1}};

  read_t r2 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = empty, .end = 0, .size = 1},
               .qual = {.b = empty, .end = 0, .size = 1}};

  AsyncIOData iodata = {.r1 = r1, .r2 = r2, .ptr = NULL,
                        .fq_offset1 = 0, .fq_offset2 = 2};

  AsyncIOInput iotask = {.file1 = NULL, .file2 = NULL,
                         .fq_offset = 0, .interleaved = false};

  CorrectAlnInput task = CORRECT_ALN_INPUT_INIT;
  task.matedir = READPAIR_FF;
  task.crt_params = params;
  memcpy(&task.files, &iotask, sizeof(AsyncIOInput));

  gen_paths_worker_seq(gen_path_wrkr, &iodata, &task);
}

void generate_paths(CorrectAlnInput *tasks, size_t num_inputs,
                    GenPathWorker *workers, size_t num_workers)
{
  size_t i, read_counter = 0;

  for(i = 0; i < num_workers; i++)
    workers[i].shared_nreads = &read_counter;

  AsyncIOInput *asyncio_tasks = ctx_malloc(num_inputs * sizeof(AsyncIOInput));
  correct_aln_input_to_asycio(asyncio_tasks, tasks, num_inputs);

  asyncio_run_pool(asyncio_tasks, num_inputs, generate_paths_worker,
                   workers, num_workers, sizeof(GenPathWorker));

  ctx_free(asyncio_tasks);

  // Merge stats into workers[0]
  for(i = 1; i < num_workers; i++)
    correct_aln_merge_stats(&workers[0].corrector, &workers[i].corrector);
}
