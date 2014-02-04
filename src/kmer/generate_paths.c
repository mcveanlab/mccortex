#include "global.h"

#include <pthread.h>
#include "msgpool.h" // pool for getting jobs

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
  AsyncIOData data; // current data
  CorrectAlnReadsTask task; // current task
  LoadingStats stats;

  dBAlignment aln;
  CorrectAlnWorker corrector;

  // Nucleotides and positions of junctions
  Nucleotide *nuc_fw, *nuc_rv;
  size_t *pos_fw, *pos_rv;
  size_t num_fw, num_rv, junc_arrsize;

  // Packed representation of path
  uint8_t *packed;
  size_t packed_memcap;
};

#define INIT_BUFLEN 1024

// Should we print all paths?
boolean gen_paths_print_contigs = false;
volatile size_t print_contig_id = 0;


size_t gen_paths_worker_est_mem(const dBGraph *db_graph)
{
  size_t job_mem, corrector_mem, junc_mem, packed_mem;
  job_mem = 1024*4; // Assume 1024 bytes per read, 2 reads, seq+qual
  corrector_mem = correct_aln_worker_est_mem(db_graph);
  junc_mem = 2 * INIT_BUFLEN * (sizeof(Nucleotide)+sizeof(size_t));
  packed_mem = INIT_BUFLEN;

  return job_mem + corrector_mem + junc_mem + packed_mem + sizeof(GenPathWorker);
}

static void _gen_paths_worker_alloc(GenPathWorker *wrkr, dBGraph *db_graph)
{
  GenPathWorker tmp = {.db_graph = db_graph, .pool = NULL};

  asynciodata_alloc(&tmp.data);
  db_alignment_alloc(&tmp.aln);
  correct_aln_worker_alloc(&tmp.corrector, db_graph);
  loading_stats_init(&tmp.stats);
  // Junction data
  tmp.junc_arrsize = INIT_BUFLEN;
  tmp.nuc_fw = malloc2(tmp.junc_arrsize * sizeof(*tmp.nuc_fw));
  tmp.nuc_rv = malloc2(tmp.junc_arrsize * sizeof(*tmp.nuc_rv));
  tmp.pos_fw = malloc2(tmp.junc_arrsize * sizeof(*tmp.pos_fw));
  tmp.pos_rv = malloc2(tmp.junc_arrsize * sizeof(*tmp.pos_rv));
  tmp.num_fw = tmp.num_rv = 0;
  // Packed path temporary space
  tmp.packed_memcap = INIT_BUFLEN;
  tmp.packed = malloc2(tmp.packed_memcap);

  memcpy(wrkr, &tmp, sizeof(GenPathWorker));
}

static void _gen_paths_worker_dealloc(GenPathWorker *wrkr)
{
  asynciodata_dealloc(&wrkr->data);
  db_alignment_dealloc(&wrkr->aln);
  correct_aln_worker_dealloc(&wrkr->corrector);
  free(wrkr->nuc_fw); free(wrkr->nuc_rv);
  free(wrkr->pos_fw); free(wrkr->pos_rv);
  free(wrkr->packed);
}


GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph)
{
  size_t i;
  GenPathWorker *workers = malloc2(n * sizeof(GenPathWorker));
  for(i = 0; i < n; i++) _gen_paths_worker_alloc(&workers[i], graph);
  return workers;
}

void gen_paths_workers_dealloc(GenPathWorker *workers, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) _gen_paths_worker_dealloc(&workers[i]);
  free(workers);
}

static inline void worker_nuc_cap(GenPathWorker *wrkr, size_t cap)
{
  if(cap > wrkr->junc_arrsize) {
    wrkr->junc_arrsize = roundup2pow(cap);
    wrkr->nuc_fw = realloc2(wrkr->nuc_fw, wrkr->junc_arrsize * sizeof(Nucleotide));
    wrkr->nuc_rv = realloc2(wrkr->nuc_rv, wrkr->junc_arrsize * sizeof(Nucleotide));
    wrkr->pos_fw = realloc2(wrkr->pos_fw, wrkr->junc_arrsize * sizeof(size_t));
    wrkr->pos_rv = realloc2(wrkr->pos_rv, wrkr->junc_arrsize * sizeof(size_t));
  }
}

static inline void worker_packed_cap(GenPathWorker *wrkr, size_t nbytes)
{
  if(nbytes > wrkr->packed_memcap) {
    wrkr->packed_memcap = roundup2pow(nbytes);
    wrkr->packed = realloc2(wrkr->packed, wrkr->packed_memcap * sizeof(Nucleotide));
  }
}

// Merge stats into workers[0]
void generate_paths_merge_stats(GenPathWorker *wrkrs, size_t num_workers)
{
  size_t i;
  for(i = 1; i < num_workers; i++)
    correct_alignment_merge_hists(&wrkrs[0].corrector, &wrkrs[i].corrector);
}


// assume nbases > 0
#define bases_in_top_byte(nbases) ((((nbases) - 1) & 3) + 1)
#define bits_in_top_byte(nbases) (bases_in_top_byte(nbases) * 2)

// Returns number of paths added
static inline size_t _juncs_to_paths(const size_t *restrict pos_pl,
                                     const size_t *restrict pos_mn,
                                     size_t num_pl, size_t num_mn,
                                     const Nucleotide *nuc_pl,
                                     const boolean pl_is_fw,
                                     const dBNode *nodes,
                                     GenPathWorker *wrkr)
{
  size_t i, num_added = 0;
  size_t ctxcol = wrkr->task.crt_params.ctxcol;
  size_t ctpcol = wrkr->task.crt_params.ctpcol;

  dBNode node;
  size_t start_mn, start_pl, pos;
  PathLen plen, plen_orient;
  boolean added;
  PathIndex pindex; // address of path once added

  #ifdef CTXVERBOSE
    char str[num_pl+1];
    for(i = 0; i < num_pl; i++) str[i] = dna_nuc_to_char(nuc_pl[i]);
    str[num_pl] = '\0';
    status("[addpath] %s %s", pl_is_fw ? "fw" : "rv", str);
  #endif

  dBGraph *db_graph = wrkr->db_graph;

  // Create packed path with four diff offsets (0..3), point to correct one
  size_t pckd_memsize = sizeof(PathLen) + (num_pl+3)/4;
  uint8_t *packed_ptrs[4], *packed_ptr;

  worker_packed_cap(wrkr, pckd_memsize*4);
  for(i = 0; i < 4; i++) packed_ptrs[i] = wrkr->packed + i*pckd_memsize;

  packed_ptr = packed_ptrs[0]+sizeof(PathLen);
  pack_bases(packed_ptr, nuc_pl, num_pl);
  packed_cpy(packed_ptrs[1]+sizeof(PathLen), packed_ptr, 1, num_pl);
  packed_cpy(packed_ptrs[2]+sizeof(PathLen), packed_ptr, 2, num_pl);
  packed_cpy(packed_ptrs[3]+sizeof(PathLen), packed_ptr, 3, num_pl);

  // pl => plus in direction
  // mn => minus against direction

  // char nstr[num_pl+1];
  // for(i = 0; i < num_pl; i++) nstr[i] = dna_nuc_to_char(nuc_pl[i]);
  // nstr[num_pl] = '\0';
  // status("_juncs_to_paths(): len %zu %s", num_pl, nstr);

  // PathLen tmplen = num_pl;
  // for(i = 0; i < 4; i++) {
  //   memcpy(packed_ptrs[i], &tmplen, sizeof(PathLen));
  //   print_path(i, packed_ptrs[i], &db_graph->pdata);
  // }

  // 0,1,2,3
  // 3,2,1

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
    start_pl -= (start_pl > 0 && pos_pl[start_pl-1] == pos);

    // Check path is not too long (MAX_PATHLEN is the limit)
    plen = (PathLen)MIN2(num_pl - start_pl, MAX_PATHLEN);
    node = pl_is_fw ? nodes[pos] : db_node_reverse(nodes[pos]);

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

    // debug: Check path before adding
    graph_path_check_valid(node, ctxcol, packed_ptr+sizeof(PathLen), plen,
                           db_graph);

    added = graph_paths_find_or_add_mt(node.key, db_graph, ctpcol,
                                      packed_ptr, plen, &pindex);
    packed_ptr[top_idx] = top_byte; // restore top byte

    #ifdef CTXVERBOSE
      printf("We %s\n", added ? "added" : "abandoned");
    #endif

    // If the path already exists, all of its subpaths also already exist
    if(!added && plen < MAX_PATHLEN) break;
    num_added++;

    // debug: check path was added correctly
    Colour cols[2] = {ctxcol, ctpcol};
    GraphPathPairing gp = {.ctxcols = cols, .ctpcols = cols+1, .n = 1};
    graph_path_check_path(node.key, pindex, &gp, db_graph);
  }

  return num_added;
}

#undef bits_in_top_byte
#undef bases_in_top_byte

static void worker_junctions_to_paths(GenPathWorker *wrkr,
                                      const dBNodeBuffer *contig)
{
  size_t num_fw = wrkr->num_fw, num_rv = wrkr->num_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide *nuc_fw = wrkr->nuc_fw, *nuc_rv = wrkr->nuc_rv;

  worker_packed_cap(wrkr, (MAX2(num_fw, num_rv)+3)/4);

  assert(num_fw > 0 && num_rv > 0);

  #ifdef CTXVERBOSE
    status("num_fw: %zu num_rv: %zu", num_fw, num_rv);
  #endif

  // Reverse pos_rv, nuc_rv
  size_t i, j, tmp_pos;
  Nucleotide tmp_nuc;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp_pos);
    SWAP(nuc_rv[i], nuc_rv[j], tmp_nuc);
  }

  // _juncs_to_paths returns the number of paths added
  size_t n;
  const dBNode *nodes = contig->data;

  n = _juncs_to_paths(pos_fw, pos_rv, num_fw, num_rv, nuc_fw, true, nodes, wrkr);
  if(n) _juncs_to_paths(pos_rv, pos_fw, num_rv, num_fw, nuc_rv, false, nodes, wrkr);
}

static void worker_contig_to_junctions(GenPathWorker *wrkr,
                                       const dBNodeBuffer *contig)
{
  // status("nodebuf: %zu", wrkr->contig.len+MAX_KMER_SIZE+1);
  worker_nuc_cap(wrkr, contig->len);

  // gen_paths_print_contigs = true;
  if(gen_paths_print_contigs) {
    pthread_mutex_lock(&biglock);
    printf(">contig%zu\n", print_contig_id++);
    db_nodes_print(contig->data, contig->len, wrkr->db_graph, stdout);
    printf("\n");
    pthread_mutex_unlock(&biglock);
  }

  // Find forks in this colour
  Edges edges;
  size_t  i, num_fw = 0, num_rv = 0;
  int indegree, outdegree;
  Nucleotide *nuc_fw = wrkr->nuc_fw, *nuc_rv = wrkr->nuc_rv;
  size_t *pos_fw = wrkr->pos_fw, *pos_rv = wrkr->pos_rv;
  Nucleotide nuc;
  BinaryKmer bkmer;

  dBGraph *db_graph = wrkr->db_graph;
  const size_t kmer_size = wrkr->db_graph->kmer_size;

  const dBNode *nodes = contig->data;
  const size_t contig_len = contig->len;
  const size_t ctxcol = wrkr->task.crt_params.ctxcol;

  for(i = 0; i < contig_len; i++)
  {
    edges = db_node_get_edges(db_graph, ctxcol, nodes[i].key);
    outdegree = edges_get_outdegree(edges, nodes[i].orient);
    indegree = edges_get_indegree(edges, nodes[i].orient);

    if(indegree > 1 && i > 0)
    {
      bkmer = db_node_get_bkmer(db_graph, nodes[i-1].key);
      nuc = nodes[i-1].orient == FORWARD
              ? dna_nuc_complement(binary_kmer_first_nuc(bkmer, kmer_size))
              : binary_kmer_last_nuc(bkmer);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
    }

    // Only adding forward junctions after num_rv > 0 seems like a nice
    // optimisation but is actually a bad idea, in the case of outdegree > 1
    // just before first indegree >1, we need that junction
    if(outdegree > 1 && i+1 < contig_len)
    {
      bkmer = db_node_get_bkmer(db_graph, nodes[i+1].key);
      nuc = db_node_get_last_nuc(bkmer, nodes[i+1].orient, kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
    }
  }

  wrkr->num_fw = num_fw;
  wrkr->num_rv = num_rv;

  if(num_fw > 0 && num_rv > 0)
    worker_junctions_to_paths(wrkr, contig);
}

static void reads_to_paths(GenPathWorker *wrkr)
{
  AsyncIOData *data = &wrkr->data;
  read_t *r1 = &data->r1, *r2 = data->r2.seq.end > 0 ? &data->r2 : NULL;

  // if(r1->seq.end < 100) {
  //   printf("1>%s\n", r1->seq.b);
  //   if(r2 != NULL) printf("2>%s\n", r2->seq.b);
  // }

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
                          wrkr->db_graph);

  // For debugging
  // db_alignment_print(&wrkr->aln, wrkr->db_graph);

  // Correct sequence errors in the alignment
  correct_alignment_init(&wrkr->corrector, &wrkr->aln, wrkr->task.crt_params);

  dBNodeBuffer *nbuf;
  while((nbuf = correct_alignment_nxt(&wrkr->corrector)) != NULL)
    worker_contig_to_junctions(wrkr, nbuf);
}

// pthread method, loop: grabs job, does processing
static void* generate_paths_worker(void *ptr)
{
  GenPathWorker *wrkr = (GenPathWorker*)ptr;

  AsyncIOData data;
  while(msgpool_read(wrkr->pool, &data, &wrkr->data)) {
    // status("read: %s %s", data.r1.name.b, data.r2.name.b);
    wrkr->data = data;
    memcpy(&wrkr->task, data.ptr, sizeof(CorrectAlnReadsTask));
    reads_to_paths(wrkr);
  }

  pthread_exit(NULL);
}

void gen_path_worker_seq(GenPathWorker *wrkr, const CorrectAlnReadsTask *task,
                         const char *seq, size_t len)
{
  read_t *r1 = &wrkr->data.r1, *r2 = &wrkr->data.r2;
  seq_read_reset(r1);
  seq_read_reset(r2);

  // Copy seq to read1
  buffer_ensure_capacity(&r1->seq, len);
  memcpy(r1->seq.b, seq, len);
  r1->seq.b[len] = '\0';
  r1->seq.end = len;

  // Reset asyncio input data
  wrkr->data.fq_offset1 = wrkr->data.fq_offset2 = 0;
  wrkr->data.ptr = NULL;

  // Copy task to worker
  memcpy(&wrkr->task, task, sizeof(CorrectAlnReadsTask));

  reads_to_paths(wrkr);
}

void generate_paths(CorrectAlnReadsTask *tasks, size_t num_inputs,
                    GenPathWorker *workers, size_t num_workers)
{
  size_t i;
  MsgPool pool;
  msgpool_alloc_spinlock(&pool, MSGPOOLRSIZE, sizeof(AsyncIOData));
  msgpool_iterate(&pool, asynciodata_pool_init, NULL);

  for(i = 0; i < num_workers; i++)
    workers[i].pool = &pool;

  AsyncIOReadTask *asyncio_tasks = malloc2(num_inputs * sizeof(AsyncIOReadTask));
  correct_reads_input_to_asycio(asyncio_tasks, tasks, num_inputs);

  asyncio_run_threads(&pool, asyncio_tasks, num_inputs,
                      generate_paths_worker, workers, num_workers);

  free(asyncio_tasks);
  msgpool_iterate(&pool, asynciodata_pool_destroy, NULL);
  msgpool_dealloc(&pool);

  // Merge gap counts into worker[0]
  generate_paths_merge_stats(workers, num_workers);
}


// Save gap size distribution
// base_fmt is the beginning of the file name - the reset is <num>.csv or something
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void gen_paths_dump_gap_sizes(const char *base_fmt,
                              const uint64_t *arr, size_t arrlen,
                              size_t kmer_size, boolean insert_sizes,
                              size_t nreads)
{
  assert(arrlen > 0);

  // Print summary statistics: min, mean, median, mode, max
  size_t i, min, max, total, ngaps = 0, mode = 0;
  max = total = arr[0];

  for(min = 0; min < arrlen && arr[min] == 0; min++) {}

  if(min == arrlen) {
    if(insert_sizes) status("No insert gaps traversed");
    else status("No seq error gaps traversed");
    return;
  }

  for(i = 1; i < arrlen; i++) {
    if(arr[i] > 0) max = i;
    if(arr[i] > arr[mode]) mode = i;
    ngaps += arr[i];
    total += arr[i] * i;
  }

  double mean = (double)total / ngaps;
  float median = find_hist_median(arr, arrlen, ngaps);

  size_t ninputs = insert_sizes ? nreads/2 : nreads;
  char ngaps_str[100], ninputs_str[100];
  ulong_to_str(ngaps, ngaps_str);
  ulong_to_str(ninputs, ninputs_str);

  status("%s size distribution: "
         "min: %zu mean: %.1f median: %.1f mode: %zu max: %zu",
         insert_sizes ? "Insert" : "Seq error gap",
         min, mean, median, mode, max);

  status("  Gaps per read%s: %s / %s [%.3f]",
         insert_sizes ? " pair" : "", ngaps_str, ninputs_str,
         (double)ngaps / ninputs);

  StrBuf *csv_dump = strbuf_new();
  FILE *fout;

  if(!futil_generate_filename(base_fmt, csv_dump)) {
    warn("Cannot dump gapsize");
    strbuf_free(csv_dump);
    return;
  }

  if((fout = fopen(csv_dump->buff, "w")) == NULL) {
    warn("Cannot dump gapsize [cannot open: %s]", csv_dump->buff);
    strbuf_free(csv_dump);
    return;
  }

  fprintf(fout, "gap_in_kmers\tbp\tcount\n");

  if(arrlen > 0)
  {
    size_t start = 0, end = arrlen-1;

    while(start < arrlen && arr[start] == 0) start++;
    while(end > start && arr[end] == 0) end--;

    for(i = start; i <= end; i++) {
      fprintf(fout, "%4zu\t%4li\t%4zu\n",
              i, (long)i-(long)kmer_size, (size_t)arr[i]);
    }
  }

  status("Contig %s sizes dumped to %s\n",
         insert_sizes ? "insert" : "gap", csv_dump->buff);

  fclose(fout);
  strbuf_free(csv_dump);
}

// Get histogram array
const uint64_t* gen_paths_get_ins_gap(GenPathWorker *worker, size_t *len)
{
  return correct_alignment_get_inshist(&worker->corrector, len);
}

const uint64_t* gen_paths_get_err_gap(GenPathWorker *worker, size_t *len)
{
  return correct_alignment_get_errhist(&worker->corrector, len);
}

void gen_paths_get_stats(const GenPathWorker *worker, size_t num_workers,
                         LoadingStats *stats)
{
  size_t i;
  for(i = 0; i < num_workers; i++)
    loading_stats_merge(stats, &worker[i].stats);
}
