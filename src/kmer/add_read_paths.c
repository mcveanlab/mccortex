#include "global.h"

#include <pthread.h>
#include <unistd.h> // usleep

#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "util.h"
#include "file_util.h"
#include "seq_reader.h"
#include "file_reader.h"
#include "add_read_paths.h"

#include "repeat_walker.h"

// #define CTXVERBOSE 1

// Don't store assembly info longer than 1000 junctions
#define MAX_PATH 100

typedef struct {
  read_t r1, r2;
  Colour ctp_col, ctx_col;
  size_t gap_limit;
  int qcutoff1, qcutoff2;
  int hp_cutoff;
} AddPathsJob;

// Temp data store for adding paths to graph
typedef struct
{
  const size_t wid; // Worker id
  dBGraph *const db_graph;
  GraphWalker wlk; // traversing covg gaps
  RepeatWalker rptwlk; // traversing covg gaps
  dBNodeBuffer nodebuf; // Contigs constructed here
  AddPathsJob job; // Incoming read
  // Biggest gap between reads we'll try to traverse
  uint64_t *insert_sizes, *gap_sizes; // length gap_limit+1
} AddPathsWorker;

struct PathsWorkerPool
{
  // Threads
  size_t num_of_threads;
  pthread_t *threads;
  AddPathsWorker *workers;

  // Shared data
  dBGraph *db_graph;
  boolean seen_pe;

  // Data currently being read
  read_t r1, r2;
  Colour ctp_col, ctx_col;
  size_t gap_limit;
};

//
// Multithreading
//
pthread_attr_t thread_attr;
pthread_mutex_t add_paths_mutex;
pthread_mutex_t reader_mutex;
pthread_mutex_t data_written_mutex, data_read_mutex;
pthread_cond_t data_written_cond, data_read_cond;

// Incoming data
volatile boolean data_waiting = false, input_ended = false;
volatile AddPathsJob next_job;

static void paths_worker_alloc(AddPathsWorker *worker, size_t wid,
                               size_t gap_limit, dBGraph *db_graph)
{
  AddPathsWorker tmp = {.wid = wid, .db_graph = db_graph};
  memcpy(worker, &tmp, sizeof(AddPathsWorker));

  db_node_buf_alloc(&worker->nodebuf, 4096);
  graph_walker_alloc(&worker->wlk);
  walker_alloc(&worker->rptwlk, db_graph->ht.capacity, 22); // 4MB

  if(!seq_read_alloc(&worker->job.r1) || !seq_read_alloc(&worker->job.r2))
    die("Out of memory");

  worker->insert_sizes = calloc2(2*(gap_limit+1), sizeof(uint64_t));
  worker->gap_sizes = worker->insert_sizes + gap_limit+1;
}

static void paths_worker_dealloc(AddPathsWorker *worker)
{
  free(worker->insert_sizes);
  walker_dealloc(&worker->rptwlk);
  graph_walker_dealloc(&worker->wlk);
  db_node_buf_dealloc(&worker->nodebuf);
  seq_read_dealloc(&worker->job.r1);
  seq_read_dealloc(&worker->job.r2);
}

static void dump_gap_sizes(const char *base_fmt, const uint64_t *arr,
                           size_t arrlen, size_t kmer_size)
{
  StrBuf *csv_dump = strbuf_new();

  if(!futil_generate_filename(base_fmt, csv_dump)) {
    warn("Cannot dump gapsize");
    return;
  }

  FILE *fh;

  if((fh = fopen(csv_dump->buff, "w")) == NULL) {
    warn("Cannot dump gapsize [cannot open: %s]", csv_dump->buff);
    strbuf_free(csv_dump);
    return;
  }

  fprintf(fh, "gap_in_kmers,bp,count\n");

  if(arrlen > 0)
  {
    size_t i, start = 0, end = arrlen-1;

    while(start < arrlen && arr[start] == 0) start++;
    while(end > start && arr[end] == 0) end--;

    for(i = start; i <= end; i++) {
      fprintf(fh, "%4zu,%4li,%4zu\n", i, (long)i-kmer_size, (size_t)arr[i]);
    }
  }

  status("Contig gap sizes dumped to %s\n", csv_dump->buff);

  fclose(fh);
  strbuf_free(csv_dump);
}

static void construct_paths(Nucleotide *nuc_fw, size_t *pos_fw, size_t num_fw,
                            Nucleotide *nuc_rv_tmp, size_t *pos_rv_tmp, size_t num_rv,
                            const dBNode *nodes, dBGraph *db_graph,
                            Colour ctp_col)
{
  hkey_t node;
  Orientation orient;
  size_t start_fw, start_rv, pos;
  PathIndex prev_index, new_index;
  PathLen plen;
  Nucleotide *bases;

  PathStore *paths = &db_graph->pdata;

  // Reverse rv
  size_t i, j;
  size_t pos_rv[MAX_PATH];
  Nucleotide nuc_rv[MAX_PATH];

  for(i = 0, j = num_rv-1; i < num_rv; i++, j--) {
    pos_rv[j] = pos_rv_tmp[i];
    nuc_rv[j] = nuc_rv_tmp[i];
  }

  // Get Lock
  pthread_mutex_lock(&add_paths_mutex);

  //
  // Generate paths going backwards through the contig
  //
  #ifdef CTXVERBOSE
    printf("==REV==\n");
    char str[MAX_KMER_SIZE+1];
  #endif

  for(start_rv = 0, start_fw = num_fw-1; start_fw != SIZE_MAX; start_fw--)
  {
    while(start_rv < num_rv && pos_rv[start_rv] > pos_fw[start_fw]) start_rv++;
    if(start_rv == num_rv) break;

    pos = pos_fw[start_fw] + 1;
    start_rv -= (start_rv > 0 && pos_rv[start_rv-1] == pos);

    bases = nuc_rv + start_rv;
    plen = num_rv - start_rv;
    node = nodes[pos].key;
    orient = rev_orient(nodes[pos].orient);
    prev_index = db_node_paths(db_graph, node);

    #ifdef CTXVERBOSE
      binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, str);
      printf(" %s:%i) start_rv: %zu start_fw: %zu {%zu}\n", str, orient,
             start_rv, start_fw, pos_fw[start_fw]);
    #endif

    new_index = path_store_add(paths, prev_index, plen, bases, orient, ctp_col);
    if(new_index != PATH_NULL) db_node_paths(db_graph, node) = new_index;
  }

  //
  // Generate forward paths
  //
  #ifdef CTXVERBOSE
    printf("==FWD==\n");
  #endif

  for(start_fw = 0, start_rv = num_rv-1; start_rv != SIZE_MAX; start_rv--)
  {
    while(start_fw < num_fw && pos_fw[start_fw] < pos_rv[start_rv]) start_fw++;
    if(start_fw == num_fw) break;

    pos = pos_rv[start_rv] - 1;
    start_fw -= (start_fw > 0 && pos_fw[start_fw-1] == pos);

    bases = nuc_fw + start_fw;
    plen = num_fw - start_fw;
    node = nodes[pos].key;
    orient = nodes[pos].orient;
    prev_index = db_node_paths(db_graph, node);

    #ifdef CTXVERBOSE
      binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, str);
      printf(" %s:%i) start_rv: %zu start_fw: %zu\n", str, orient,
             start_rv, start_fw);
    #endif

    new_index = path_store_add(paths, prev_index, plen, bases, orient, ctp_col);
    if(new_index != PATH_NULL) db_node_paths(db_graph, node) = new_index;
  }

  // Free lock
  pthread_mutex_unlock(&add_paths_mutex);
}

static void add_read_path(const dBNode *nodes, size_t len,
                          dBGraph *graph, Colour ctx_col, Colour ctp_col)
{
  if(len < 3) return;

  #ifdef CTXVERBOSE
    printf("contig: ");
    db_nodes_print(nodes, len, graph, stdout);
    printf("\n");
  #endif

  // Find forks in this colour
  Edges edges;
  size_t  i, j, k, num_fw = 0, num_rv = 0, indegree, outdegree, addfw, addrv;
  Nucleotide nuc_fw[MAX_PATH], nuc_rv[MAX_PATH];
  size_t pos_fw[MAX_PATH], pos_rv[MAX_PATH];
  Nucleotide nuc;
  BinaryKmer bkmer;
  size_t last_fw_num = 0;

  for(i = 0; i < len; i++)
  {
    edges = db_node_edges(graph, ctx_col, nodes[i].key);
    outdegree = edges_get_outdegree(edges, nodes[i].orient);
    indegree = edges_get_indegree(edges, nodes[i].orient);

    addfw = (outdegree > 1 && i+1 < len);
    addrv = (indegree > 1 && i > 0);

    if(addfw) {
      bkmer = db_node_bkmer(graph, nodes[i+1].key);
      nuc = db_node_last_nuc(bkmer, nodes[i+1].orient, graph->kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
    }
    if(addrv) {
      bkmer = db_node_bkmer(graph, nodes[i-1].key);
      nuc = nodes[i-1].orient == FORWARD
              ? binary_nuc_complement(binary_kmer_first_nuc(bkmer, graph->kmer_size))
              : binary_kmer_last_nuc(bkmer);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
    }

    if(num_fw == MAX_PATH || num_rv == MAX_PATH)
    {
      size_t cutoff = i + 1 - MAX_PATH;
      for(j = 0; j < num_fw && pos_fw[j] <= cutoff; j++) {}
      for(k = 0; k < num_rv && pos_rv[k] <= cutoff; k++) {}
      if(k > 0 && num_fw > 0) {
        construct_paths(nuc_fw, pos_fw, num_fw, nuc_rv, pos_rv, num_rv,
                        nodes, graph, ctp_col);
        last_fw_num = num_fw-j;
      }
      if(j > 0) {
        num_fw -= j;
        memmove(nuc_fw, nuc_fw+j, num_fw * sizeof(*nuc_fw));
        memmove(pos_fw, pos_fw+j, num_fw * sizeof(*pos_fw));
      }
      if(k > 0) {
        num_rv -= k;
        memmove(nuc_rv, nuc_rv+k, num_rv * sizeof(*nuc_rv));
        memmove(pos_rv, pos_rv+k, num_rv * sizeof(*pos_rv));
      }
    }
  }

  if(num_fw > last_fw_num && num_rv > 0) {
    construct_paths(nuc_fw, pos_fw, num_fw, nuc_rv, pos_rv, num_rv,
                    nodes, graph, ctp_col);
  }
}

// fill in gap in read node1==>---<==node2
// Adds at most GAP_LIMIT nodes to nodebuf
// If successful: traverses gap and adds new nodes including node2/orient2
//    returns total number of nodes added (>= 0)
// If unsucessful: doesn't add anything to the nodebuf, returns -1
static int traverse_gap(dBNodeBuffer *nodebuf,
                        hkey_t node2, Orientation orient2,
                        size_t gap_limit,
                        Colour ctxcol, Colour ctpcol,
                        GraphWalker *wlk, RepeatWalker *rptwlk,
                        const dBGraph *db_graph)
{
  hkey_t node1 = nodebuf->data[nodebuf->len-1].key;
  Orientation orient1 = nodebuf->data[nodebuf->len-1].orient;

  #ifdef CTXVERBOSE
    char tmp1[MAX_KMER_SIZE+1], tmp2[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_bkmer(db_graph, node1), db_graph->kmer_size, tmp1);
    binary_kmer_to_str(db_node_bkmer(db_graph, node2), db_graph->kmer_size, tmp2);
    printf("traverse gap: %s:%i -> %s:%i\n", tmp1, orient1, tmp2, orient2);
  #endif

  // First and last node already match
  if(node1 == node2 && orient1 == orient2) return 0;

  // Ensure capacity
  db_node_buf_ensure_capacity(nodebuf, nodebuf->len + gap_limit);

  dBNode *nodes = nodebuf->data + nodebuf->len;
  size_t pos = 0;
  Nucleotide lost_nuc;

  // Walk from left -> right
  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node1, orient1);
  lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

  // need to call db_node_has_col only if more than one colour loaded
  while(pos < gap_limit && graph_traverse(wlk) &&
        walker_attempt_traverse(rptwlk, wlk, wlk->node, wlk->orient, wlk->bkmer))
  {
    graph_walker_node_add_counter_paths(wlk, lost_nuc);
    lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

    nodes[pos].key = wlk->node;
    nodes[pos].orient = wlk->orient;
    pos++;

    if(wlk->node == node2 && wlk->orient == orient2)
      break;
  }

  graph_walker_finish(wlk);
  walker_fast_clear(rptwlk, nodes, pos);

  if(wlk->node == node2 && wlk->orient == orient2) {
    nodebuf->len += pos;
    return pos;
  }

  // Walk from right -> left
  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node2, opposite_orientation(orient2));

  pos = gap_limit-1;
  nodes[pos].key = node2;
  nodes[pos].orient = orient2;
  // pos is now the index at which we last added a node

  Orientation orient;
  boolean success = false;
  lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

  // need to call db_node_has_col only if more than one colour loaded
  while(pos > 0 && graph_traverse(wlk) &&
        walker_attempt_traverse(rptwlk, wlk, wlk->node, wlk->orient, wlk->bkmer))
  {
    graph_walker_node_add_counter_paths(wlk, lost_nuc);
    orient = opposite_orientation(wlk->orient);

    if(wlk->node == node1 && orient == orient1) {
      success = true;
      break;
    }

    lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

    pos--;
    nodes[pos].key = wlk->node;
    nodes[pos].orient = orient;
  }

  graph_walker_finish(wlk);

  walker_fast_clear(rptwlk, nodes+pos, gap_limit-pos);

  if(success)
  {
    #ifdef CTXVERBOSE
      printf(" traverse success\n");
    #endif
    size_t num = gap_limit - pos;
    memmove(nodes, nodes + gap_limit - num, num * sizeof(dBNode));
    nodebuf->len += num;
    return num;
  }

  return -1;
}

void read_to_path(AddPathsWorker *worker)
{
  dBNodeBuffer *nodebuf = &worker->nodebuf;
  dBGraph *db_graph = worker->db_graph;
  AddPathsJob *job = &worker->job;
  const size_t ctx_col = job->ctx_col;

  size_t i, r2_start = 0;
  int r2_offset = -1;
  nodebuf->len = 0;

  get_nodes_from_read(&job->r1, job->qcutoff1, job->hp_cutoff, db_graph, nodebuf);

  if(job->r2.seq.end >= db_graph->kmer_size)
  {
    seq_read_reverse_complement(&job->r2);

    if(nodebuf->len > 0)
    {
      // Insert gap
      db_node_buf_ensure_capacity(nodebuf, nodebuf->len+1);
      nodebuf->data[nodebuf->len].key = HASH_NOT_FOUND;
      nodebuf->data[nodebuf->len].orient = FORWARD;
      nodebuf->len++;
    }

    r2_start = nodebuf->len;
    r2_offset = get_nodes_from_read(&job->r2, job->qcutoff2, job->hp_cutoff,
                                    db_graph, nodebuf);
  }

  if(nodebuf->len == 0) return;

  // Check for gaps
  for(i = 0; i < nodebuf->len && nodebuf->data[i].key != HASH_NOT_FOUND; i++) {}

  if(i == nodebuf->len)
  {
    // No gaps in contig
    add_read_path(nodebuf->data, nodebuf->len, db_graph,
                  ctx_col, job->ctp_col);
  }
  else
  {
    const size_t end = nodebuf->len;

    db_node_buf_ensure_capacity(nodebuf, nodebuf->len + nodebuf->len);

    hkey_t node, prev_node;
    Orientation orient;

    node = prev_node = nodebuf->data[0].key;
    assert(nodebuf->data[0].key != HASH_NOT_FOUND);

    for(i = 0; i < end; i++, prev_node = node)
    {
      node = nodebuf->data[i].key;
      orient = nodebuf->data[i].orient;

      #ifdef CTXVERBOSE
        char str[MAX_KMER_SIZE+1+7] = {0};
        if(node == HASH_NOT_FOUND) strcpy(str, "(none)");
        else {
          BinaryKmer bkmer = db_node_bkmer(db_graph, node);
          binary_kmer_to_str(bkmer, db_graph->kmer_size, str);
        }
        printf("  node:%zu %zu %s\n", i, (size_t)node, str);
      #endif

      if(node != HASH_NOT_FOUND)
      {
        if(prev_node == HASH_NOT_FOUND)
        {
          // Can we branch the gap from the prev contig?
          int gapsize = traverse_gap(nodebuf, node, orient, job->gap_limit,
                                     ctx_col, job->ctp_col,
                                     &worker->wlk, &worker->rptwlk, db_graph);

          if(gapsize == -1)
          {
            // Failed to bridge gap
            add_read_path(nodebuf->data+end, nodebuf->len-end, db_graph,
                          ctx_col, job->ctp_col);

            nodebuf->data[end].key = node;
            nodebuf->data[end].orient = orient;
            nodebuf->len = end+1;
          }
          else
          {
            // Update stats (gapsize is <= GAP_LIMIT)
            if(i == r2_start) {
              // Kmer is start of second read
              int gap = gapsize - r2_offset;
              worker->insert_sizes[gap < 0 ? 0 : gap]++;
            }
            else {
              worker->gap_sizes[gapsize]++;
            }
          }
        }
        else {
          nodebuf->data[nodebuf->len].key = node;
          nodebuf->data[nodebuf->len].orient = orient;
          nodebuf->len++;
        }
      }
    }

    add_read_path(nodebuf->data+end, nodebuf->len-end, db_graph,
                  ctx_col, job->ctp_col);
  }
}

void* add_paths_thread(void *ptr)
{
  AddPathsWorker *worker = (AddPathsWorker*)ptr;
  int got_job;

  while(data_waiting || !input_ended)
  {
    got_job = 0;

    pthread_mutex_lock(&reader_mutex);

    // wait until data are ready
    pthread_mutex_lock(&data_written_mutex);
    while(!data_waiting && !input_ended)
      pthread_cond_wait(&data_written_cond, &data_written_mutex);
    pthread_mutex_unlock(&data_written_mutex);

    if(data_waiting) {
      AddPathsJob tmp_job;
      SWAP(worker->job, next_job, tmp_job);
      data_waiting = false;
      got_job = 1;
    }

    pthread_mutex_lock(&data_read_mutex);
    pthread_cond_signal(&data_read_cond);
    pthread_mutex_unlock(&data_read_mutex);

    pthread_mutex_unlock(&reader_mutex);

    // Do work
    if(got_job) read_to_path(worker);
  }

  pthread_exit(NULL);
}

// This function is passed to parse_filelist to load paths from sequence data
static void load_paths(read_t *r1, read_t *r2,
                       int fq_offset1, int fq_offset2,
                       const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                       void *ptr)
{
  // Don't bother checking for duplicates
  (void)stats;

  PathsWorkerPool *pool = (PathsWorkerPool*)ptr;

  // printf("READ1: %s\n", r1->seq.b);
  // if(r2 != NULL) printf("READ2: %s\n", r2->seq.b);

  int qcutoff1 = prefs->quality_cutoff;
  int qcutoff2 = prefs->quality_cutoff;

  if(prefs->quality_cutoff > 0)
  {
    qcutoff1 += fq_offset1;
    qcutoff2 += fq_offset2;
  }

  pthread_mutex_lock(&data_read_mutex);
  while(data_waiting) pthread_cond_wait(&data_read_cond, &data_read_mutex);
  pthread_mutex_unlock(&data_read_mutex);

  next_job.qcutoff1 = qcutoff1;
  next_job.qcutoff2 = qcutoff2;
  next_job.hp_cutoff = prefs->homopolymer_cutoff;
  next_job.ctp_col = pool->ctp_col;
  next_job.ctx_col = pool->ctx_col;
  next_job.gap_limit = pool->gap_limit;

  read_t tmp_read;
  SWAP(*r1, next_job.r1, tmp_read);
  if(r2 != NULL) SWAP(*r2, next_job.r2, tmp_read);
  else {
    next_job.r2.seq.end = 0;
    next_job.r2.seq.b[0] = '\0';
  }

  data_waiting = true;

  // Pass to a thread
  pthread_mutex_lock(&data_written_mutex);
  pthread_cond_signal(&data_written_cond);
  pthread_mutex_unlock(&data_written_mutex);

  // Update stats
  stats->total_good_reads += 1 + (r2 != NULL);
}

PathsWorkerPool* paths_worker_pool_new(size_t num_of_threads,
                                       dBGraph *db_graph, size_t gap_limit)
                                       // FILE *output_handle)
{
  PathsWorkerPool *pool = calloc2(1, sizeof(PathsWorkerPool));

  size_t i;
  pool->workers = malloc2(sizeof(AddPathsWorker) * num_of_threads);
  pool->num_of_threads = num_of_threads;
  pool->db_graph = db_graph;
  pool->seen_pe = false;

  if(!seq_read_alloc(&pool->r1) || !seq_read_alloc(&pool->r2))
    die("Out of memory");

  for(i = 0; i < num_of_threads; i++) {
    paths_worker_alloc(pool->workers + i, i, gap_limit, db_graph);
  }

  // Start threads
  pool->threads = malloc2(num_of_threads * sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

  if(pthread_mutex_init(&add_paths_mutex, NULL) != 0) die("mutex init failed");
  if(pthread_mutex_init(&reader_mutex, NULL) != 0) die("mutex init failed");
  if(pthread_mutex_init(&data_written_mutex, NULL) != 0) die("mutex init failed");
  if(pthread_mutex_init(&data_read_mutex, NULL) != 0) die("mutex init failed");

  if(pthread_cond_init(&data_written_cond, NULL) != 0) die("pthread_cond init failed");
  if(pthread_cond_init(&data_read_cond, NULL) != 0) die("pthread_cond init failed");

  if(seq_read_alloc((read_t*)&next_job.r1) == NULL ||
     seq_read_alloc((read_t*)&next_job.r2) == NULL)
    die("Out of memory");

  data_waiting = input_ended = false;
  // fout = output_handle;

  int rc;
  for(i = 0; i < num_of_threads; i++)
  {
    rc = pthread_create(pool->threads+i, &thread_attr, add_paths_thread,
                        (void*)(pool->workers+i));
    if(rc != 0) die("Creating thread failed");
  }

  return pool;
}

void paths_worker_pool_dealloc(PathsWorkerPool *pool)
{
  input_ended = true;
  status("Waiting for threads to finish...");

  // Release waiting worker threads
  size_t i, j, gap_limit = pool->gap_limit;
  for(i = 0; i < pool->num_of_threads; i++) {
    pthread_mutex_lock(&data_written_mutex);
    pthread_cond_signal(&data_written_cond);
    pthread_mutex_unlock(&data_written_mutex);
  }

  // sleep(5);

  int rc;
  for(i = 0; i < pool->num_of_threads; i++) {
    rc = pthread_join(pool->threads[i], NULL);
    if(rc != 0) die("Joining thread failed");
  }

  // Sum insert sizes and gap sizes into worker 0
  uint64_t *insert_sizes = pool->workers[0].insert_sizes;
  uint64_t *gap_sizes = pool->workers[0].gap_sizes;

  for(i = 1; i < pool->num_of_threads; i++) {
    for(j = 0; j < gap_limit+1; j++) {
      insert_sizes[j] += pool->workers[i].insert_sizes[j];
      gap_sizes[j] += pool->workers[i].gap_sizes[j];
    }
    paths_worker_dealloc(&(pool->workers[i]));
  }

  // Print mp gap size / insert stats to a file
  size_t kmer_size = pool->db_graph->kmer_size;
  dump_gap_sizes("gap_sizes.%u.csv", gap_sizes, gap_limit+1, kmer_size);

  if(pool->seen_pe)
    dump_gap_sizes("mp_sizes.%u.csv", insert_sizes, gap_limit+1, kmer_size);

  paths_worker_dealloc(&(pool->workers[0]));

  free(pool->threads);
  free(pool->workers);

  seq_read_dealloc(&pool->r1);
  seq_read_dealloc(&pool->r2);

  seq_read_dealloc((read_t*)&next_job.r1);
  seq_read_dealloc((read_t*)&next_job.r2);

  pthread_cond_destroy(&data_written_cond);
  pthread_cond_destroy(&data_read_cond);

  pthread_mutex_destroy(&add_paths_mutex);
  pthread_mutex_destroy(&reader_mutex);
  pthread_mutex_destroy(&data_written_mutex);
  pthread_mutex_destroy(&data_read_mutex);
  pthread_attr_destroy(&thread_attr);

  free(pool);
}

void add_read_paths_to_graph(PathsWorkerPool *pool,
                             seq_file_t *sf1, seq_file_t *sf2,
                             size_t gap_limit, size_t ctx_col, size_t ctp_col,
                             const SeqLoadingPrefs *prefs,
                             SeqLoadingStats *stats)
{
  pool->gap_limit = gap_limit;
  pool->ctx_col = ctx_col;
  pool->ctp_col = ctp_col;

  if(sf2 == NULL)
    seq_parse_se_sf(sf1, &pool->r1, &pool->r2, prefs, stats, &load_paths, pool);
  else
    seq_parse_pe_sf(sf1, sf2, &pool->r1, &pool->r2, prefs, stats, &load_paths, pool);
}
