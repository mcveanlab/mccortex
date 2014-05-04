#include "global.h"
#include "breakpoint_caller.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "kmer_occur.h"
#include "graph_crawler.h"

// Require 5 kmers on the reference before and after breakpoint
#define REQ_REF_NKMERS 5

typedef struct
{
  // Specific to this instance
  const size_t threadid, nthreads;

  // Temporary memory used by this instance
  GraphCrawler crawlers[2]; // [0] => FORWARD, [1] => REVERSE
  dBNodeBuffer nbuf;

  // 5p flank
  KOccurRunBuffer korun_buf;
  dBNodeBuffer flank5pbuf;

  // Passed to all instances
  const KOGraph kograph;
  const dBGraph *db_graph;
  gzFile gzout;
  pthread_mutex_t *const out_lock;
  size_t *callid;
} BreakpointCaller;

static BreakpointCaller* brkpt_callers_new(size_t num_callers, gzFile gzout,
                                           const KOGraph kograph,
                                           const dBGraph *db_graph)
{
  BreakpointCaller *callers = ctx_malloc(num_callers * sizeof(BreakpointCaller));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  size_t *callid = ctx_calloc(1, sizeof(size_t));

  size_t i;
  for(i = 0; i < num_callers; i++)
  {
    BreakpointCaller tmp = {.threadid = i,
                            .nthreads = num_callers,
                            .kograph = kograph,
                            .db_graph = db_graph,
                            .gzout = gzout,
                            .out_lock = out_lock,
                            .callid = callid};

    memcpy(&callers[i], &tmp, sizeof(BreakpointCaller));

    db_node_buf_alloc(&callers[i].nbuf, 1024);
    db_node_buf_alloc(&callers[i].flank5pbuf, 1024);
    kmer_run_buf_alloc(&callers[i].korun_buf, 512);
    graph_crawler_alloc(&callers[i].crawlers[0], db_graph);
    graph_crawler_alloc(&callers[i].crawlers[1], db_graph);
  }

  return callers;
}

static void brkpt_callers_destroy(BreakpointCaller *callers, size_t num_callers)
{
  size_t i;
  for(i = 0; i < num_callers; i++) {
    db_node_buf_dealloc(&callers[i].nbuf);
    db_node_buf_dealloc(&callers[i].flank5pbuf);
    kmer_run_buf_dealloc(&callers[i].korun_buf);
    graph_crawler_dealloc(&callers[i].crawlers[0]);
    graph_crawler_dealloc(&callers[i].crawlers[1]);
  }
  pthread_mutex_destroy(callers[0].out_lock);
  ctx_free(callers[0].out_lock);
  ctx_free(callers[0].callid);
  ctx_free(callers);
}

static inline void ko_gzprint(gzFile gzout, size_t kmer_size,
                              KOGraph kograph, KOccurRun ko)
{
  const char strand[] = {'+','-'};
  const char *chrom = kograph_chrom(kograph,ko).name;
  // get end coord as inclusive coord
  size_t start, end;
  if(ko.next > ko.start) {
    start = ko.start;
    end = ko.next-1+kmer_size-1;
  } else {
    start = ko.start+kmer_size-1;
    end = ko.next+1;
  }
  // +1 to coords to convert to 1-based
  gzprintf(gzout, "%s:%zu-%zu:%c", chrom, start+1, end+1, strand[ko.strand]);
}

static void process_contig(BreakpointCaller *caller,
                           dBNodeBuffer *flank5p,
                           const KOccurRunBuffer *korun_buf,
                           dBNodeBuffer *nbuf,
                           const uint32_t *cols, size_t ncols)
{
  gzFile gzout = caller->gzout;
  KOGraph kograph = caller->kograph;

  // Find first place we meet the ref
  // DEV: add requirement for REQ_REF_NKMERS kmers
  size_t i, end;
  const size_t kmer_size = caller->db_graph->kmer_size;

  for(i = 0; i < nbuf->len; i++) {
    if(kograph_num(caller->kograph, nbuf->data[i].key) > 0) break;
  }
  end = i;

  if(end == nbuf->len) {
    // status("Never met ref");
    return; // we never met the ref again
  }

  size_t callid = __sync_fetch_and_add((volatile size_t*)caller->callid, 1);

  pthread_mutex_lock(caller->out_lock);

  gzprintf(gzout, ">call.%zu.5pflank cols=%zu", callid, cols[0]);
  for(i = 1; i < ncols; i++) gzprintf(gzout, ",%zu", cols[i]);

  // Print kmer occurances
  gzprintf(gzout, " chr=");
  for(i = 0; i < korun_buf->len; i++) {
    if(i > 0) gzputc(gzout, ',');
    ko_gzprint(gzout, kmer_size, kograph, korun_buf->data[i]);
  }

  gzputc(gzout, '\n');
  db_nodes_gzprint(flank5p->data, flank5p->len, caller->db_graph, gzout);
  gzputc(gzout, '\n');

  gzprintf(gzout, ">call.%zu.3pflank\n", callid);
  db_nodes_gzprint_cont(nbuf->data+end, nbuf->len-end, caller->db_graph, gzout);
  gzputc(gzout, '\n');

  gzprintf(gzout, ">call.%zu.path\n", callid);
  db_nodes_gzprint_cont(nbuf->data, end, caller->db_graph, gzout);
  gzprintf(gzout, "\n\n");

  pthread_mutex_unlock(caller->out_lock);
}

// Traverse backwards from node1 -> node0
static void traverse_5pflank(GraphCrawler *crawler, dBNode node0, dBNode node1)
{
  const dBGraph *db_graph = crawler->cache.db_graph;
  dBNode nodes[4];
  Nucleotide bases[4];
  size_t i, num_next;
  BinaryKmer bkmer0 = db_node_get_bkmer(db_graph, node0.key);

  num_next = db_graph_next_nodes(db_graph, bkmer0, node0.orient,
                                 db_node_edges(db_graph, node0.key, 0),
                                 nodes, bases);

  for(i = 0; i < num_next && !db_nodes_are_equal(nodes[i],node1); i++) {}

  // Go backwards to get 5p flank
  // NULL means loop from 0..(ncols-1)
  graph_crawler_fetch(crawler, node0, nodes, bases, i, num_next,
                      REQ_REF_NKMERS, NULL, db_graph->num_of_cols);
}

// Walk the graph remembering the last time we met the ref
// When traversal fails, dump sequence up to last meeting with the ref
static void follow_break(BreakpointCaller *caller, dBNode node)
{
  size_t i, j, k, num_next;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  const dBGraph *db_graph = caller->db_graph;

  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);
  Edges edges = db_node_get_edges(db_graph, node.key, 0);

  num_next = db_graph_next_nodes(db_graph, bkey, node.orient, edges,
                                 next_nodes, next_nucs);

  // Filter out next nodes in the reference
  for(i = 0, j = 0; i < num_next; i++) {
    if(kograph_num(caller->kograph, next_nodes[i].key) == 0) {
      next_nodes[j] = next_nodes[i];
      next_nucs[j] = next_nucs[i];
      j++;
    }
  }

  // Abandon if all options are in ref or none are
  if(j == num_next || j == 0) return;
  num_next = j;

  // Follow all paths not in ref, in all colours
  GraphCrawler *fw_crawler = &caller->crawlers[node.orient];
  GraphCrawler *rv_crawler = &caller->crawlers[!node.orient];
  dBNodeBuffer *nbuf = &caller->nbuf, *flank5pbuf = &caller->flank5pbuf;
  KOccurRunBuffer *korun_buf = &caller->korun_buf; // Kmer occurance runs

  for(i = 0; i < num_next; i++)
  {
    // Go backwards to get 5p flank
    traverse_5pflank(rv_crawler, db_node_reverse(next_nodes[i]),
                     db_node_reverse(node));

    // Loop over the flanks we got
    for(j = 0; j < rv_crawler->num_paths; j++)
    {
      // Get 5p flank
      db_node_buf_reset(flank5pbuf);
      graph_crawler_get_path_nodes(rv_crawler, j, flank5pbuf);
      db_nodes_reverse_complement(flank5pbuf->data, flank5pbuf->len);

      // Check this flank is in the ref
      kograph_filter_stretch(caller->kograph, flank5pbuf->data, flank5pbuf->len,
                             korun_buf);

      if(korun_buf->len > 0)
      {
        // Only traverse in the colours we have a flank for
        GCMultiColPath *flankpath = &rv_crawler->multicol_paths[j];
        graph_crawler_fetch(fw_crawler, node,
                            next_nodes, next_nucs, j, num_next, -1,
                            flankpath->cols, flankpath->num_cols);

        // Assemble contigs - fetch forwards for each 5p flank
        for(k = 0; k < fw_crawler->num_paths; k++) {
          db_node_buf_reset(nbuf);
          graph_crawler_get_path_nodes(fw_crawler, k, nbuf);
          GCMultiColPath *multicol_path = &fw_crawler->multicol_paths[k];
          process_contig(caller, flank5pbuf, korun_buf,
                         nbuf, multicol_path->cols, multicol_path->num_cols);
        }
      }
    }
  }
}

void breakpoint_caller_node(hkey_t hkey, BreakpointCaller *caller)
{
  Edges edges;

  graph_crawler_reset(&caller->crawlers[0]);
  graph_crawler_reset(&caller->crawlers[1]);

  // check node is in the ref
  if(kograph_num(caller->kograph, hkey) > 0) {
    edges = db_node_get_edges(caller->db_graph, hkey, 0);
    if(edges_get_outdegree(edges, FORWARD) > 1) {
      follow_break(caller, (dBNode){.key = hkey, .orient = FORWARD});
    }
    if(edges_get_outdegree(edges, REVERSE) > 1) {
      follow_break(caller, (dBNode){.key = hkey, .orient = REVERSE});
    }
  }
}

static void breakpoint_caller(void *ptr)
{
  BreakpointCaller *caller = (BreakpointCaller*)ptr;
  ctx_assert(caller->db_graph->num_edge_cols == 1);

  HASH_ITERATE_PART(&caller->db_graph->ht, caller->threadid, caller->nthreads,
                    breakpoint_caller_node, caller);
}

static void breakpoints_print_header(gzFile gzout, const CmdArgs *args,
                                     const char **seq_paths, size_t nseq_paths)
{
  char datestr[9], cwd[PATH_MAX + 1];
  size_t i;

  ctx_assert(nseq_paths > 0);

  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  gzprintf(gzout, "##fileFormat=CtxBreakpointsv0.1\n##fileDate=%s\n", datestr);
  gzprintf(gzout, "##cmd=\"%s\"\n", args->cmdline);
  if(futil_get_current_dir(cwd) != NULL) gzprintf(gzout, "##wkdir=%s\n", cwd);

  gzprintf(gzout, "##reference=%s", seq_paths[0]);
  for(i = 1; i < nseq_paths; i++) gzprintf(gzout, ":%s", seq_paths[i]);
  gzputc(gzout, '\n');

  gzprintf(gzout, "##ctxVersion=\""VERSION_STATUS_STR"\"\n");
  gzprintf(gzout, "##ctxKmerSize=%i\n", MAX_KMER_SIZE);
}

void breakpoints_call(size_t num_of_threads,
                      const read_t *reads, size_t num_reads,
                      gzFile gzout, const char *out_path,
                      const char **seq_paths, size_t num_seq_paths,
                      const CmdArgs *args, dBGraph *db_graph)
{
  KOGraph kograph = kograph_create(reads, num_reads, true,
                                   num_of_threads, db_graph);
  BreakpointCaller *callers = brkpt_callers_new(num_of_threads, gzout,
                                                kograph, db_graph);

  status("Running BreakpointCaller with %zu thread%s, output to: %s",
         num_of_threads, num_of_threads == 1 ? "" : "s", out_path);

  breakpoints_print_header(gzout, args, seq_paths, num_seq_paths);

  util_run_threads(callers, num_of_threads, sizeof(callers[0]),
                   num_of_threads, breakpoint_caller);

  brkpt_callers_destroy(callers, num_of_threads);
  kograph_free(kograph);
}
