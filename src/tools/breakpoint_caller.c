#include "global.h"
#include "breakpoint_caller.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "multicol_walker.h"
#include "kmer_occur.h"

typedef struct
{
  // Specific to this instance
  const size_t threadid, nthreads;
  const size_t ncols; // Number of colours to call in

  // Temporary memory used by this instance
  size_t *cols, *cols_used, *col_lengths; // arrays of length ncols
  MulticolWalker walker;
  dBNodeBuffer nbuf;

  // Passed to all instances
  const KOGraph kograph;
  const dBGraph *db_graph;
  gzFile gzout;
  pthread_mutex_t *const out_lock;
} BreakpointCaller;

static BreakpointCaller* brkpt_callers_new(size_t num_callers, gzFile gzout,
                                           const KOGraph kograph,
                                           const dBGraph *db_graph)
{
  size_t i, ncols = db_graph->num_of_cols;

  BreakpointCaller *callers = ctx_malloc(num_callers * sizeof(BreakpointCaller));

  // Each caller uses 3 lists of max length ncols
  size_t *colour_lists = ctx_malloc(num_callers * ncols * 3 * sizeof(size_t));

  pthread_mutex_t *out_lock = ctx_malloc(sizeof(pthread_mutex_t));
  if(pthread_mutex_init(out_lock, NULL) != 0) die("mutex init failed");

  for(i = 0; i < num_callers; i++)
  {
    BreakpointCaller tmp = {.threadid = i, .nthreads = num_callers,
                            .cols        = colour_lists,
                            .cols_used   = colour_lists+ncols,
                            .col_lengths = colour_lists+ncols+ncols,
                            .kograph = kograph,
                            .db_graph = db_graph,
                            .gzout = gzout,
                            .out_lock = out_lock};

    memcpy(&callers[i], &tmp, sizeof(BreakpointCaller));

    colour_lists += ncols*3;

    db_node_buf_alloc(&callers[i].nbuf, 1024);
    multicol_walker_alloc(&callers[i].walker, db_graph);
  }

  return callers;
}

static void brkpt_callers_destroy(BreakpointCaller *callers, size_t num_callers)
{
  size_t i;
  for(i = 0; i < num_callers; i++) {
    db_node_buf_dealloc(&callers[i].nbuf);
    multicol_walker_dealloc(&callers[i].walker);
  }
  ctx_free(callers[0].cols);
  pthread_mutex_destroy(callers[0].out_lock);
  ctx_free(callers[0].out_lock);
  ctx_free(callers);
}

static void process_contig(BreakpointCaller *caller, dBNodeBuffer *nbuf)
{
  gzFile gzout = caller->gzout;

  // Work backwards to find last place we met the ref
  // nbuf[0] is ref node
  // nbuf[1] is first node not in ref
  ctx_assert(kograph_num(caller->kograph, nbuf->data[0].key) > 0);
  ctx_assert(kograph_num(caller->kograph, nbuf->data[1].key) == 0);

  // AGGGCGTTAGCGGGTTGGAG
  printf("we got: ");
  db_nodes_print(nbuf->data, nbuf->len, caller->db_graph, stdout);
  printf("\n");

  size_t i, end;
  for(i = nbuf->len-1; kograph_num(caller->kograph, nbuf->data[i].key) == 0; i--);

  if(i == 0) {
    status("Never met ref");
    return; // we never met the ref again
  }

  // Found second ref meeting - find first node in block that is in ref
  end = i;
  while(kograph_num(caller->kograph, nbuf->data[end-1].key) > 0) end--;

  ctx_assert(end > 1); // nbuf[1] should not be in the ref

  pthread_mutex_lock(caller->out_lock);

  gzprintf(gzout, ">call.5pflank\n");
  // DEV
  gzprintf(gzout, "\n");

  gzprintf(gzout, ">call.3pflank\n");
  db_nodes_gzprint_cont(nbuf->data+end, nbuf->len - end, caller->db_graph, gzout);

  gzprintf(gzout, ">call.path\n");
  db_nodes_gzprint_cont(nbuf->data, end, caller->db_graph, gzout);
  // DEV

  pthread_mutex_unlock(caller->out_lock);
}

// DEV: for speed could remember where colours split off
static void assemble_contig(BreakpointCaller *caller, dBNodeBuffer *nbuf)
{
  // Follow path using multicol_walker_assemble_contig
  size_t i, ncols = caller->db_graph->num_of_cols, ncols_used;

  // Reset list of colours we are considering
  for(i = 0; i < ncols; i++) caller->cols[i] = i;

  while(ncols > 0)
  {
    ncols_used = multicol_walker_assemble_contig(&caller->walker,
                                                 caller->cols, caller->ncols,
                                                 caller->cols_used,
                                                 caller->col_lengths, nbuf);

    printf("ncols: %zu ncols_used: %zu\n", ncols, ncols_used);

    // Print with path
    process_contig(caller, &caller->nbuf);

    // Remove the used colours from list of colours to use
    ncols = multicol_walker_rem_cols(caller->cols, caller->ncols,
                                     caller->cols_used, ncols_used);
  }
}

// Walk the graph remembering the last time we met the ref
// When traversal fails, dump sequence up to last meeting with the ref
static void follow_break(BreakpointCaller *caller, dBNode node)
{
  size_t i, j, num_next;
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

  status("Possumable fork!");

  num_next = j;

  // 2. Follow all paths not in ref, in all colours
  dBNodeBuffer *nbuf = &caller->nbuf;

  for(i = 0; i < num_next; i++)
  {
    db_node_buf_reset(nbuf);
    nbuf->data[0] = node;
    nbuf->data[1] = next_nodes[i];
    nbuf->len = 2;
    assemble_contig(caller, nbuf);
  }
}

void breakpoint_caller_node(hkey_t hkey, BreakpointCaller *caller)
{
  Edges edges;

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
  gzprintf(gzout, "##ctxKmerSize=%i\n", MAX_KMER_SIZE);
}

void breakpoints_call(size_t num_of_threads,
                      const read_t *reads, size_t num_reads,
                      gzFile gzout, const char *out_path,
                      const CmdArgs *args, dBGraph *db_graph)
{
  KOGraph kograph = kograph_create(reads, num_reads, true,
                                   num_of_threads, db_graph);
  BreakpointCaller *callers = brkpt_callers_new(num_of_threads, gzout,
                                                kograph, db_graph);

  status("Running BreakpointCaller... output to: %s", out_path);

  breakpoints_print_header(gzout, args);

  util_run_threads(callers, num_of_threads, sizeof(callers[0]),
                   num_of_threads, breakpoint_caller);

  brkpt_callers_destroy(callers, num_of_threads);
  kograph_free(kograph);
}
