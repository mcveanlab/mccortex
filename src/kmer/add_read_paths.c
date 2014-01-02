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

boolean print_traversed_inserts = true;

boolean add_read_paths_mutex_setup = false;
pthread_mutex_t add_paths_mutex;

void add_read_paths_init()
{
  if(!add_read_paths_mutex_setup) {
    add_read_paths_mutex_setup = true;
    if(pthread_mutex_init(&add_paths_mutex, NULL) != 0) die("mutex init failed");
  }
}

void add_read_paths_cleanup()
{
  if(add_read_paths_mutex_setup) {
    add_read_paths_mutex_setup = false;
    pthread_mutex_destroy(&add_paths_mutex);
  }
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
  boolean added, rv_paths_added = false;

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

    new_index = path_store_find_or_add(paths, prev_index, plen, bases,
                                       orient, ctp_col, &added);

    // If the path already exists, all of its subpaths also already exist
    if(!added) break;
    db_node_paths(db_graph, node) = new_index;
    rv_paths_added = true;
  }

  // If we added no reverse paths, there are no forward paths to add either
  // (they are all already in the PathStore)
  if(rv_paths_added)
  {
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

      new_index = path_store_find_or_add(paths, prev_index, plen, bases,
                                         orient, ctp_col, &added);

      // If the path already exists, all of its subpaths also already exist
      if(!added) break;
      db_node_paths(db_graph, node) = new_index;
    }
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
    printf(" [len: %zu]\n", len);
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
              ? dna_nuc_complement(binary_kmer_first_nuc(bkmer, graph->kmer_size))
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
// Adds at most ins_gap_max nodes to nodebuf
// If successful: traverses gap and adds new nodes including node2/orient2
//    returns total number of nodes added (>= 0)
// If unsucessful: doesn't add anything to the nodebuf, returns -1
static int traverse_gap(dBNodeBuffer *nodebuf,
                        hkey_t node2, Orientation orient2,
                        Colour ctxcol, Colour ctpcol,
                        size_t ins_gap_min, size_t ins_gap_max,
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
  db_node_buf_ensure_capacity(nodebuf, nodebuf->len + ins_gap_max+1);

  dBNode *nodes = nodebuf->data + nodebuf->len;
  size_t pos = 0;
  Nucleotide lost_nuc;

  // Walk from left -> right
  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node1, orient1);
  lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

  // need to call db_node_has_col only if more than one colour loaded
  while(pos <= ins_gap_max && graph_traverse(wlk) &&
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

  if(wlk->node == node2 && wlk->orient == orient2 && pos >= ins_gap_min) {
    nodebuf->len += pos;
    return pos;
  }

  // Walk from right -> left
  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node2, opposite_orientation(orient2));

  pos = ins_gap_max;
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

  walker_fast_clear(rptwlk, nodes+pos, ins_gap_max-pos);

  if(success && pos >= ins_gap_min)
  {
    #ifdef CTXVERBOSE
      printf(" traverse success\n");
    #endif
    size_t num = ins_gap_max - pos;
    memmove(nodes, nodes + ins_gap_max - num, num * sizeof(dBNode));
    nodebuf->len += num;
    return num;
  }

  return -1;
}

void add_read_paths(const AddPathsJob *job, dBNodeBuffer *nodebuf,
                    GraphWalker *wlk, RepeatWalker *rptwlk,
                    uint64_t *insert_sizes, uint64_t *gap_sizes,
                    dBGraph *db_graph)
{
  // dBNodeBuffer *nodebuf = &worker->nodebuf;
  // dBGraph *db_graph = worker->db_graph;
  // AddPathsJob *job = &worker->job;
  // const size_t ctx_col = job->ctx_col;

  size_t i, r2_start = 0;
  int r2_offset = -1;
  nodebuf->len = 0;

  seq_nodes_from_read(&job->r1, job->qcutoff1, job->hp_cutoff, db_graph, nodebuf);

  if(job->r2.seq.end >= db_graph->kmer_size)
  {
    if(nodebuf->len > 0)
    {
      // Insert gap
      db_node_buf_ensure_capacity(nodebuf, nodebuf->len+1);
      nodebuf->data[nodebuf->len].key = HASH_NOT_FOUND;
      nodebuf->data[nodebuf->len].orient = FORWARD;
      nodebuf->len++;
    }

    r2_start = nodebuf->len;
    r2_offset = seq_nodes_from_read(&job->r2, job->qcutoff2, job->hp_cutoff,
                                    db_graph, nodebuf);
  }

  if(nodebuf->len == 0) return;

  // Check for gaps
  for(i = 0; i < nodebuf->len && nodebuf->data[i].key != HASH_NOT_FOUND; i++) {}

  if(i == nodebuf->len)
  {
    // No gaps in contig
    add_read_path(nodebuf->data, nodebuf->len, db_graph,
                  job->ctx_col, job->ctp_col);
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
          int gapsize = traverse_gap(nodebuf, node, orient,
                                     job->ctx_col, job->ctp_col,
                                     job->ins_gap_min, job->ins_gap_max,
                                     wlk, rptwlk, db_graph);

          if(gapsize == -1)
          {
            // Failed to bridge gap
            add_read_path(nodebuf->data+end, nodebuf->len-end, db_graph,
                          job->ctx_col, job->ctp_col);

            nodebuf->data[end].key = node;
            nodebuf->data[end].orient = orient;
            nodebuf->len = end+1;
          }
          else
          {
            // Update stats (gapsize is <= ins_gap_max)
            if(i == r2_start)
            {
              // Kmer is start of second read
              int gap = gapsize - r2_offset;
              insert_sizes[gap < 0 ? 0 : gap]++;

              if(print_traversed_inserts) {
                printf(">GappedContig\n");
                db_nodes_print(nodebuf->data+end, nodebuf->len-end, db_graph, stdout);
                printf("\n");
              }
            }
            else {
              gap_sizes[gapsize]++;
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
                  job->ctx_col, job->ctp_col);
  }
}

void dump_gap_sizes(const char *base_fmt, const uint64_t *arr, size_t arrlen,
                    size_t kmer_size, boolean insert_sizes, size_t nreads)
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

  char ngaps_str[100], nreads_str[100];
  ulong_to_str(ngaps, ngaps_str);
  ulong_to_str(nreads, nreads_str);

  status("%s size distribution: "
         "min: %zu mean: %.1f median: %.1f mode: %zu max: %zu; n=%s / %s [%zu%%]",
         insert_sizes ? "Insert" : "Seq error gap",
         min, mean, median, mode, max, ngaps_str, nreads_str,
         (size_t)((100.0 * ngaps) / nreads + 0.5));

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
    size_t i, start = 0, end = arrlen-1;

    while(start < arrlen && arr[start] == 0) start++;
    while(end > start && arr[end] == 0) end--;

    for(i = start; i <= end; i++) {
      fprintf(fout, "%4zu\t%4li\t%4zu\n", i, (long)i-kmer_size, (size_t)arr[i]);
    }
  }

  status("Contig %s sizes dumped to %s\n",
         insert_sizes ? "insert" : "gap", csv_dump->buff);

  fclose(fout);
  strbuf_free(csv_dump);
}
