#include "global.h"
#include <pthread.h>

#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "util.h"
#include "file_util.h"
#include "seq_reader.h"
#include "file_reader.h"

// Maximum distance we wish to assemble across
#define MAX_SHADE_CONTIG 1000

// Biggest gap between reads we'll try to traverse
#define GAP_LIMIT 1000

// Temp data store for adding paths to graph
struct AddPaths {
  dBNodeBuffer list;
  path_t path;
  GraphWalker wlk;
  uint64_t *visited;
  // uint64_t insert_sizes[GAP_LIMIT+1] = {0};
  // uint64_t gap_sizes[GAP_LIMIT+1] = {0};
};

// DEV: move these to AddPaths
uint64_t insert_sizes[GAP_LIMIT+1] = {0};
uint64_t gap_sizes[GAP_LIMIT+1] = {0};

pthread_mutex_t add_paths_mutex;

static void dump_gap_sizes(const char *base_fmt, uint64_t *arr, size_t arrlen,
                           uint32_t kmer_size)
{
  StrBuf *csv_dump = strbuf_new();

  if(!file_reader_generate_filename(base_fmt, csv_dump)) {
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

  printf("  Contig gap sizes dumped to %s\n", csv_dump->buff);

  fclose(fh);
  strbuf_free(csv_dump);
}

static void add_read_path(const dBNode *nodes, size_t len,
                          dBGraph *graph, Colour colour,
                          path_t *path)
{
  if(len < 3) return;

  binary_paths_t *paths = &graph->pdata;
  uint32_t kmer_size = graph->kmer_size;

  Edges edges[len];
  int indegree[len], outdegree[len];
  Nucleotide nuc_fw[len], nuc_rv[len];
  uint32_t pos_fw[len], pos_rv[len];

  // Find forks
  size_t i, j, num_fw = 0, num_rv = 0;
  Nucleotide nuc;

  for(i = 0; i < len; i++) {
    edges[i] = graph->edges[nodes[i].node];
    outdegree[i] = edges_get_outdegree(edges[i], nodes[i].orient);
    indegree[i] = edges_get_indegree(edges[i], nodes[i].orient);
  }

  #ifdef DEBUG
    char str[100];
    BinaryKmer bkmer;
    db_graph_oriented_bkmer(graph, nodes[0].node, nodes[0].orient, bkmer);
    binary_kmer_to_str(bkmer, kmer_size, str);
    printf("%s", str);
    for(i = 1; i < len; i++) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph, nodes[i].node);
      char c = db_node_last_nuc(bkmerptr, nodes[i].orient, kmer_size);
      putc(binary_nuc_to_char(c), stdout);
    }
    printf("\n");
  #endif

  for(i = 0; i < len; i++)
  {
    if(i+1 < len && outdegree[i] > 1) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph, nodes[i+1].node);
      nuc = db_node_last_nuc(bkmerptr, nodes[i+1].orient, kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
      // Ensure this nodes is in this colour
      db_node_set_col(graph, nodes[i+1].node, colour);
    }
    if(i > 0 && indegree[i] > 1) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph, nodes[i-1].node);
      nuc = nodes[i-1].orient == forward
              ? binary_nuc_complement(binary_kmer_first_nuc(bkmerptr, kmer_size))
              : binary_kmer_last_nuc(bkmerptr);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
      // Ensure this nodes is in this colour
      db_node_set_col(graph, nodes[i-1].node, colour);
    }
  }

  if(num_rv == 0 || num_fw == 0) return;

  // Reverse rv
  size_t tmp_pos;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp_pos);
    SWAP(nuc_rv[i], nuc_rv[j], nuc);
  }

  #ifdef DEBUG
    printf("fw ");
    for(i = 0; i < num_fw; i++)
      printf(" %i:%c", pos_fw[i], binary_nuc_to_char(nuc_fw[i]));
    printf("\n");

    printf("rv ");
    for(i = 0; i < num_rv; i++)
      printf(" %i:%c", pos_rv[i], binary_nuc_to_char(nuc_rv[i]));
    printf("\n");
  #endif

  path_init(path);
  bitset_set(path->core.colours, colour);

  // to add a path
  hkey_t node;
  Orientation orient;
  uint64_t pindex;

  // Store this temporarily
  Nucleotide *tmp_bases = path->bases;
  size_t start_fw, start_rv;

  // .//\/

  //
  // Generate paths going backwards through the contig
  //
  #ifdef DEBUG
    printf("==REV==\n");
  #endif

  // Get Lock
  pthread_mutex_lock(&add_paths_mutex);


  BinaryKmer tmpbkmer;
  binary_kmer_from_str("CTGATAATTTCACTCTCACAATGTAGGGGAA", 31, tmpbkmer);

  for(start_rv = 0, start_fw = num_fw-1; ; start_fw--)
  {
    while(start_rv < num_rv && pos_rv[start_rv] > pos_fw[start_fw]) start_rv++;
    if(start_rv == num_rv) break;

    size_t pos = pos_fw[start_fw] + 1;
    if(start_rv > 0 && pos_rv[start_rv-1] == pos) start_rv--;

    node = nodes[pos].node;
    orient = rev_orient(nodes[pos].orient);
    path->core.prev = db_node_paths(graph, node, orient);
    path->bases = nuc_rv + start_rv;
    path->core.len = num_rv - start_rv;

    #ifdef DEBUG
      binary_kmer_to_str(db_node_bkmer(graph, node), graph->kmer_size, str);
      printf(" %s:%i) start_rv: %zu start_fw: %zu {%i}\n", str, orient,
             start_rv, start_fw, pos_fw[start_fw]);
    #endif

    if((pindex = binary_paths_add(paths, path, colour)) != PATH_NULL) {
      if(binary_kmers_are_equal(db_node_bkmer(graph, node), tmpbkmer)) {
        char contig[1000];
        db_nodes_to_str(nodes, len, graph, contig);
        printf("hit: %s [%zu]\n", contig, len);
      }
      db_node_paths(graph, node, orient) = pindex;
    }

    if(start_fw == 0) break;
    // break;
  }

  //
  // Generate forward paths
  //
  #ifdef DEBUG
    printf("==FWD==\n");
  #endif

  for(start_fw = 0, start_rv = num_rv-1; ; start_rv--)
  {
    while(start_fw < num_fw && pos_fw[start_fw] < pos_rv[start_rv]) start_fw++;
    if(start_fw == num_fw) break;

    size_t pos = pos_rv[start_rv] - 1;
    if(start_fw > 0 && pos_fw[start_fw-1] == pos) start_fw--;

    node = nodes[pos].node;
    orient = nodes[pos].orient;
    path->core.prev = db_node_paths(graph, node, orient);
    path->bases = nuc_fw + start_fw;
    path->core.len = num_fw - start_fw;

    #ifdef DEBUG
      binary_kmer_to_str(db_node_bkmer(graph, node), graph->kmer_size, str);
      printf(" %s:%i) start_rv: %zu start_fw: %zu\n", str, orient,
             start_rv, start_fw);
    #endif

    if((pindex = binary_paths_add(paths, path, colour)) != PATH_NULL)
      db_node_paths(graph, node, orient) = pindex;

    if(start_rv == 0) break;
    // break;
  }

  // Free lock
  pthread_mutex_unlock(&add_paths_mutex);

  // Restore saved
  path->bases = tmp_bases;
}

// fill in gap in read node1==>---<==node2
// Adds at most GAP_LIMIT nodes to list
// If successful: traverses gap and adds new nodes including node2/orient2
//    returns total number of nodes added (>= 0)
// If unsucessful: doesn't add anything to the list, returns -1
static int traverse_gap(dBNodeBuffer *list,
                        hkey_t node2, Orientation orient2,
                        const dBGraph *db_graph, uint64_t *visited,
                        Colour colour, GraphWalker *wlk)
{
  hkey_t node1 = list->data[list->len-1].node;
  Orientation orient1 = list->data[list->len-1].orient;

  #ifdef DEBUG
    char tmp1[100], tmp2[100];
    binary_kmer_to_str(db_node_bkmer(db_graph, node1), db_graph->kmer_size, tmp1);
    binary_kmer_to_str(db_node_bkmer(db_graph, node2), db_graph->kmer_size, tmp2);
    printf("traverse gap: %s:%i -> %s:%i\n", tmp1, orient1, tmp2, orient2);
  #endif

  // First and last node already match
  if(node1 == node2 && orient1 == orient2) return 0;

  // Ensure capacity
  db_node_buf_ensure_capacity(list, list->len + GAP_LIMIT);

  dBNode *nodes = list->data + list->len;

  // Walk from left -> right
  graph_walker_init(wlk, db_graph, colour, node1, orient1);
  db_node_set_traversed(visited, wlk->node, wlk->orient);

  size_t i, pos = 0;

  while(pos < GAP_LIMIT && graph_traverse(wlk) &&
        !db_node_has_traversed(visited, wlk->node, wlk->orient))
  {
    db_node_set_traversed(visited, wlk->node, wlk->orient);

    nodes[pos].node = wlk->node;
    nodes[pos].orient = wlk->orient;
    pos++;

    if(wlk->node == node2 && wlk->orient == orient2)
      break;
  }

  graph_walker_finish(wlk);

  db_node_fast_clear_traversed(visited, node1);
  for(i = 0; i < pos; i++)
    db_node_fast_clear_traversed(visited, nodes[i].node);

  if(wlk->node == node2 && wlk->orient == orient2) {
    list->len += pos;
    return pos;
  }

  // Walk from right -> left
  graph_walker_init(wlk, db_graph, colour, node2, opposite_orientation(orient2));
  db_node_set_traversed(visited, wlk->node, wlk->orient);

  pos = GAP_LIMIT-1;
  nodes[pos].node = node2;
  nodes[pos].orient = orient2;
  // pos is now the index at which we last added a node

  boolean success = false;

  while(pos > 0 && graph_traverse(wlk) &&
        !db_node_has_traversed(visited, wlk->node, wlk->orient))
  {
    Orientation orient = opposite_orientation(wlk->orient);

    if(wlk->node == node1 && orient == orient1) {
      success = true;
      break;
    }

    db_node_set_traversed(visited, wlk->node, wlk->orient);

    pos--;
    nodes[pos].node = wlk->node;
    nodes[pos].orient = orient;
  }

  graph_walker_finish(wlk);

  for(i = pos; i < GAP_LIMIT; i++)
    db_node_fast_clear_traversed(visited, nodes[i].node);

  if(success)
  {
    size_t num = GAP_LIMIT - pos;
    memmove(nodes, nodes + GAP_LIMIT - num, num * sizeof(dBNode));
    list->len += num;
    return num;
  }

  return -1;
}


// This function is passed to parse_filelist to load paths from sequence data
static void load_paths(read_t *r1, read_t *r2,
                       int fq_offset1, int fq_offset2,
                       SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                       void *ptr)
{
  // Don't bother checking for duplicates
  (void)stats;

  // print_debug = (strncmp("CTGTCCACAGGTTATCTGCAGTACTGAAGAG", r1->seq.b, 31) == 0);
  // if(print_debug) {
  //   printf("READ1: %s\n", r1->seq.b);
  //   printf("READ2: %s\n", r2->seq.b);
  // }

  struct AddPaths *add_paths = (struct AddPaths*)ptr;
  dBNodeBuffer *list = &add_paths->list;
  path_t *path = &add_paths->path;
  GraphWalker *wlk = &add_paths->wlk;
  uint64_t *visited = add_paths->visited;
  dBGraph *db_graph = prefs->db_graph;
  Colour colour = prefs->into_colour;

  int qcutoff1 = prefs->quality_cutoff;
  int qcutoff2 = prefs->quality_cutoff;

  if(prefs->quality_cutoff > 0)
  {
    qcutoff1 += fq_offset1;
    qcutoff2 += fq_offset2;
  }

  int hp_cutoff = prefs->homopolymer_cutoff;
  size_t i, r2_start = 0;
  int r2_offset = -1;
  list->len = 0;
  boolean useful_path_info = false;
  hkey_t node;

  get_nodes_from_read(r1, qcutoff1, hp_cutoff, db_graph, list);

  if(r2 != NULL && r2->seq.end >= db_graph->kmer_size)
  {
    seq_read_reverse_complement(r2);

    // Insert gap
    if(list->len > 0)
    {
      list->data[list->len].node = HASH_NOT_FOUND;
      list->data[list->len].orient = forward;
      list->len++;
      r2_start = list->len;
    }

    r2_offset = get_nodes_from_read(r2, qcutoff2, hp_cutoff, db_graph, list);
  }

  // DEV: do hash look up in new thread
  // DEV: pass list to thread here

  for(i = 1; i < list->len; i++) {
    if((node = list->data[i].node) != HASH_NOT_FOUND &&
       edges_get_indegree(db_graph->edges[node], list->data[i].orient) > 0)
    {
      useful_path_info = true;
      break;
    }
  }

  if(r2 != NULL && !useful_path_info) {
    for(i = r2_start; i < list->len - 1; i++) {
      if((node = list->data[i].node) != HASH_NOT_FOUND &&
         edges_get_outdegree(db_graph->edges[node], list->data[i].orient) > 0)
      {
        useful_path_info = true;
        break;
      }
    }
  }

  if(list->len == 0 || !useful_path_info) return;

  for(i = 0; i < list->len && list->data[i].node != HASH_NOT_FOUND; i++) {}
  boolean no_gaps = (i == list->len);

  if(no_gaps)
  {
    add_read_path(list->data, list->len, db_graph, colour, path);
  }
  else
  {
    // Append nodes onto end of list whilst attempting to close gaps
    const size_t end = list->len;

    db_node_buf_ensure_capacity(list, list->len + list->len);

    hkey_t node = HASH_NOT_FOUND, prev_node = HASH_NOT_FOUND;
    Orientation orient;

    for(i = 0; i < end; i++, prev_node = node)
    {
      // printf("%zu / %zu\n", i, end);
      node = list->data[i].node;
      orient = list->data[i].orient;

      #ifdef DEBUG
        char str[100];
        printf("node: %zu / %zu\n", (size_t)node, (size_t)db_graph->ht.capacity);
        // if(node == HASH_NOT_FOUND) strcpy(str, "(none)");
        // else binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, str);
        printf("  node:%zu %zu %s\n", i, (size_t)node, str);
      #endif

      if(node != HASH_NOT_FOUND)
      {
        if(prev_node == HASH_NOT_FOUND && list->len > end)
        {
          // Can we branch the gap from the prev contig?
          int gapsize = traverse_gap(list, node, orient, db_graph, visited,
                                     colour, wlk);

          if(gapsize == -1)
          {
            // printf("  traverse failed\n");
            // Failed to bridge gap
            add_read_path(list->data+end, list->len-end, db_graph, colour, path);

            list->data[end].node = node;
            list->data[end].orient = orient;
            list->len = end+1;
          }
          else
          {
            // printf("  traverse good\n");
            // Update stats (gapsize is <= GAP_LIMIT)
            if(i == r2_start) {
              int gap = gapsize - r2_offset;
              insert_sizes[gap < 0 ? 0 : gap]++;
            }
            else {
              gap_sizes[gapsize]++;
            }
          }
        }
        else {
          list->data[list->len].node = node;
          list->data[list->len].orient = orient;
          list->len++;
        }
      }
    }

    add_read_path(list->data+end, list->len-end, db_graph, colour, path);
  }

  // print_debug = 0;
}

void add_read_paths_to_graph(const char *se_list,
                             const char *pe_list1, const char *pe_list2,
                             Colour seq_colour,
                             const char *colour_list,
                             Colour col_list_first_colour,
                             SeqLoadingPrefs prefs)//,
                             // int num_of_threads)
{
  // Reset values we don't want in SeqLoadingPrefs
  prefs.remove_dups_se = false;
  prefs.remove_dups_pe = false;
  prefs.load_binaries = false;
  prefs.update_ginfo = false;
  prefs.into_colour = seq_colour;

  if(pthread_mutex_init(&add_paths_mutex, NULL) != 0)
    die("mutex init failed");

  // DEV: create threads
  // PTHREAD_CREATE_JOINABLE
  // sig_atomic_t flag on each thread

  SeqLoadingStats *stats = seq_loading_stats_create(0);

  // struct AddPaths *tmpdata = malloc(sizeof(struct AddPaths));
  // db_node_buf_alloc(&tmpdata->list, 4096);
  // path_alloc(&tmpdata->path);
  // graph_walker_alloc(&tmpdata->wlk);

  struct AddPaths tmpdata;
  db_node_buf_alloc(&tmpdata.list, 4096);
  path_alloc(&tmpdata.path);
  graph_walker_alloc(&tmpdata.wlk);

  size_t visited_bytes = 2 * round_bits_to_words64(prefs.db_graph->ht.capacity);
  tmpdata.visited = calloc(visited_bytes, sizeof(uint64_t));

  if(tmpdata.visited == NULL) die("Out of memory");

  // load se data
  if(se_list != NULL)
  {
    parse_filelists(se_list, NULL, READ_FALIST, &prefs, stats,
                    &load_paths, &tmpdata);
  }

  // load pe data
  if(pe_list1 != NULL && pe_list1[0] != '\0')
  {
    parse_filelists(pe_list1, pe_list2, READ_FALIST,
                    &prefs, stats, &load_paths, &tmpdata);
  }

  // Load colour list
  prefs.into_colour = col_list_first_colour;

  if(colour_list != NULL)
  {
    parse_filelists(colour_list, NULL, READ_COLOURLIST,
                    &prefs, stats, &load_paths, &tmpdata);
  }

  graph_walker_dealloc(&tmpdata.wlk);
  path_dealloc(&tmpdata.path);
  db_node_buf_dealloc(&tmpdata.list);
  free(tmpdata.visited);
  // free(tmpdata);

  seq_loading_stats_free(stats);

  // DEV: merge threads (join)
  // DEV: merge gap_sizes / insert_sizes

  pthread_mutex_destroy(&add_paths_mutex);

  message("Paths added\n");

  // Print mp gap size / insert stats to a file
  uint32_t kmer_size = prefs.db_graph->kmer_size;
  dump_gap_sizes("gap_sizes.%u.csv", gap_sizes, GAP_LIMIT+1, kmer_size);
  dump_gap_sizes("mp_sizes.%u.csv", insert_sizes, GAP_LIMIT+1, kmer_size);

  char mem_used_str[100], num_paths_str[100];
  binary_paths_t *pdata = &prefs.db_graph->pdata;
  size_t paths_mem_used = pdata->next - pdata->store;
  bytes_to_str(paths_mem_used, 1, mem_used_str);
  ulong_to_str(pdata->num_paths, num_paths_str);
  message("Currently %s used for %s paths\n\n", mem_used_str, num_paths_str);
}
