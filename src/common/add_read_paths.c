#include "global.h"
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

// Store stats here
static uint64_t insert_sizes[GAP_LIMIT] = {0};
static uint64_t gap_sizes[GAP_LIMIT] = {0};

// Temp data store for adding paths to graph
struct AddPaths {
  dBNodeBuffer list;
  path_t path;
  GraphWalker wlk;
  uint64_t *visited;
};

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

  printf("Contig gap sizes dumped to %s\n", csv_dump->buff);

  fclose(fh);
  strbuf_free(csv_dump);
}


static void add_read_path(const dBNodeBuffer *list,
                          dBGraph *graph, Colour colour,
                          path_t *path)
{
  if(list->len < 3) return;

  const hkey_t *nodes = list->nodes;
  const Orientation *orients = list->orients;

  uint32_t kmer_size = graph->kmer_size;

  Edges edges[list->len];
  int indegree[list->len], outdegree[list->len];
  Nucleotide nuc_fw[list->len], nuc_rv[list->len];
  uint32_t pos_fw[list->len], pos_rv[list->len];

  // Find forks
  size_t i, j, num_fw = 0, num_rv = 0;
  Nucleotide nuc;

  for(i = 0; i < list->len; i++) {
    edges[i] = graph->edges[nodes[i]];
    outdegree[i] = edges_get_outdegree(edges[i], orients[i]);
    indegree[i] = edges_get_indegree(edges[i], orients[i]);
  }

  #ifdef DEBUG
    char str[100];
    BinaryKmer bkmer;
    db_graph_oriented_bkmer(graph, nodes[0], orients[0], bkmer);
    binary_kmer_to_str(bkmer, kmer_size, str);
    printf("%s", str);
    for(i = 1; i < list->len; i++) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph, nodes[i]);
      char c = db_node_last_nuc(bkmerptr, orients[i], kmer_size);
      putc(binary_nuc_to_char(c), stdout);
    }
    printf("\n");
  #endif

  for(i = 0; i < list->len; i++)
  {
    if(i+1 < list->len && outdegree[i] > 1) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph,nodes[i+1]);
      nuc = db_node_last_nuc(bkmerptr, orients[i+1], kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = i;
      // Ensure this nodes is in this colour
      db_node_set_col(graph, nodes[i+1], colour);
    }
    if(i > 0 && indegree[i] > 1) {
      ConstBinaryKmerPtr bkmerptr = db_node_bkmer(graph, nodes[i-1]);
      nuc = orients[i-1] == forward
              ? binary_nuc_complement(binary_kmer_first_nuc(bkmerptr, kmer_size))
              : binary_kmer_last_nuc(bkmerptr);
      nuc_rv[num_rv] = nuc;
      pos_rv[num_rv++] = i;
      // Ensure this nodes is in this colour
      db_node_set_col(graph, nodes[i-1], colour);
    }
  }

  if(num_rv == 0 || num_fw == 0) return;

  // Reverse rv
  size_t tmp;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp);
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

  // Path data
  binary_paths_t *paths = &graph->pdata;

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

  for(start_rv = 0, start_fw = num_fw-1; ; start_fw--)
  {
    while(start_rv < num_rv && pos_rv[start_rv] > pos_fw[start_fw]) start_rv++;
    if(start_rv == num_rv) break;

    size_t pos = pos_fw[start_fw] + 1;
    if(start_rv > 0 && pos_rv[start_rv-1] == pos) start_rv--;

    node = nodes[pos];
    orient = rev_orient(orients[pos]);
    path->core.prev = db_node_paths(graph, node, orient);
    path->bases = nuc_rv + start_rv;
    path->core.len = num_rv - start_rv;

    #ifdef DEBUG
      binary_kmer_to_str(db_node_bkmer(graph, node), graph->kmer_size, str);
      printf(" %s:%i) start_rv: %zu start_fw: %zu {%i}\n", str, orient,
             start_rv, start_fw, pos_fw[start_fw]);
    #endif

    if((pindex = binary_paths_add(paths, path, colour)) != PATH_NULL)
      db_node_paths(graph, node, orient) = pindex;

    if(start_fw == 0) break;
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

    node = nodes[pos];
    orient = orients[pos];
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
  }

  // Restore saved
  path->bases = tmp_bases;
}

// fill in gap in read node1==>---<==node2
//
// If successful: traverses gap and adds new nodes including node2/orient2
//    returns total number of nodes added (>= 0)
// If unsucessful: doesn't add anything to the list, returns -1
static int traverse_gap(dBNodeBuffer *list,
                        hkey_t node2, Orientation orient2,
                        const dBGraph *db_graph, uint64_t *visited, int colour,
                        GraphWalker *wlk)
{
  hkey_t node1 = list->nodes[list->len-1];
  Orientation orient1 = list->orients[list->len-1];

  // First and last node already match
  if(node1 == node2 && orient1 == orient2) return 0;

  // Ensure capacity
  db_node_buf_ensure_capacity(list, list->len + GAP_LIMIT + 1);

  hkey_t *nlist = list->nodes + list->len;
  Orientation *olist = list->orients + list->len;

  #ifdef DEBUG
    char tmp1[100], tmp2[100];
    binary_kmer_to_str(db_node_bkmer(db_graph, node1), db_graph->kmer_size, tmp1);
    binary_kmer_to_str(db_node_bkmer(db_graph, node2), db_graph->kmer_size, tmp2);
    printf("traverse gap: %s:%i -> %s:%i\n", tmp1, orient1, tmp2, orient2);
  #endif

  // Walk from left -> right
  graph_walker_init(wlk, db_graph, colour, node1, orient1);
  db_node_set_traversed(visited, wlk->node, wlk->orient);

  int i, pos = 0;

  while(pos < GAP_LIMIT && graph_traverse(wlk) &&
        !db_node_has_traversed(visited, wlk->node, wlk->orient))
  {
    db_node_set_traversed(visited, wlk->node, wlk->orient);

    nlist[pos] = wlk->node;
    olist[pos] = wlk->orient;
    pos++;

    if(wlk->node == node2 && wlk->orient == orient2)
      break;
  }

  db_node_fast_clear_traversed(visited, node1);
  for(i = 0; i < pos; i++)
    db_node_fast_clear_traversed(visited, nlist[i]);

  if(wlk->node == node2 && wlk->orient == orient2) {
    list->len += pos;
    return pos;
  }

  int num_left = pos;

  graph_walker_finish(wlk);

  // Look for the last node
  hkey_t tgt_n = list->nodes[list->len+num_left-1];
  Orientation tgt_o = opposite_orientation(list->orients[list->len+num_left-1]);

  // Walk from right to num_left
  graph_walker_init(wlk, db_graph, colour, node2, opposite_orientation(orient2));
  db_node_set_traversed(visited, wlk->node, wlk->orient);

  pos = GAP_LIMIT-1;
  nlist[pos] = node2;
  olist[pos] = orient2;
  pos--;
  // pos is the next index at which to add a node

  boolean success = false;

  while(pos >= num_left && graph_traverse(wlk) &&
        !db_node_has_traversed(visited, wlk->node, wlk->orient))
  {
    if(wlk->node == tgt_n && wlk->orient == tgt_o)
    {
      success = true;
      break;
    }

    db_node_set_traversed(visited, wlk->node, wlk->orient);

    nlist[pos] = wlk->node;
    olist[pos] = opposite_orientation(wlk->orient);
    pos--;
  }

  graph_walker_finish(wlk);

  for(i = pos+1; i < GAP_LIMIT; i++)
    db_node_fast_clear_traversed(visited, nlist[i]);

  if(success)
  {
    int num_right = GAP_LIMIT - 1 - pos;
    memmove(nlist+num_left, nlist + pos + 1, num_right * sizeof(hkey_t));
    memmove(olist+num_left, olist + pos + 1, num_right * sizeof(Orientation));
    int num_kmers_added = num_left + num_right;
    list->len += num_kmers_added;
    return num_kmers_added;
  }

  return -1;
}

// Returns offset to continue from
// stop if contig_len-offset < kmer_size
// mp_first_kmer is the number of bp into the second mp of the first contig
// or -1 if not mate pair read or not first contig from mp read
static size_t parse_contig(dBNodeBuffer *list, dBGraph *graph, uint64_t *visited,
                           Colour colour, const char *contig, size_t contig_len,
                           int mp_first_kmer, path_t *path,
                           GraphWalker *wlk)
{
  uint32_t kmer_size = graph->kmer_size;

  #ifdef DEBUG
    printf("loading contig: %.*s [%i]\n", (int)contig_len, contig, (int)contig_len);
  #endif

  // Find the first kmer from the contig that is in the graph
  BinaryKmer bkmer, tmp_key;
  hkey_t node, prev_node = HASH_NOT_FOUND;
  Orientation orientation;

  binary_kmer_from_str(contig, kmer_size, bkmer);
  binary_kmer_right_shift_one_base(bkmer);

  size_t next_base = kmer_size-1;
  while(next_base < contig_len)
  {
    Nucleotide nuc = binary_nuc_from_char(contig[next_base++]);
    binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    db_node_get_key(bkmer, kmer_size, tmp_key);
    node = hash_table_find(&graph->ht, tmp_key);

    // char tmp[100];
    // binary_kmer_to_str(bkmer, kmer_size, tmp);
    // printf("%s [%zu]\n", tmp, (size_t)node);

    if(node != HASH_NOT_FOUND)
    {
      // Did we find any nodes from this contig in the graph?
      orientation = db_node_get_orientation(bkmer, tmp_key);

      if(list->len > 0 && prev_node == HASH_NOT_FOUND)
      {
        // Can we branch the gap from the prev contig?
        int gapsize = traverse_gap(list, node, orientation, graph, visited,
                                   colour, wlk);

        if(gapsize == -1)
        {
          // Failed to bridge gap
          add_read_path(list, graph, colour, path);

          list->nodes[0] = node;
          list->orients[0] = orientation;
          list->len = 1;
        }
        else
        {
          // Update stats (gapsize is already constrained to be <= GAP_LIMIT 
          if(mp_first_kmer != -1)
            insert_sizes[mp_first_kmer <= gapsize ? gapsize-mp_first_kmer : 0]++;
          else
            gap_sizes[gapsize]++;
        }
      }
      else
      {
        list->nodes[list->len] = node;
        list->orients[list->len] = orientation;
        list->len++;
      }
    }
    prev_node = node;
  }

  return (next_base+1) - kmer_size;
}


// This function is passed to parse_filelist to load paths from sequence data
void load_paths(read_t *r1, read_t *r2,
                int fq_offset1, int fq_offset2,
                SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                void *ptr)
{
  // Don't bother checking for duplicates
  (void)stats;
  struct AddPaths *add_paths = (struct AddPaths*)ptr;

  dBNodeBuffer *list = &add_paths->list;
  path_t *path = &add_paths->path;
  GraphWalker *wlk = &add_paths->wlk;
  uint64_t *visited = add_paths->visited;

  dBGraph *graph = prefs->db_graph;

  if(r1->seq.end < graph->kmer_size) return;

  int qcutoff1 = prefs->quality_cutoff;
  int qcutoff2 = prefs->quality_cutoff;

  if(prefs->quality_cutoff > 0)
  {
    qcutoff1 += fq_offset1;
    qcutoff2 += fq_offset2;
  }

  int hp_cutoff = prefs->homopolymer_cutoff;

  size_t contig_start, contig_end, search_start = 0;

  while((contig_start = seq_contig_start(r1, search_start, graph->kmer_size,
                                         qcutoff1, hp_cutoff)) < r1->seq.end)
  {
    contig_end = seq_contig_end(r1, contig_start, graph->kmer_size,
                                qcutoff1, hp_cutoff, &search_start);

    size_t contig_len = contig_end - contig_start;
    parse_contig(list, graph, visited, prefs->into_colour,
                 r1->seq.b + contig_start, contig_len, -1, path, wlk);
  }

  if(r2 != NULL && r2->seq.end >= graph->kmer_size)
  {
    seq_read_reverse_complement(r2);

    search_start = 0;

    while((contig_start = seq_contig_start(r2, search_start, graph->kmer_size,
                                           qcutoff2, hp_cutoff)) < r2->seq.end)
    {
      // mp_first_kmer is -1 for all but the first iteration of this loop
      int mp_first_kmer = (search_start == 0 ? (int)contig_start : -1);
      contig_end = seq_contig_end(r2, contig_start, graph->kmer_size,
                                  qcutoff2, hp_cutoff, &search_start);

      size_t contig_len = contig_end - contig_start;
      parse_contig(list, graph, visited, prefs->into_colour,
                   r2->seq.b + contig_start, contig_len, mp_first_kmer, path, wlk);
    }
  }

  // load the current supercontig
  add_read_path(list, graph, prefs->into_colour, path);
  list->len = 0;
}

void add_read_paths_to_graph(const char *se_list,
                             const char *pe_list1, const char *pe_list2,
                             Colour seq_colour,
                             const char *colour_list,
                             Colour col_list_first_colour,
                             SeqLoadingPrefs prefs)
{
  // Reset values we don't want in SeqLoadingPrefs
  prefs.remove_dups_se = false;
  prefs.remove_dups_pe = false;
  prefs.load_binaries = false;
  prefs.update_ginfo = false;
  prefs.into_colour = seq_colour;

  SeqLoadingStats *stats = seq_loading_stats_create(0);

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

  seq_loading_stats_free(stats);

  // Print mp gap size / insert stats to a file
  uint32_t kmer_size = prefs.db_graph->kmer_size;
  dump_gap_sizes("gap_sizes.%u.csv", gap_sizes, GAP_LIMIT, kmer_size);
  dump_gap_sizes("mp_sizes.%u.csv", insert_sizes, GAP_LIMIT, kmer_size);

  binary_paths_t *pdata = &prefs.db_graph->pdata;
  size_t bytes_used = pdata->next - pdata->store;
  char mem_used[100];
  bytes_to_str(bytes_used, 1, mem_used);

  message("Paths added\n");
  message("Currently %s used for %zu paths\n", mem_used, pdata->num_paths);

  free(tmpdata.visited);
}
