#include "global.h"
#include "db_graph.h"
#include "read_paths.h"
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

struct AddShades {
  dBNodeBuffer list;
  path_t path;
};

static void dump_gap_sizes(const char *base_fmt, uint64_t *arr, size_t arrlen,
                           uint32_t kmer_size)
{
  StrBuf *csv_dump = file_reader_generate_filename(base_fmt);
  FILE *fh;

  if(csv_dump == NULL) {
    warn("Cannot dump gapsize");
    return;
  }

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

  printf("shade contig gap sizes dumped to %s\n", csv_dump->buff);

  fclose(fh);
  strbuf_free(csv_dump);
}


// Assume colour bitfield is set up
static void add_read_path(const dBNodeBuffer *list,
                          dBGraph *graph, Colour colour,
                          path_t *path)
{
  if(list->len < 3) return;

  const hkey_t *nodes = list->nodes;
  const Orientation *orients = list->orients;

  uint32_t kmer_size = graph->kmer_size;

  // Path data
  binary_paths_t *paths = &graph->pdata;
  uint64_t *kmer_paths = graph->kmer_paths;

  // Edges
  Edges *edges = graph->edges;

  Nucleotide nuc_fw[list->len], nuc_rv[list->len];
  uint32_t pos_fw[list->len], pos_rv[list->len];

  // Find forks
  size_t pos, i, j, num_fw = 0, num_rv = 0;
  Nucleotide nuc;

  for(pos = 1; pos+1 < list->len; pos++)
  {
    hkey_t hkey = nodes[pos];
    Edges e = edges[hkey];
    Orientation or = orients[pos];

    if(edges_get_outdegree(e, or) > 0) {
      nuc = db_node_last_nuc(db_graph_bkmer(graph,nodes[pos+1]), or, kmer_size);
      nuc_fw[num_fw] = nuc;
      pos_fw[num_fw++] = pos;
    }
    if(edges_get_indegree(e, or) > 0) {
      nuc = db_node_first_nuc(db_graph_bkmer(graph,nodes[pos-1]), !or, kmer_size);
      nuc_rv[num_rv] = binary_nucleotide_complement(nuc);
      pos_rv[num_rv++] = pos;
    }
  }

  // Unsure path is long enough
  // Add is faster than MAX2
  size_t len = round_bits_to_bytes((num_fw+num_rv) * 2);
  if(len > path->cmpcap) {
    path->cmpcap = ROUNDUP2POW(path->cmpcap);
    path->cmpctseq = realloc(path->cmpctseq, path->cmpcap * sizeof(uint8_t));
  }

  // Reverse rv
  size_t tmp;
  for(i = 0, j = num_rv-1; i < j; i++, j--) {
    SWAP(pos_rv[i], pos_rv[j], tmp);
    SWAP(nuc_rv[i], nuc_rv[j], nuc);
  }

  // to add a path
  hkey_t hkey;
  uint64_t newpath;

  Nucleotide *bases = path->bases;

  // Generate paths
  for(i = 0; i < num_rv; i++) {
    for(j = 0; j < num_fw && pos_fw[j] <= pos_rv[i]; j++)
    {
      // start = pos_rv[i]-1;
      // end = pos_fw[j]+1;

      // Add forward paths
      hkey = nodes[pos_rv[i]-1];

      for(len = 1; j+len <= num_fw; len++)
      {
        path->core.prev = kmer_paths[hkey];
        path->core.len = len;
        path->bases = nuc_fw + j;
        newpath = binary_paths_add(paths, path, colour);
        if(newpath != PATH_NULL) kmer_paths[hkey] = newpath;
      }

      // Add reverse paths
      hkey = nodes[pos_fw[j]+1];

      for(len = 1; len <= i+1; len++)
      {
        path->core.prev = kmer_paths[hkey];
        path->core.len = len;
        path->bases = nuc_rv + i;
        newpath = binary_paths_add(paths, path, colour);
        if(newpath != PATH_NULL) kmer_paths[hkey] = newpath;
      }
    }
  }

  path->bases = bases;
}

// fill in gap in read node1==>---<==node2
//
// If successful: traverses gap and adds new nodes including node2/orient2
//    returns total number of nodes added (>= 0)
// If unsucessful: doesn't add anything to the list, returns -1
static int traverse_gap(dBNodeBuffer *list,
                        hkey_t node2, Orientation orient2,
                        dBGraph *db_graph, int colour)
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
    binary_kmer_to_str(node1->kmer, db_graph->kmer_size, tmp1);
    binary_kmer_to_str(node2->kmer, db_graph->kmer_size, tmp2);
    printf("traverse gap: %s:%i -> %s:%i\n", tmp1, orient1, tmp2, orient2);
  #endif

  GraphWalker wlk;

  // Walk from left -> right
  graph_walker_init(&wlk, db_graph, colour, node1, orient1);
  db_node_set_traversed(db_graph, wlk.node, wlk.orient);

  int i, pos = 0;

  while(pos < GAP_LIMIT && graph_traverse(&wlk) &&
        !db_node_has_traversed(db_graph, wlk.node, wlk.orient))
  {
    db_node_set_traversed(db_graph, wlk.node, wlk.orient);

    nlist[pos] = wlk.node;
    olist[pos] = wlk.orient;
    pos++;

    if(wlk.node == node2 && wlk.orient == orient2)
      break;
  }

  db_node_fast_clear_traversed(db_graph, node1);
  for(i = 0; i < pos; i++)
    db_node_fast_clear_traversed(db_graph, nlist[i]);

  if(wlk.node == node2 && wlk.orient == orient2) {
    list->len += pos;
    return pos;
  }

  int num_left = pos;

  // Look for the last node
  hkey_t tgt_n = list->nodes[list->len+num_left-1];
  Orientation tgt_o = opposite_orientation(list->orients[list->len+num_left-1]);

  // Walk from right to num_left
  graph_walker_init(&wlk, db_graph, colour, node2, opposite_orientation(orient2));
  db_node_set_traversed(db_graph, wlk.node, wlk.orient);

  pos = GAP_LIMIT-1;
  nlist[pos] = node2;
  olist[pos] = orient2;
  pos--;
  // pos is the next index at which to add a node

  boolean success = false;

  while(pos >= num_left && graph_traverse(&wlk) &&
        !db_node_has_traversed(db_graph, wlk.node, wlk.orient))
  {
    if(wlk.node == tgt_n && wlk.orient == tgt_o)
    {
      success = true;
      break;
    }

    db_node_set_traversed(db_graph, wlk.node, wlk.orient);

    nlist[pos] = wlk.node;
    olist[pos] = opposite_orientation(wlk.orient);
    pos--;
  }

  for(i = pos+1; i < GAP_LIMIT; i++)
    db_node_fast_clear_traversed(db_graph, nlist[i]);

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
static size_t parse_contig(dBNodeBuffer *list, dBGraph *graph, Colour colour,
                           const char *contig, size_t contig_len,
                           int mp_first_kmer, path_t *path)
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
    Nucleotide nuc = char_to_binary_nucleotide(contig[next_base++]);
    binary_kmer_left_shift_one_base(bkmer, kmer_size);
    binary_kmer_set_last_nuc(bkmer, nuc);
    db_node_get_key(bkmer, kmer_size, tmp_key);
    node = hash_table_find(&graph->ht, tmp_key);

    if(node != HASH_NOT_FOUND)
    {
      // Did we find any nodes from this contig in the graph?
      orientation = db_node_get_orientation(bkmer, tmp_key);

      if(list->len > 0 && prev_node == HASH_NOT_FOUND)
      {
        // Can we branch the gap from the prev contig?
        int gapsize = traverse_gap(list, node, orientation, graph, colour);

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
void load_shades(read_t *r1, read_t *r2,
                 int fq_offset1, int fq_offset2,
                 SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                 void *ptr)
{
  // Don't bother checking for duplicates
  (void)stats;
  struct AddShades *add_shades = (struct AddShades*)ptr;

  dBNodeBuffer *list = &add_shades->list;
  path_t *path = &add_shades->path;

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
    parse_contig(list, graph, prefs->into_colour,
                 r1->seq.b + contig_start, contig_len, -1, path);
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
      parse_contig(list, graph, prefs->into_colour,
                   r2->seq.b + contig_start, contig_len, mp_first_kmer, path);
    }
  }

  // load the current supercontig
  add_read_path(list, graph, prefs->into_colour, path);
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

  struct AddShades tmp_pair;
  db_node_buf_alloc(&tmp_pair.list, 4096);
  path_alloc(&tmp_pair.path);

  // load se data
  if(se_list != NULL)
  {
    parse_filelists(se_list, NULL, READ_FALIST, &prefs, stats,
                    &load_shades, &tmp_pair);
  }

  // load pe data
  if(pe_list1 != NULL && pe_list1[0] != '\0')
  {
    parse_filelists(pe_list1, pe_list2, READ_FALIST,
                    &prefs, stats, &load_shades, &tmp_pair);
  }

  // Load colour list
  prefs.into_colour = col_list_first_colour;

  if(colour_list != NULL)
  {
    parse_filelists(colour_list, NULL, READ_COLOURLIST,
                    &prefs, stats, &load_shades, &tmp_pair);
  }

  path_dealloc(&tmp_pair.path);
  db_node_buf_dealloc(&tmp_pair.list);

  seq_loading_stats_free(stats);

  // Print mp gap size / insert stats to a file
  uint32_t kmer_size = prefs.db_graph->kmer_size;
  dump_gap_sizes("gap_sizes.%u.csv", gap_sizes, GAP_LIMIT, kmer_size);
  dump_gap_sizes("mp_sizes.%u.csv", insert_sizes, GAP_LIMIT, kmer_size);

  message("Paths added");
}
