#include "global.h"

#include <stddef.h> // defines ptrdiff_t
#include <time.h> // printing datetime

#include "khash.h"

// cortex_var headers
#include "shaded_caller.h"
#include "db_graph.h"
#include "util.h"
#include "file_util.h"
#include "file_reader.h"
#include "seq_reader.h"

// Hash functions
#include "lookup3.h"
#include "city.h"

#include "binary_paths.h"
#include "graph_walker.h"

#ifdef DEBUG
#define DEBUG_CALLER 1
#endif

#define MAX_FLANK_KMERS 1000
#define MAX_ALLELE_KMERS 1000
#define MAX_WALK_BACK_KMERS 10

#define NODE_BUFFER_SIZE ((MAX_ALLELE_KMERS) * (NUM_OF_COLOURS) * 4)

// Print absolute path to a file
static void print_filepath_abs(gzFile *out, const char *name, const char *file)
{
  char absolute_path[PATH_MAX + 1];
  char *abs_path = realpath(file, absolute_path);

  if(abs_path == NULL)
    warn("Cannot get absolute path: %s\n", file);

  gzprintf(out, "##%s=%s\n", name, abs_path);
}

static void print_calling_header(dBGraph *db_graph, gzFile *out,
                                 const char* out_file)
                                 // const CmdLine *cmd)
{
  const GraphInfo *ginfo = &db_graph->ginfo;

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));
  char cwd[PATH_MAX + 1];

  gzprintf(out, "##fileformat=CTXv1.0\n");

  // Cortex details
  // gzprintf(out, "##ctxCmd=%s\n", cmd->cmdstr);

  if(file_reader_get_current_dir(cwd) != NULL)
    gzprintf(out, "##ctxCwd=%s\n", cwd);
  
  gzprintf(out, "##ctxDate=%s\n", datestr);
  gzprintf(out, "##ctxVersion=<version=%d.%d.%d.%d,kmer_size=%i,"
              "compiled_shades=%i,compiled_colours=%i>\n",
          VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION,
          db_graph->kmer_size, NUM_OF_SHADES, NUM_OF_COLOURS);

  // Print input sources
  // if(cmd->input_ctx_binary)
  //   print_filepath_abs(out, "ctx_load_binary", cmd->ctx_binary_path);
  // if(cmd->input_se_pe_lists && cmd->se_list != NULL)
  //   print_filepath_abs(out, "ctx_se_list", cmd->se_list);
  // if(cmd->input_se_pe_lists && cmd->pe_list_lh_mates != NULL) {
  //   print_filepath_abs(out, "ctx_pe_list1", cmd->pe_list_lh_mates);
  //   print_filepath_abs(out, "ctx_pe_list2", cmd->pe_list_rh_mates);
  // }
  // if(cmd->input_colours)
  //   print_filepath_abs(out, "ctx_colourlist", cmd->colour_list);

  print_filepath_abs(out, "ctx_bubbles_file", out_file);

  // Print ref colour if there is one
  // if(cmd->ref_colour != -1)
  //   gzprintf(out, "##ctx_refcol=%i\n", cmd->ref_colour);

  // Print colours we're calling in
  gzprintf(out, "##ctx_calling_in_colours=%i\n", ginfo->num_of_colours_loaded);

  StrBuf *sample_name = strbuf_new();

  // Print sample names
  Colour col;
  for(col = 0; col < ginfo->num_of_colours_loaded; col++)
  {
    const ErrorCleaning *ec = ginfo->cleaning + col;
    printf("Colour %u\n", col);

    // Find and replace double quotes with single quotes
    char *sname = ginfo->sample_names[col]->buff;

    if(strcmp(sname, "undefined") == 0 || strchr(sname, '\t') != NULL ||
       strchr(sname, ' ') != NULL || strchr(sname, '\r') != NULL ||
       strchr(sname, '\n') != NULL)
    {
      strbuf_reset(sample_name);
      strbuf_sprintf(sample_name, "sample%u", col);
    }
    else {
      strbuf_set(sample_name, ginfo->sample_names[col]->buff);
    }

    // gzprintf(out, "%s", sample_name->buff);
    // gzprintf(out, "%s", ginfo->sample_names[col]->buff);
    // // gzprintf(out, "%s", col);
    // gzprintf(out, "%zu", (size_t)ginfo->mean_read_length[col]);
    // gzprintf(out, "%zu", (size_t)ginfo->total_sequence[col]);
    // gzprintf(out, "%Lf", ginfo->seq_err[col]);
    // gzprintf(out, "%s", ec->tip_clipping ? "yes" : "no");
    // gzprintf(out, "%u", ec->remv_low_cov_sups_thresh);
    // gzprintf(out, "%u", ec->remv_low_cov_nodes_thresh);
    // gzprintf(out, "%s", ec->cleaned_against_graph_name->buff);

    gzprintf(out, "##SAMPLE=<ID=%s,name=\"%s\",colour=%i,"
                  "meanreadlen=%zu,totalseqloaded=%zu,"
                  "seqerror=%Lf,tipclipped=%s,removelowcovgsupernodes=%u,"
                  "removelowcovgkmer=%u,cleanedagainstgraph=%s>\n",
             sample_name->buff, ginfo->sample_names[col]->buff, col,
             (size_t)ginfo->mean_read_length[col],
             (size_t)ginfo->total_sequence[col], ginfo->seq_err[col],
             ec->tip_clipping ? "yes" : "no", ec->remv_low_cov_sups_thresh,
             ec->remv_low_cov_nodes_thresh,
             ec->cleaned_against_graph_name->buff);
  }

  strbuf_free(sample_name);
}

static void print_branch(hkey_t *nodes, Orientation *orients, size_t len,
                         char print_first_kmer, dBGraph *db_graph, gzFile out)
{
  size_t i;
  Nucleotide nuc;

  uint32_t kmer_size = db_graph->kmer_size;

  // Print sequence
  if(print_first_kmer)
  {
    BinaryKmer bkmer;
    char tmp[MAX_KMER_SIZE+1];
    db_graph_oriented_bkmer(db_graph, nodes[0], orients[0], bkmer);
    binary_kmer_to_str(bkmer, kmer_size, tmp);
    gzputs(out, tmp);
    i = 1;
  }
  else i = 0;

  for(; i < len; i++)
  {
    nuc = db_node_last_nuc(db_graph_bkmer(db_graph, nodes[i]),
                           orients[i], kmer_size);
    gzputc(out, binary_nuc_to_char(nuc));
  }

  gzputc(out, '\n');

  size_t col_loaded = db_graph->ginfo.num_of_colours_loaded;

  // Print covgs
  // if(len == 0) {
    for(i = 0; i < col_loaded; i++) gzputc(out, '\n');
  // }
  // else {
  //   for(i = 0; i < col_len; i++)
  //   {
  //     int colour = cols[i];
  //     gzprintf(out, "%i", nodes[0]->coverage[colour]);
  //     for(i = 1; i < len; i++)
  //       gzprintf(out, " %i", nodes[i]->coverage[colour]);
  //     gzputc(out, '\n');
  //   }
  // }
}

//
// Data types and structs for calling variants with the shaded caller
//
typedef Orientation SuperOrientation;

typedef struct CallerSupernode CallerSupernode;
typedef struct SupernodePath SupernodePath;
typedef struct SupernodePathPos SupernodePathPos;

struct CallerSupernode
{
  hkey_t *nodes;
  Orientation *orients;

  // ShadeSet *shades;
  size_t num_of_nodes;

  int num_prev, num_next;
  hkey_t prev_nodes[4], next_nodes[4];
  Orientation prev_orients[4], next_orients[4];

  SupernodePathPos *first_pathpos;
};

struct SupernodePath
{
  int colour;
  CallerSupernode *supernodes[MAX_ALLELE_KMERS];
  SuperOrientation superorients[MAX_ALLELE_KMERS];
  size_t length;
};

struct SupernodePathPos
{
  SupernodePath *path;
  size_t pos;
  SupernodePathPos *next;
};

// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_MAP_INIT_INT64(supnode_hsh, CallerSupernode*)

static void print_bubble(gzFile out, size_t bnum, dBGraph *db_graph,
                         SupernodePathPos **spp_arr, size_t num_of_paths,
                         hkey_t *flank5pe,
                         Orientation *flank5po,
                         size_t flank5pkmers);

#define supernode_get_orientation(snode,node,or) \
  ((node) == (snode)->nodes[0] && (or) == (snode)->orients[0] ? forward : reverse)

static void reverse_node_list(hkey_t *nlist, Orientation *olist, size_t len)
{
  if(len == 0) return;
  else if(len == 1)
  {
    olist[0] = opposite_orientation(olist[0]);
  }
  else
  {
    size_t i, j;
    hkey_t tmp_n;
    Orientation tmp_o;
    for(i = 0, j = len-1; i <= j; i++, j--)
    {
      SWAP(nlist[i], nlist[j], tmp_n);
      tmp_o = olist[i];
      olist[i] = opposite_orientation(olist[j]);
      olist[j] = opposite_orientation(tmp_o);
    }
  }
}

static void naturalise_supernode(hkey_t *nlist, Orientation *olist, size_t len)
{
  // Sort supernode into forward orientation
  if(len == 1)
    olist[0] = forward;
  else if(nlist[0] > nlist[len-1])
    reverse_node_list(nlist, olist, len);
}

// Returns the number of nodes added
// adds no more than `limit`
static size_t fetch_supernode(hkey_t node, Orientation or,
                              hkey_t *nlist, Orientation *olist,
                              size_t limit, const dBGraph *db_graph,
                              char *out_of_space)
{
#ifdef DEBUG_CALLER
  char tmpstr[100];
  binary_kmer_to_str(db_graph_bkmer(db_graph, node), db_graph->kmer_size, tmpstr);
  printf("  fetch %s:%i\n", tmpstr, (int)or);
#endif

  nlist[0] = node;
  olist[0] = or;

  size_t num_nodes = 1;
  Nucleotide nuc;

  const Edges *edges = db_graph->edges;

  while(edges_has_precisely_one_edge(edges[node], or, &nuc))
  {
    #ifdef DEBUG_CALLER
      char tmp[100];
      binary_kmer_to_str(db_graph_bkmer(db_graph, node), db_graph->kmer_size, tmp);
      printf(">%s:%i nuc:%c\n", tmp, or, binary_nuc_to_char(nuc));
    #endif

    db_graph_next_node_orient(db_graph, db_graph_bkmer(db_graph, node), nuc, or,
                              &node, &or);

    if(edges_has_precisely_one_edge(edges[node], rev_orient(or), &nuc))
    {
      if(num_nodes == limit) { *out_of_space = 1; break; }
      if(node == nlist[0]) break; // don't create a loop

      nlist[num_nodes] = node;
      olist[num_nodes] = or;
      num_nodes++;
    }
    else break;
  }

  // char tmp[1000];
  // nodes_to_str(nlist, olist, num_nodes, db_graph->kmer_size, tmp);
  // printf("   supernode: %s\n", tmp);

  return num_nodes;
}

// Remove mark traversed and reset shades
static void supernode_reset(dBGraph *db_graph, CallerSupernode *snode)
{
  db_node_fast_clear_traversed(db_graph, snode->nodes[0]);
  db_node_fast_clear_traversed(db_graph, snode->nodes[snode->num_of_nodes-1]);
}

// Returns 0 on failure, otherwise snode->num_of_nodes
static size_t create_supernode(hkey_t node, Orientation or,
                               CallerSupernode *snode,
                               size_t limit, const dBGraph *db_graph)
{
#ifdef DEBUG_CALLER
  char tmpstr[100];
  binary_kmer_to_str(db_graph_bkmer(db_graph, node), db_graph->kmer_size, tmpstr);
  printf(" create %s:%i\n", tmpstr, (int)or);
#endif

  // extend path
  char out_of_space = 0;
  snode->num_of_nodes = fetch_supernode(node, or, snode->nodes, snode->orients,
                                        limit, db_graph, &out_of_space);

  if(out_of_space)
    return 0;

  naturalise_supernode(snode->nodes, snode->orients, snode->num_of_nodes);

  snode->first_pathpos = NULL;
  snode->num_prev = 0;
  snode->num_next = 0;

  // memset(snode->shades, 0, sizeof(ShadeSet) * snode->num_of_nodes);

  Nucleotide nuc;
  Edges union_edges;
  hkey_t first_node, last_node;
  Orientation first_or, last_or;

  first_node = snode->nodes[0];
  first_or = opposite_orientation(snode->orients[0]);
  last_node = snode->nodes[snode->num_of_nodes-1];
  last_or = snode->orients[snode->num_of_nodes-1];

  // Prev nodes
  union_edges = db_graph->edges[first_node];

  for(nuc = 0; nuc < 4; nuc++) {
    if(edges_has_edge(union_edges, nuc, first_or)) {
      db_graph_next_node_orient(db_graph, db_graph_bkmer(db_graph,first_node),
                                nuc, first_or,
                                snode->prev_nodes + snode->num_prev,
                                snode->prev_orients + snode->num_prev);
      snode->num_prev++;
    }
  }

  // Next nodes
  union_edges = db_graph->edges[last_node];

  for(nuc = 0; nuc < 4; nuc++) {
    if(edges_has_edge(union_edges, nuc, last_or)) {
      db_graph_next_node_orient(db_graph, db_graph_bkmer(db_graph,last_node),
                                nuc, last_or,
                                snode->next_nodes + snode->num_next,
                                snode->next_orients + snode->num_next);
      snode->num_next++;
    }
  }

#ifdef DEBUG_CALLER
  char tmpstr1[100], tmpstr2[100];
  binary_kmer_to_str(db_graph_bkmer(db_graph, first_node), db_graph->kmer_size, tmpstr1);
  binary_kmer_to_str(db_graph_bkmer(db_graph, last_node), db_graph->kmer_size, tmpstr2);
  printf("   ( first %s:%i [%i]; last: %s:%i [%i] )\n",
         tmpstr1, (int)first_or, snode->num_prev,
         tmpstr2, (int)last_or, snode->num_next);
#endif

  return snode->num_of_nodes;
}

// Returns 1 if successfully walked across supernode, 0 otherwise
static boolean walk_supernode(GraphWalker *wlk, CallerSupernode *snode,
                              SuperOrientation snorient)
{
  size_t first = 0, last = snode->num_of_nodes-1, tmp;
  hkey_t first_node, last_node;
  Orientation first_orient, last_orient;
  BinaryKmer first_bkmer, last_bkmer;

  // Only need to traverse the last nodes of a supernode
  if(snorient == forward) {
    SWAP(first, last, tmp);
  }

  first_node = snode->nodes[first];
  first_orient = snode->orients[first];
  last_node = snode->nodes[last];
  last_orient = snode->orients[last];

  // Check if we are happy traversing
  if(wlk->num_paths == 0)
  {
    // We only need to mark nodes as visited if we are 
    if(db_node_has_traversed(wlk->db_graph, first_node, first_orient))
      return false;

    db_node_set_traversed(wlk->db_graph, first_node, first_orient);
  }

  // Force traversal of first and last nodes
  db_graph_oriented_bkmer(wlk->db_graph, first_node, first_orient, first_bkmer);
  graph_traverse_force(wlk, first_node, first_bkmer, first_orient, true);

  if(last > first) {
    db_graph_oriented_bkmer(wlk->db_graph, last_node, last_orient, last_bkmer);
    graph_traverse_force(wlk, last_node, last_bkmer, last_orient, false);
  }

  return true;
}

// Constructs a path of supernodes (SupernodePath)
// returns number of supernodes loaded
static void load_allele_path(hkey_t node, Orientation or,
                             SupernodePath *path,
                             khash_t(supnode_hsh) *snode_hash,
                             GraphWalker *wlk, // walker set to go
                             // these 4 params are tmp memory
                             hkey_t *node_store, Orientation *or_store,
                             CallerSupernode *snode_store,
                             SupernodePathPos *snodepos_store,
                             size_t *node_count_ptr,
                             size_t *snode_count_ptr,
                             size_t *snodepos_count_ptr)
{
  dBGraph *db_graph = wlk->db_graph;

  #ifdef DEBUG_CALLER
    char tmp[100];
    binary_kmer_to_str(db_graph_bkmer(db_graph,node), db_graph->kmer_size, tmp);
    printf(" load_allele_path: %s:%i\n", tmp, or);
  #endif

  size_t node_count = *node_count_ptr;
  size_t snode_count = *snode_count_ptr;
  size_t snodepos_count = *snodepos_count_ptr;

  size_t supindx, kmers_in_path = 0;
  for(supindx = 0; ; supindx++)
  {
    // Find or add supernode
    int hashret;
    khiter_t k = kh_put(supnode_hsh, snode_hash, (uint64_t)node, &hashret);
    CallerSupernode *snode;

    boolean supernode_already_exists = (hashret == 0);

    if(supernode_already_exists)
    {
      // already in supernode hash table
      snode = kh_value(snode_hash, k);
    }
    else
    {
      // add to supernode hash table
      snode = snode_store + snode_count;
      snode->nodes = node_store + node_count;
      snode->orients = or_store + node_count;

      size_t node_buf_space = NODE_BUFFER_SIZE - node_count;
      size_t sn_len = create_supernode(node, or, snode, node_buf_space, db_graph);

      // We fail if we cannot extend this supernode anymore
      // (num of nodes > node_buf_space)
      if(sn_len == 0)
        break;

      kh_value(snode_hash, k) = snode;

      // Add end node
      if(sn_len > 1)
      {
        hkey_t node2 = snode->nodes[node == snode->nodes[0] ? sn_len-1 : 0];
        k = kh_put(supnode_hsh, snode_hash, (uint64_t)node2, &hashret);
        kh_value(snode_hash, k) = snode;
      }

      snode_count++;
      node_count += snode->num_of_nodes;
    }

    SuperOrientation snorient = supernode_get_orientation(snode, node, or);

    if(kmers_in_path + snode->num_of_nodes > MAX_ALLELE_KMERS)
      break;

    // Can we walk along this new supernode?
    // Update graph walker by walking along supernode
    if(!walk_supernode(wlk, snode, snorient)) {
      // Remove mark traversed and reset shades
      supernode_reset(db_graph, snode);
      break;
    }

    kmers_in_path += snode->num_of_nodes;

    // Create SupernodePathPos
    SupernodePathPos *pathpos = snodepos_store + snodepos_count++;
    pathpos->path = path;
    pathpos->pos = supindx;
    pathpos->next = NULL;

    // Add pathpos to front of stack (linkedlist)
    pathpos->next = snode->first_pathpos;
    snode->first_pathpos = pathpos;

    path->supernodes[supindx] = snode;
    path->superorients[supindx] = snorient;

    size_t i, num_edges;
    hkey_t *next_nodes;
    Orientation *next_orients;

    if(snorient == forward) {
      num_edges = snode->num_next;
      next_nodes = snode->next_nodes;
      next_orients = snode->next_orients;
    }
    else {
      num_edges = snode->num_prev;
      next_nodes = snode->prev_nodes;
      next_orients = snode->prev_orients;
    }

    // Get oriented bkmers
    // BinaryKmer next_bkmers[4];
    // for(i = 0; i < num_edges; i++) {
    //   db_graph_oriented_bkmer(db_graph, next_nodes[i], next_orients[i],
    //                           next_bkmers[i]);
    // }

    // if(!graph_traverse_nodes(wlk, num_edges, next_nodes, next_bkmers, next_orients))
    //   return;

    // Get last bases
    Nucleotide next_bases[4];
    for(i = 0; i < num_edges; i++) {
      next_bases[i] = db_node_last_nuc(db_graph_bkmer(db_graph, next_nodes[i]),
                                       next_orients[i], db_graph->kmer_size);
    }

    int nxt_idx = graph_walker_choose(wlk, num_edges, next_bases);
    if(nxt_idx == -1) break;

    // printf("TRAVERSED %zu %i\n", (size_t)node, or);

    node = next_nodes[nxt_idx];
    or = next_orients[nxt_idx];
  }
  // printf("DONE\n");

  *node_count_ptr = node_count;
  *snode_count_ptr = snode_count;
  *snodepos_count_ptr = snodepos_count;
}

//
// Sort Paths
//

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2)
{
  const SupernodePathPos * const pp1 = p1;
  const SupernodePathPos * const pp2 = p2;

  ptrdiff_t cmp = (ptrdiff_t)pp1->pos - pp2->pos;

  if(cmp != 0)
    return cmp;

  int i, len = pp1->pos;
  for(i = 0; i < len; i++)
  {
    cmp = pp1->path->superorients[i] - pp2->path->superorients[i];

    if(cmp != 0)
      return cmp;

    cmp = pp1->path->supernodes[i] - pp2->path->supernodes[i];

    if(cmp != 0)
      return (cmp > 0 ? 1 : -1);
  }

  return 0;
}

#define supernode_pathpos_equal(a,b) (cmp_snpath_pos(a,b) == 0)

static uint64_t supernode_pathpos_hash(SupernodePathPos *spp)
{
  uint32_t hsh, len = spp->pos + 1;
  size_t snode_size = sizeof(CallerSupernode*) * len;
  size_t sorients_size = sizeof(SuperOrientation) * len;

#ifdef CITY_HASH
  // Use Google's CityHash
  hsh = CityHash32((char*)spp->path->supernodes, snode_size);
  hsh ^= CityHash32((char*)spp->path->superorients, sorients_size);
#else
  // Use Bob Jenkin's lookup3
  hsh = hashlittle(spp->path->supernodes, snode_size, 0);
  hsh = hashlittle(spp->path->superorients, sorients_size, hsh);
#endif

  return hsh;
}

// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_INIT(snpps_hsh, SupernodePathPos*, char, 0, supernode_pathpos_hash, supernode_pathpos_equal)

static void remove_snpath_pos_dupes(SupernodePathPos **results, size_t *arrlen,
                                    khash_t(snpps_hsh) *spphash)
{
  kh_clear(snpps_hsh, spphash);

  size_t i, j, len = *arrlen;
  for(i = 0, j = 0; i < len; i++) {
    int hashret;
    kh_put(snpps_hsh, spphash, results[i], &hashret);
    if(hashret != 0) results[j++] = results[i];
  }
  *arrlen = j;
}

static void get_prev_supernode_pos(SupernodePathPos *spp,
                                   CallerSupernode **snode,
                                   SuperOrientation *snorient)
{
  if(spp->pos == 0) {
    *snode = NULL;
    *snorient = forward;
  } else {
    *snode = spp->path->supernodes[spp->pos-1];
    *snorient = spp->path->superorients[spp->pos-1];
  }
}

static char is_bubble_flank(SupernodePathPos **spp_arr, int num)
{
  if(num == 0) return 0;

  int i;
  CallerSupernode *snode0, *snode1;
  Orientation snorient0, snorient1;

  // Check that at least one path has a different first supernode
  snode0 = spp_arr[0]->path->supernodes[0];
  snorient0 = spp_arr[0]->path->superorients[0];

  for(i = 1; i < num; i++) {
    if(snode0 != spp_arr[i]->path->supernodes[0] ||
       snorient0 != spp_arr[i]->path->superorients[0]) {
      break;
    }
  }

  if(i == num) return 0;

  // Check that at least one path has a different last supernode
  get_prev_supernode_pos(spp_arr[0], &snode0, &snorient0);

  for(i = 1; i < num; i++) {
    get_prev_supernode_pos(spp_arr[i], &snode1, &snorient1);
    if(snode0 != snode1 || snorient0 != snorient1) return 1;
  }

  return 0;
}

static void find_bubbles(hkey_t fork_n, Orientation fork_o,
                         dBGraph *db_graph, GraphWalker *wlk,
                         khash_t(supnode_hsh) *snode_hash,
                         khash_t(snpps_hsh) *spp_hash,
                         hkey_t *node_store, Orientation *or_store,
                         SupernodePath *paths,
                         CallerSupernode *snode_store,
                         SupernodePathPos *snodepos_store,
                         gzFile out, size_t *bnum)
{
  size_t snode_count = 0, snodepos_count = 0, node_count = 0;

  // Clear the hash table of supernodes
  kh_clear(supnode_hsh, snode_hash);

  hkey_t nodes[4];
  BinaryKmer bkmers[4];
  Orientation orients[4];
  Edges merged_edges = db_graph->edges[fork_n];
  size_t num_next;

  num_next = db_graph_next_nodes_orient(db_graph,
                                        db_graph_bkmer(db_graph, fork_n),
                                        merged_edges, fork_o,
                                        nodes, bkmers, orients);

#ifdef DEBUG_CALLER
  char tmpstr[100];
  binary_kmer_to_str(db_graph_bkmer(db_graph, fork_n), db_graph->kmer_size, tmpstr);
  printf("fork %s:%i out-degree:%i\n", tmpstr, (int)fork_o, (int)num_next);
#endif
  // loop over alleles, then colours
  size_t i, k, supindx, num_of_paths = 0;
  Colour col, col_loaded = db_graph->ginfo.num_of_colours_loaded;
  for(i = 0; i < num_next; i++)
  {
    for(col = 0; col < col_loaded; col++)
    {
      SupernodePath *path = paths + num_of_paths++;

      // See if we can walk back to pick up shades for this allele/colour
      hkey_t preallele_nodes[MAX_WALK_BACK_KMERS];
      Orientation preallele_or[MAX_WALK_BACK_KMERS];

      Orientation fork_opp = opposite_orientation(fork_o);
      Orientation allele_opp = opposite_orientation(orients[i]);

      // DEV: move this to function graph_walker_init_context
      preallele_nodes[0] = nodes[i];
      preallele_or[0] = orients[i];
      preallele_nodes[1] = fork_n;
      preallele_or[1] = fork_o;

      // Walk backwards over the first kmer of this allele
      graph_walker_init(wlk, db_graph, col, nodes[i], allele_opp);
      db_node_set_traversed(db_graph, wlk->node, wlk->orient);

      // then backwards over the forking kmer (last kmer of the 5p flank)
      BinaryKmer tmpbkmer[4];
      graph_traverse_nodes(wlk, 1, &fork_n, tmpbkmer, &fork_opp);
      db_node_set_traversed(db_graph, wlk->node, wlk->orient);

      k = 2;
      while(k < MAX_WALK_BACK_KMERS && graph_traverse(wlk) &&
            !db_node_has_traversed(db_graph, wlk->node, wlk->orient))
      {
        db_node_set_traversed(db_graph, wlk->node, wlk->orient);
        preallele_nodes[k] = wlk->node;
        preallele_or[k] = opposite_orientation(wlk->orient);
        k++;
      }

      #ifdef DEBUG
        printf("Finished backtracking (%zu nodes)\n", k);
      #endif

      size_t num_prev_kmers = k;

      graph_walker_finish(wlk);

      k--;
      graph_walker_init(wlk, db_graph, col, preallele_nodes[k], preallele_or[k]);

      // Walk back over the kmers (but not onto the first kmer of the allele)
      for(k--; k > 0; k--)
        graph_traverse_nodes(wlk, 1, preallele_nodes+k, tmpbkmer, preallele_or+k);

      // Remove marks on all kmers
      for(k = 0; k < num_prev_kmers; k++)
        db_node_fast_clear_traversed(db_graph, preallele_nodes[k]);

      // Constructs a path of supernodes (SupernodePath)
      load_allele_path(nodes[i], orients[i], path, snode_hash, wlk,
                       node_store, or_store, snode_store, snodepos_store,
                       &node_count, &snode_count, &snodepos_count);

      graph_walker_finish(wlk);

      // Remove mark traversed and reset shades
      for(supindx = 0; supindx < snode_count; supindx++)
        supernode_reset(db_graph, snode_store + supindx);
    }
  }

  hkey_t flank5pe[MAX_FLANK_KMERS];
  Orientation flank5po[MAX_FLANK_KMERS];
  size_t flank5pkmers = 0;

  // Loop over supernodes checking if they are 3p flanks
  for(i = 0; i < snode_count; i++)
  {
    #ifdef DEBUG_CALLER
    char tmpsup[100];
    binary_kmer_to_str(db_graph_bkmer(db_graph, snode_store[i].nodes[0]), db_graph->kmer_size, tmpsup);
    printf("check supernode: %s\n", tmpsup);
    #endif

    SupernodePathPos *pp = snode_store[i].first_pathpos;
    // pp may be null if the node could not be traversed by any path
    // e.g. supernode is also 5p flank and no shades available
    if(pp != NULL && pp->next != NULL)
    {
      // possible 3p flank (i.e. bubble end)
      // there a 4 possible allele paths per colour
      // each path can hit a node at most NUM_OF_SHADES times
      SupernodePathPos *spp_forward[(NUM_OF_SHADES+1) * NUM_OF_COLOURS * 4];
      SupernodePathPos *spp_reverse[(NUM_OF_SHADES+1) * NUM_OF_COLOURS * 4];
      size_t num_forward = 0, num_reverse = 0;

      do
      {
        if(pp->path->superorients[pp->pos] == forward) {
          spp_forward[num_forward++] = pp;
        }
        else {
          spp_reverse[num_reverse++] = pp;
        }
      } while((pp = pp->next) != NULL);

      #ifdef DEBUG_CALLER
      printf(" num_forward: %zu; num_reverse: %zu\n", num_forward, num_reverse);
      #endif

      boolean is_fw_flank = is_bubble_flank(spp_forward, num_forward);
      boolean is_rv_flank = is_bubble_flank(spp_reverse, num_reverse);

      if((is_fw_flank || is_rv_flank) && flank5pkmers == 0)
      {
        // Get 5p flank
        // we can ignore the value of out_of_space
        char out_of_space;
        flank5pkmers = fetch_supernode(fork_n, opposite_orientation(fork_o),
                                       flank5pe, flank5po, MAX_FLANK_KMERS,
                                       db_graph, &out_of_space);

        // Reverse and change orients
        reverse_node_list(flank5pe, flank5po, flank5pkmers);
      }

      // Check if actual bubbles in fw / rv
      if(is_fw_flank)
      {
        remove_snpath_pos_dupes(spp_forward, &num_forward, spp_hash);
        print_bubble(out, *bnum, db_graph,
                     spp_forward, num_forward,
                     flank5pe, flank5po, flank5pkmers);
        (*bnum)++;
      }

      if(is_rv_flank)
      {
        remove_snpath_pos_dupes(spp_reverse, &num_reverse, spp_hash);
        print_bubble(out, *bnum, db_graph,
                     spp_reverse, num_reverse,
                     flank5pe, flank5po, flank5pkmers);
        (*bnum)++;
      }
    }
  }
}

// Returns number of nodes
// if maxlen is zero it is ignored
static size_t suppathpos_to_list(const SupernodePathPos *snodepathpos,
                                 size_t start, size_t suplen,
                                 hkey_t *nlist, Orientation *olist,
                                 size_t maxlen)
{
  SupernodePath *path = snodepathpos->path;

  size_t pos, i, j = 0, end = start+suplen;

  for(pos = start; pos < end && j < maxlen; pos++)
  {
    CallerSupernode *snode = path->supernodes[pos];

    if(path->superorients[pos] == forward)
    {
      size_t len = MIN2(snode->num_of_nodes, maxlen-j);
      memcpy(nlist+j, snode->nodes, len*sizeof(hkey_t));
      memcpy(olist+j, snode->orients, len*sizeof(Orientation));
      j += len;
    }
    else
    {
      for(i = snode->num_of_nodes - 1; j < maxlen; i--) {
        nlist[j] = snode->nodes[i];
        olist[j] = opposite_orientation(snode->orients[i]);
        j++;
        if(i == 0) break;
      }
    }
  }
  return j;
}


static void print_bubble(gzFile out, size_t bnum, dBGraph *db_graph,
                         SupernodePathPos **spp_arr, size_t num_of_paths,
                         hkey_t *flank5pe, Orientation *flank5po,
                         size_t flank5pkmers)
{
  // tmp variables
  hkey_t tmp_e[MAX_ALLELE_KMERS];
  Orientation tmp_o[MAX_ALLELE_KMERS];
  size_t i, num_kmers;

  // 5p flank
  gzprintf(out, ">var_%lu_5p_flank length=%zu\n", bnum, flank5pkmers);
  print_branch(flank5pe, flank5po, flank5pkmers, 1, db_graph, out);

  // 3p flank
  num_kmers = suppathpos_to_list(spp_arr[0], spp_arr[0]->pos, 1,
                                 tmp_e, tmp_o, MAX_FLANK_KMERS);
  gzprintf(out, ">var_%lu_3p_flank length=%zu\n", bnum, num_kmers);
  print_branch(tmp_e, tmp_o, num_kmers, 0, db_graph, out);

  // Print alleles
  for(i = 0; i < num_of_paths; i++)
  {
    num_kmers = suppathpos_to_list(spp_arr[i], 0, spp_arr[i]->pos,
                                   tmp_e, tmp_o, MAX_ALLELE_KMERS);
    gzprintf(out, ">var_%zu_branch_%zu length=%zu\n", bnum, i, num_kmers);
    print_branch(tmp_e, tmp_o, num_kmers, 0, db_graph, out);
  }

  gzputc(out, '\n');
}


static void shaded_bubble_caller(hkey_t node,
                                 dBGraph *db_graph, GraphWalker *wlk,
                                 khash_t(supnode_hsh) *snode_hash,
                                 khash_t(snpps_hsh) *spp_hash,
                                 hkey_t *node_store, Orientation *or_store,
                                 SupernodePath *paths,
                                 CallerSupernode *snode_store,
                                 SupernodePathPos *snodepos_store,
                                 gzFile out, size_t *bnum)
{
  Edges edges = db_graph->edges[node];

  if(edges_get_outdegree(edges, forward) > 1)
  {
    find_bubbles(node, forward, db_graph, wlk,
                 snode_hash, spp_hash, node_store, or_store,
                 paths, snode_store, snodepos_store,
                 out, bnum);
  }
  if(edges_get_outdegree(edges, reverse) > 1)
  {
    find_bubbles(node, reverse, db_graph, wlk,
                 snode_hash, spp_hash, node_store, or_store,
                 paths, snode_store, snodepos_store,
                 out, bnum);
  }
}


void invoke_shaded_bubble_caller(dBGraph *db_graph, const char* out_file)
                                 // const CmdLine *cmd)
{
  // Open output file
  gzFile out = gzopen(out_file, "w");

  if(out == NULL)
    die("Cannot open paths bubble caller output file: %s", out_file);

  // Print header
  print_calling_header(db_graph, out, out_file);//, cmd);

  // Array to re-use
  int cols4 = 4 * NUM_OF_COLOURS;

  SupernodePath *paths = malloc(cols4 * sizeof(SupernodePath));
  CallerSupernode *snodes = malloc(NODE_BUFFER_SIZE * sizeof(CallerSupernode));
  SupernodePathPos *snodespos = malloc(NODE_BUFFER_SIZE * sizeof(SupernodePathPos));

  khash_t(supnode_hsh) *snode_hash = kh_init(supnode_hsh);
  khash_t(snpps_hsh) *spp_hash = kh_init(snpps_hsh);

  hkey_t *nlist = malloc(NODE_BUFFER_SIZE * sizeof(hkey_t));
  Orientation *olist = malloc(NODE_BUFFER_SIZE * sizeof(Orientation));

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  size_t num_of_bubbles = 0;

  // Iterate over nodes in the hash table
  // for each one, check if it has degree > 1
  HASH_TRAVERSE(&db_graph->ht, shaded_bubble_caller,
                db_graph, &wlk,
                snode_hash, spp_hash, nlist, olist,
                paths, snodes, snodespos,
                out, &num_of_bubbles);

  // Test calling from a particular node
  // hkey_t tmp = hash_table_find_str("ACCATCCCCTAGGGCTTCTCCACACTCAACT", db_graph);
  // find_bubbles(tmp, reverse, db_graph,
  //                snode_hash, spp_hash, nlist, olist,
  //                paths, snodes, snodespos,
  //                out, &num_of_bubbles);

  gzclose(out);

  message("%zu bubbles called with Paths-Bubble-Caller\n", num_of_bubbles);
  message("  saved to: %s\n", out_file);

  graph_walker_dealloc(&wlk);

  free(paths);
  free(snodes);
  free(snodespos);
  kh_destroy(supnode_hsh, snode_hash);
  kh_destroy(snpps_hsh, spp_hash);
  free(nlist);
  free(olist);
}
