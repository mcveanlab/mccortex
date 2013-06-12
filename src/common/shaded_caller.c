#include "global.h"

#include <stddef.h> // defines ptrdiff_t
#include <time.h> // printing datetime
#include <pthread.h> // multithreading

#include "khash.h"

// cortex_var headers
#include "shaded_caller.h"
#include "db_graph.h"
#include "db_node.h"
#include "util.h"
#include "file_util.h"
#include "file_reader.h"
#include "seq_reader.h"
#include "binary_paths.h"
#include "graph_walker.h"
#include "caller_supernode.h"

// Hash functions
#include "lookup3.h"
#include "city.h"

#ifdef DEBUG
#define DEBUG_CALLER 1
#endif

// Print absolute path to a file
static void print_filepath_abs(gzFile *out, const char *name, const char *file)
{
  char absolute_path[PATH_MAX + 1];
  char *abs_path = realpath(file, absolute_path);

  if(abs_path == NULL)
    warn("Cannot get absolute path: %s\n", file);

  gzprintf(out, "##%s=%s\n", name, abs_path);
}

static void print_calling_header(const dBGraph *db_graph, gzFile *out,
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
              "compiled_colours=%i>\n",
          VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION,
          db_graph->kmer_size, NUM_OF_COLOURS);

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

  print_filepath_abs(out, "ctxBubblesFile", out_file);
  gzprintf(out, "##ctxKmerSize=%u\n", db_graph->kmer_size);

  // Print ref colour if there is one
  // if(cmd->ref_colour != -1)
  //   gzprintf(out, "##ctx_refcol=%i\n", cmd->ref_colour);

  // Print colours we're calling in
  gzprintf(out, "##ctxNumCallingUsedInColours=%i\n", ginfo->num_of_colours_loaded);

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
                         char print_first_kmer, const dBGraph *db_graph,
                         gzFile out)
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
    nuc = db_node_last_nuc(db_node_bkmer(db_graph, nodes[i]),
                           orients[i], kmer_size);
    gzputc(out, binary_nuc_to_char(nuc));
  }

  gzputc(out, '\n');

  Colour col, cols_loaded = db_graph->ginfo.num_of_colours_loaded;

  // Print (fake) covgs
  if(len > 0) {
    for(col = 0; col < cols_loaded; col++) {
      gzprintf(out, "1");
      for(i = 1; i < len; i++) {
        gzprintf(out, " 1");
      }
      gzputc(out, '\n');
    }
  } else gzputc(out, '\n');
}


// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_MAP_INIT_INT64(supnode_hsh, CallerSupernode*)

static void print_bubble(gzFile out, size_t bnum, const dBGraph *db_graph,
                         SupernodePathPos **spp_arr, size_t num_of_paths,
                         hkey_t *flank5pe, Orientation *flank5po,
                         size_t flank5pkmers);

// Returns 1 if successfully walked across supernode, 0 otherwise
static void walk_supernode_end(GraphWalker *wlk, CallerSupernode *snode,
                               SuperOrientation snorient)
{
  // Only need to traverse the first and last nodes of a supernode
  size_t last = snode->num_of_nodes-1;
  hkey_t last_node;
  Orientation last_orient;
  BinaryKmer last_bkmer;

  if(last > 0) {
    if(snorient == forward) {
      last_node = snode->nodes[last];
      last_orient = snode->orients[last];
    } else {
      last_node = snode->nodes[0];
      last_orient = rev_orient(snode->orients[0]);
    }

    db_graph_oriented_bkmer(wlk->db_graph, last_node, last_orient, last_bkmer);
    graph_traverse_force_jump(wlk, last_node, last_bkmer, false);
  }
}


// Constructs a path of supernodes (SupernodePath)
// returns number of supernodes loaded
static void load_allele_path(hkey_t node, Orientation or,
                             SupernodePath *path,
                             khash_t(supnode_hsh) *snode_hash,
                             GraphWalker *wlk, // walker set to go
                             uint64_t *visited,
                             // these 4 params are tmp memory
                             hkey_t *node_store, Orientation *or_store,
                             CallerSupernode *snode_store,
                             SupernodePathPos *snodepos_store,
                             size_t *node_count_ptr,
                             size_t *snode_count_ptr,
                             size_t *snodepos_count_ptr)
{
  const dBGraph *db_graph = wlk->db_graph;

  #ifdef DEBUG_CALLER
    char tmp[100];
    binary_kmer_to_str(db_node_bkmer(db_graph,node), db_graph->kmer_size, tmp);
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

    // Walk along supernode (always succedes)
    walk_supernode_end(wlk, snode, snorient);

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

    // Get last bases
    Nucleotide next_bases[4];
    for(i = 0; i < num_edges; i++) {
      next_bases[i] = db_node_last_nuc(db_node_bkmer(db_graph, next_nodes[i]),
                                       next_orients[i], db_graph->kmer_size);
    }

    int nxt_idx = graph_walker_choose(wlk, num_edges, next_nodes, next_bases);
    if(nxt_idx == -1) break;

    node = next_nodes[nxt_idx];
    or = next_orients[nxt_idx];
    Nucleotide base = next_bases[nxt_idx];

    // Check if we are happy traversing
    // if(wlk->num_paths == 0) {
      // We only need to mark nodes as visited if we have no paths
      if(db_node_has_traversed(visited, node, or)) break;
      db_node_set_traversed(visited, node, or);
    // }

    graph_traverse_force(wlk, node, base, num_edges > 1);
  }
  // printf("DONE\n");

  *node_count_ptr = node_count;
  *snode_count_ptr = snode_count;
  *snodepos_count_ptr = snodepos_count;
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
                         const dBGraph *db_graph, GraphWalker *wlk,
                         uint64_t *visited,
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
  size_t i, num_next;

  num_next = db_graph_next_nodes_orient(db_graph,
                                        db_node_bkmer(db_graph, fork_n),
                                        db_node_edges(db_graph, fork_n), fork_o,
                                        nodes, bkmers);

  for(i = 0; i < num_next; i++) {
    orients[i] = db_node_get_orientation(db_node_bkmer(db_graph, nodes[i]),
                                         bkmers[i]);
  }

  #ifdef DEBUG_CALLER
    char tmpstr[100];
    binary_kmer_to_str(db_node_bkmer(db_graph, fork_n), db_graph->kmer_size, tmpstr);
    printf("fork %s:%i out-degree:%i\n", tmpstr, (int)fork_o, (int)num_next);
  #endif

  // loop over alleles, then colours
  size_t supindx, num_of_paths = 0;
  Colour colour, colours_loaded = db_graph->ginfo.num_of_colours_loaded;
  for(i = 0; i < num_next; i++)
  {
    for(colour = 0; colour < colours_loaded; colour++)
    {
      SupernodePath *path = paths + num_of_paths++;

      // See if we can walk back to pick up paths for this allele/colour
      graph_init_context(wlk, db_graph, visited, colour, nodes[i], orients[i]);

      // Constructs a path of supernodes (SupernodePath)
      load_allele_path(nodes[i], orients[i], path, snode_hash, wlk,
                       visited,
                       node_store, or_store, snode_store, snodepos_store,
                       &node_count, &snode_count, &snodepos_count);

      graph_walker_finish(wlk);

      // Remove mark traversed and reset shades
      for(supindx = 0; supindx < snode_count; supindx++)
      {
        CallerSupernode *snode = snode_store + supindx;
        db_node_fast_clear_traversed(visited, snode->nodes[0]);
        db_node_fast_clear_traversed(visited, snode->nodes[snode->num_of_nodes-1]);
      }
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
    binary_kmer_to_str(db_node_bkmer(db_graph, snode_store[i].nodes[0]), db_graph->kmer_size, tmpsup);
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
      SupernodePathPos *spp_forward[MAX_ALLELE_KMERS * NUM_OF_COLOURS * 4];
      SupernodePathPos *spp_reverse[MAX_ALLELE_KMERS * NUM_OF_COLOURS * 4];
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
        boolean out_of_space;
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


static void print_bubble(gzFile out, size_t bnum, const dBGraph *db_graph,
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


struct caller_region_t
{
  const dBGraph *db_graph;
  const gzFile out;
  const uint64_t start_hkey, end_hkey;
  size_t num_of_bubbles;
};

void* shaded_bubble_caller(void *args)
{
  struct caller_region_t *tdata = (struct caller_region_t *)args;
  const dBGraph *db_graph = tdata->db_graph;

  // Arrays to re-use
  int cols4 = 4 * NUM_OF_COLOURS;

  size_t node_bits = round_bits_to_words64(db_graph->ht.capacity);
  uint64_t *visited = calloc(2*node_bits, sizeof(uint64_t));
  // uint64_t *pickedpaths = calloc(node_bits, sizeof(uint64_t));

  SupernodePath *snode_paths = malloc(cols4 * sizeof(SupernodePath));
  CallerSupernode *snode_store = malloc(NODE_BUFFER_SIZE * sizeof(CallerSupernode));
  SupernodePathPos *snodepos_store = malloc(NODE_BUFFER_SIZE * sizeof(SupernodePathPos));

  khash_t(supnode_hsh) *snode_hash = kh_init(supnode_hsh);
  khash_t(snpps_hsh) *spp_hash = kh_init(snpps_hsh);

  hkey_t *node_store = malloc(NODE_BUFFER_SIZE * sizeof(hkey_t));
  Orientation *or_store = malloc(NODE_BUFFER_SIZE * sizeof(Orientation));

  if(visited == NULL || // pickedpaths == NULL ||
     snode_paths == NULL || snode_store == NULL || snodepos_store == NULL ||
     snode_hash == NULL || spp_hash == NULL ||
     node_store == NULL || or_store == NULL)
  {
    die("Out of memory");
  }

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  BinaryKmer *table = db_graph->ht.table;
  BinaryKmer *ptr = table + tdata->start_hkey;
  BinaryKmer *end = table + tdata->end_hkey;

  for(; ptr < end; ptr++) {
    if(HASH_ENTRY_ASSIGNED(*ptr)) {
      hkey_t node = ptr - table;
      Edges edges = db_graph->edges[node];
      if(edges_get_outdegree(edges, forward) > 1) {
        find_bubbles(node, forward, db_graph, &wlk, visited, // pickedpaths,
                     snode_hash, spp_hash, node_store, or_store,
                     snode_paths, snode_store, snodepos_store,
                     tdata->out, &tdata->num_of_bubbles);
      }
      if(edges_get_outdegree(edges, reverse) > 1) {
        find_bubbles(node, reverse, db_graph, &wlk, visited, // pickedpaths,
                     snode_hash, spp_hash, node_store, or_store,
                     snode_paths, snode_store, snodepos_store,
                     tdata->out, &tdata->num_of_bubbles);
      }
    }
  }

  graph_walker_dealloc(&wlk);

  free(visited);
  // free(pickedpaths);

  free(snode_paths);
  free(snode_store);
  free(snodepos_store);
  kh_destroy(supnode_hsh, snode_hash);
  kh_destroy(snpps_hsh, spp_hash);
  free(node_store);
  free(or_store);

  return NULL;
}

void invoke_shaded_bubble_caller(const dBGraph *db_graph, const char* out_file,
                                 int num_threads)
{
  // Open output file
  gzFile out = gzopen(out_file, "w");

  if(out == NULL)
    die("Cannot open paths bubble caller output file: %s", out_file);

  // Print header
  print_calling_header(db_graph, out, out_file);

  size_t num_of_bubbles = 0;

  if(num_threads <= 1)
  {
    struct caller_region_t tdata = {db_graph, out, 0, db_graph->ht.capacity, 0};
    shaded_bubble_caller((void*)&tdata);
    num_of_bubbles += tdata.num_of_bubbles;
  }
  else
  {
    pthread_t threads[num_threads];
    struct caller_region_t tdata[num_threads];
    char *paths[num_threads];
    int rc, i;

    StrBuf *tmppath = strbuf_new();

    uint64_t capacity = db_graph->ht.capacity, start = 0, end;
    gzFile tmpout;

    for(i = 0; i < num_threads; i++, start = end)
    {
      end = i == num_threads-1 ? capacity : start + capacity/num_threads;
      strbuf_set(tmppath, out_file);
      strbuf_sprintf(tmppath, ".%i", i);
      paths[i] = strbuf_dup(tmppath);
      tmpout = gzopen(paths[i], "w");
      if(tmpout == NULL) die("Cannot write to tmp output: %s", paths[i]);

      struct caller_region_t tmptdata = {db_graph, tmpout, start, end, 0};
      memcpy(tdata+i, &tmptdata, sizeof(struct caller_region_t));

      rc = pthread_create(threads+i, NULL, shaded_bubble_caller, (void*)(tdata+i));
      if(rc != 0) die("Creating thread failed");
    }

    /* wait for all threads to complete */
    for(i = 0; i < num_threads; i++) {
      rc = pthread_join(threads[i], NULL);
      if(rc != 0) die("Joining thread failed");
      gzclose(tdata[i].out);

      // merge file
      gzFile in = gzopen(paths[i], "r");
      if(in == NULL) die("Cannot merge file: %s", paths[i]);
      #define MEGABYTE (1<<20)
      char data[MEGABYTE];
      int len;
      while((len = gzread(in, data, MEGABYTE)) > 0) gzwrite(out, data, len);
      gzclose(in);
      unlink(paths[i]);

      free(paths[i]);
      num_of_bubbles += tdata[i].num_of_bubbles;
    }
  
    strbuf_free(tmppath);
  }

  message("%zu bubbles called with Paths-Bubble-Caller\n", num_of_bubbles);
  message("  saved to: %s\n", out_file);

  gzclose(out);
}
