#include "global.h"

#include <stddef.h> // defines ptrdiff_t
#include <time.h> // printing datetime
#include <pthread.h> // multithreading

#include "khash.h"

// cortex_var headers
#include "bubble_caller.h"
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
  gzprintf(out, "##ctxVersion=<version=%s,MAXK=%i>\n",
           CTXVERSIONSTR, MAX_KMER_SIZE);

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
  gzprintf(out, "##ctxNumCallingUsedInColours=%i\n", db_graph->num_of_cols_used);

  StrBuf *sample_name = strbuf_new();

  // Print sample names
  Colour col;
  for(col = 0; col < db_graph->num_of_cols_used; col++)
  {
    const GraphInfo *ginfo = db_graph->ginfo + col;
    const ErrorCleaning *ec = &ginfo->cleaning;

    // Find and replace double quotes with single quotes
    char *sname = ginfo->sample_name.buff;

    if(strcmp(sname, "undefined") == 0 || strchr(sname, '\t') != NULL ||
       strchr(sname, ' ') != NULL || strchr(sname, '\r') != NULL ||
       strchr(sname, '\n') != NULL)
    {
      strbuf_reset(sample_name);
      strbuf_sprintf(sample_name, "sample%u", col);
    }
    else {
      strbuf_set(sample_name, ginfo->sample_name.buff);
    }

    gzprintf(out, "##SAMPLE=<ID=%s,name=\"%s\",colour=%i,"
                  "meanreadlen=%zu,totalseqloaded=%zu,"
                  "seqerror=%Lf,tipclipped=%s,removelowcovgsupernodes=%u,"
                  "removelowcovgkmer=%u,cleanedagainstgraph=%s>\n",
             sample_name->buff, ginfo->sample_name.buff, col,
             (size_t)ginfo->mean_read_length, (size_t)ginfo->total_sequence,
             ginfo->seq_err,
             ec->tip_clipping ? "yes" : "no", ec->remv_low_cov_sups_thresh,
             ec->remv_low_cov_nodes_thresh,
             ec->cleaned_against_graph_name.buff);
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
}


// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_MAP_INIT_INT64(supnode_hsh, CallerSupernode*)

static void print_bubble(gzFile out, size_t bnum, const dBGraph *db_graph,
                         SupernodePathPos **spp_arr, size_t num_of_paths,
                         hkey_t *flank5pe, Orientation *flank5po,
                         size_t flank5pkmers, size_t threadid);

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
    if(snorient == FORWARD) {
      last_node = snode->nodes[last];
      last_orient = snode->orients[last];
    } else {
      last_node = snode->nodes[0];
      last_orient = rev_orient(snode->orients[0]);
    }

    db_graph_oriented_bkmer(wlk->db_graph, last_node, last_orient, last_bkmer);
    graph_traverse_force_jump(wlk, last_node, last_bkmer, false);
    // don't need counter paths here (we're at the end of a supernode)
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
    char tmp[MAX_KMER_SIZE+1];
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
    CallerSupernode *snode;
    int hashret;
    khiter_t k = kh_put(supnode_hsh, snode_hash, (uint64_t)node, &hashret);
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

      size_t node_buf_space = NODE_BUFSIZE(db_graph->num_of_cols) - node_count;
      size_t sn_len = caller_supernode_create(node, or, snode, node_buf_space,
                                              db_graph);

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

    if(snorient == FORWARD) {
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

    Nucleotide lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);
    graph_traverse_force(wlk, node, base, num_edges > 1);
    graph_walker_node_add_counter_paths(wlk, wlk->node, wlk->orient, lost_nuc);
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
    *snorient = FORWARD;
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
                         gzFile out, size_t *bnum, size_t threadid)
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
    char tmpstr[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_bkmer(db_graph, fork_n), db_graph->kmer_size, tmpstr);
    printf("fork %s:%i out-degree:%i\n", tmpstr, (int)fork_o, (int)num_next);
  #endif

  // loop over alleles, then colours
  size_t supindx, num_of_paths = 0;
  Colour colour, colours_loaded = db_graph->num_of_cols_used;
  Nucleotide next_nuc;

  for(i = 0; i < num_next; i++)
  {
    next_nuc = binary_kmer_last_nuc(bkmers[i]);

    for(colour = 0; colour < colours_loaded; colour++)
    {
      // if(db_node_has_col(db_graph, nodes[i], colour)) {
        SupernodePath *path = paths + num_of_paths++;

        graph_walker_init(wlk, db_graph, colour, fork_n, fork_o);
        Nucleotide lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

        // graph_walker_init_context(wlk, db_graph, visited, colour, fork_n, fork_o);
        graph_traverse_force(wlk, nodes[i], next_nuc, true);
        graph_walker_node_add_counter_paths(wlk, wlk->node, wlk->orient, lost_nuc);

        // Constructs a path of supernodes (SupernodePath)
        load_allele_path(nodes[i], orients[i], path, snode_hash, wlk, visited,
                         node_store, or_store, snode_store, snodepos_store,
                         &node_count, &snode_count, &snodepos_count);

        graph_walker_finish(wlk);

        // Remove mark traversed and reset shades
        for(supindx = 0; supindx < snode_count; supindx++)
        {
          CallerSupernode *snode = snode_store + supindx;
          size_t last = snode->num_of_nodes-1;
          db_node_fast_clear_traversed(visited, snode->nodes[0]);
          db_node_fast_clear_traversed(visited, snode->nodes[last]);
        }
      // }
    }
  }

  hkey_t flank5pe[MAX_FLANK_KMERS];
  Orientation flank5po[MAX_FLANK_KMERS];
  size_t flank5pkmers = 0;

  // Loop over supernodes checking if they are 3p flanks
  for(i = 0; i < snode_count; i++)
  {
    #ifdef DEBUG_CALLER
      char tmpsup[MAX_KMER_SIZE+1];
      ConstBinaryKmerPtr bptr = db_node_bkmer(db_graph, snode_store[i].nodes[0]);
      binary_kmer_to_str(bptr, db_graph->kmer_size, tmpsup);
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
      SupernodePathPos *spp_forward[NODE_BUFSIZE(db_graph->num_of_cols)];
      SupernodePathPos *spp_reverse[NODE_BUFSIZE(db_graph->num_of_cols)];
      size_t num_forward = 0, num_reverse = 0;

      do
      {
        if(pp->path->superorients[pp->pos] == FORWARD) {
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
        flank5pkmers = supernode_traverse(fork_n, opposite_orientation(fork_o),
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
                     flank5pe, flank5po, flank5pkmers, threadid);
        (*bnum)++;
      }

      if(is_rv_flank)
      {
        remove_snpath_pos_dupes(spp_reverse, &num_reverse, spp_hash);
        print_bubble(out, *bnum, db_graph,
                     spp_reverse, num_reverse,
                     flank5pe, flank5po, flank5pkmers, threadid);
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

    if(path->superorients[pos] == FORWARD)
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
                         size_t flank5pkmers, size_t threadid)
{
  // tmp variables
  hkey_t tmp_e[MAX_ALLELE_KMERS];
  Orientation tmp_o[MAX_ALLELE_KMERS];
  size_t i, num_kmers;

  // 5p flank
  gzprintf(out, ">var_%zu.%lu_5p_flank length=%zu\n", threadid, bnum, flank5pkmers);
  print_branch(flank5pe, flank5po, flank5pkmers, 1, db_graph, out);

  // 3p flank
  num_kmers = suppathpos_to_list(spp_arr[0], spp_arr[0]->pos, 1,
                                 tmp_e, tmp_o, MAX_FLANK_KMERS);
  gzprintf(out, ">var_%zu.%lu_3p_flank length=%zu\n", threadid, bnum, num_kmers);
  print_branch(tmp_e, tmp_o, num_kmers, 0, db_graph, out);

  // Print alleles
  for(i = 0; i < num_of_paths; i++)
  {
    num_kmers = suppathpos_to_list(spp_arr[i], 0, spp_arr[i]->pos,
                                   tmp_e, tmp_o, MAX_ALLELE_KMERS);
    gzprintf(out, ">var_%zu.%zu_branch_%zu length=%zu\n", threadid, bnum, i, num_kmers);
    print_branch(tmp_e, tmp_o, num_kmers, 0, db_graph, out);
  }

  gzputc(out, '\n');
}


struct caller_region_t
{
  const dBGraph *db_graph;
  const gzFile out;
  const uint64_t start_hkey, end_hkey;
  const size_t threadid;
  size_t num_of_bubbles;
};

void* bubble_caller(void *args)
{
  struct caller_region_t *tdata = (struct caller_region_t *)args;
  const dBGraph *db_graph = tdata->db_graph;

  // Arrays to re-use
  int cols4 = 4 * db_graph->num_of_cols;

  size_t node_bits = round_bits_to_words64(db_graph->ht.capacity);
  uint64_t *visited = calloc(2*node_bits, sizeof(uint64_t));

  size_t buf_size = NODE_BUFSIZE(db_graph->num_of_cols);

  SupernodePath *snode_paths = malloc(cols4 * sizeof(SupernodePath));
  CallerSupernode *snode_store = malloc(buf_size * sizeof(CallerSupernode));
  SupernodePathPos *snodepos_store = malloc(buf_size * sizeof(SupernodePathPos));

  khash_t(supnode_hsh) *snode_hash = kh_init(supnode_hsh);
  khash_t(snpps_hsh) *spp_hash = kh_init(snpps_hsh);

  hkey_t *node_store = malloc(buf_size * sizeof(hkey_t));
  Orientation *or_store = malloc(buf_size * sizeof(Orientation));

  if(visited == NULL || snode_paths == NULL ||
     snode_store == NULL || snodepos_store == NULL ||
     snode_hash == NULL || spp_hash == NULL ||
     node_store == NULL || or_store == NULL)
  {
    die("Out of memory");
  }

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  // BinaryKmer tmpkmer;
  // hkey_t node;
  // binary_kmer_from_str("CTAAAACTGGAACCCAAAACAATGGGAGTGA", db_graph->kmer_size, tmpkmer);
  // node = hash_table_find(&db_graph->ht, tmpkmer);
  // find_bubbles(node, forward, db_graph, &wlk, visited,
  //              snode_hash, spp_hash, node_store, or_store,
  //              snode_paths, snode_store, snodepos_store,
  //              tdata->out, &tdata->num_of_bubbles);

  // binary_kmer_from_str("ACCCAAAACAATGGGAGTGATGTGCTAAAAC", db_graph->kmer_size, tmpkmer);
  // node = hash_table_find(&db_graph->ht, tmpkmer);
  // find_bubbles(node, REVERSE, db_graph, &wlk, visited,
  //              snode_hash, spp_hash, node_store, or_store,
  //              snode_paths, snode_store, snodepos_store,
  //              tdata->out, &tdata->num_of_bubbles);

  
  BinaryKmer *table = db_graph->ht.table;
  BinaryKmer *ptr = table + tdata->start_hkey;
  BinaryKmer *end = table + tdata->end_hkey;

  for(; ptr < end; ptr++) {
    if(HASH_ENTRY_ASSIGNED(*ptr)) {
      hkey_t node = ptr - table;
      Edges edges = db_graph->edges[node];
      if(edges_get_outdegree(edges, FORWARD) > 1) {
        find_bubbles(node, FORWARD, db_graph, &wlk, visited,
                     snode_hash, spp_hash, node_store, or_store,
                     snode_paths, snode_store, snodepos_store,
                     tdata->out, &tdata->num_of_bubbles, tdata->threadid);
      }
      if(edges_get_outdegree(edges, REVERSE) > 1) {
        find_bubbles(node, REVERSE, db_graph, &wlk, visited,
                     snode_hash, spp_hash, node_store, or_store,
                     snode_paths, snode_store, snodepos_store,
                     tdata->out, &tdata->num_of_bubbles, tdata->threadid);
      }
    }
  }
  

  graph_walker_dealloc(&wlk);

  free(visited);
  free(snode_paths);
  free(snode_store);
  free(snodepos_store);
  kh_destroy(supnode_hsh, snode_hash);
  kh_destroy(snpps_hsh, spp_hash);
  free(node_store);
  free(or_store);

  return NULL;
}

void invoke_bubble_caller(const dBGraph *db_graph, const char* out_file,
                                 int num_threads, char **tmp_paths)
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
    struct caller_region_t tdata = {db_graph, out, 0, db_graph->ht.capacity, 0, 0};
    bubble_caller((void*)&tdata);
    num_of_bubbles += tdata.num_of_bubbles;
  }
  else
  {
    pthread_t threads[num_threads];
    pthread_attr_t thread_attr;
    struct caller_region_t tdata[num_threads];

    uint64_t capacity = db_graph->ht.capacity, start = 0, end;
    gzFile tmpout;
    int i;

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&thread_attr);
    pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

    for(i = 0; i < num_threads; i++, start = end)
    {
      end = i == num_threads-1 ? capacity : start + capacity/num_threads;

      tmpout = gzopen(tmp_paths[i], "w");
      if(tmpout == NULL) die("Cannot write to tmp output: %s", tmp_paths[i]);

      struct caller_region_t tmptdata = {db_graph, tmpout, start, end, i, 0};
      memcpy(tdata+i, &tmptdata, sizeof(struct caller_region_t));

      int rc = pthread_create(threads+i, &thread_attr, bubble_caller,
                              (void*)(tdata+i));

      if(rc != 0) die("Creating thread failed");
    }

    /* wait for all threads to complete */
    pthread_attr_destroy(&thread_attr);

    for(i = 0; i < num_threads; i++)
    {
      int rc = pthread_join(threads[i], NULL);
      if(rc != 0) die("Joining thread failed");
      gzclose(tdata[i].out);

      // merge file
      gzFile in = gzopen(tmp_paths[i], "r");
      if(in == NULL) die("Cannot merge file: %s", tmp_paths[i]);
      #define MEGABYTE (1<<20)
      char data[MEGABYTE];
      int len;
      while((len = gzread(in, data, MEGABYTE)) > 0) gzwrite(out, data, len);
      gzclose(in);

      num_of_bubbles += tdata[i].num_of_bubbles;
    }
  }

  char num_bubbles_str[100];
  ulong_to_str(num_of_bubbles, num_bubbles_str);
  message("%s bubbles called with Paths-Bubble-Caller\n", num_bubbles_str);
  message("  saved to: %s\n", out_file);

  gzclose(out);
}
