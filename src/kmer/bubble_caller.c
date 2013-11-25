#include "global.h"

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
#include "cmd.h"
#include "seq_reader.h"
#include "path_store.h"
#include "graph_walker.h"
#include "caller_supernode.h"
#include "supernode.h"

// Hash functions
// #include "lookup3.h"
// #include "city.h"

#ifdef CTXVERBOSE
#define DEBUG_CALLER 1
#endif

// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_MAP_INIT_INT64(supnode_hsh, CallerSupernode*)


// Actually hash map of hkey_t -> CallerSupernode*
// the first and last node of a supernode are mapped to the CallerSupernode
KHASH_INIT(snpps_hsh, SupernodePathPos*, char, 0, supernode_pathpos_hash, supernode_pathpos_equal)

// Print absolute path to a file
static void print_filepath_abs(gzFile out, const char *name, const char *file)
{
  char absolute_path[PATH_MAX + 1];
  char *abs_path = realpath(file, absolute_path);

  if(abs_path == NULL)
    warn("Cannot get absolute path: %s\n", file);

  gzprintf(out, "##%s=%s\n", name, abs_path);
}

static void print_calling_header(const dBGraph *db_graph, gzFile out,
                                 const char* out_file, const CmdArgs *args)
{
  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));
  char cwd[PATH_MAX + 1];

  gzprintf(out, "##fileformat=CTXv1.0\n");

  // Cortex details
  gzprintf(out, "##ctxCmd=%s\n", args->cmdline);
  if(futil_get_current_dir(cwd) != NULL) gzprintf(out, "##ctxCwd=%s\n", cwd);

  gzprintf(out, "##ctxDate=%s\n", datestr);
  gzprintf(out, "##ctxVersion=<version=%s,MAXK=%i>\n",
           CTXVERSIONSTR, MAX_KMER_SIZE);

  print_filepath_abs(out, "ctxBubblesFile", out_file);
  gzprintf(out, "##ctxKmerSize=%u\n", db_graph->kmer_size);

  // Print colours we're calling in
  gzprintf(out, "##ctxNumColoursUsedInCalling=%i\n", db_graph->num_of_cols_used);

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
      strbuf_sprintf(sample_name, "sample%zu", col);
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
             ec->intersection_name.buff);
  }

  strbuf_free(sample_name);
}

static void print_branch(hkey_t *nodes, Orientation *orients, size_t len,
                         boolean print_first_kmer, const dBGraph *db_graph,
                         gzFile out)
{
  size_t i = print_first_kmer, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;

  if(print_first_kmer) {
    char tmp[MAX_KMER_SIZE+1];
    bkmer = db_graph_oriented_bkmer(db_graph, nodes[0], orients[0]);
    binary_kmer_to_str(bkmer, kmer_size, tmp);
    gzputs(out, tmp);
  }

  // i = 1 if print_first_kmer, otherwise 0
  for(; i < len; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i]);
    nuc = db_node_last_nuc(bkmer, orients[i], kmer_size);
    gzputc(out, binary_nuc_to_char(nuc));
  }

  gzputc(out, '\n');
}

// Returns number of nodes
// if maxlen is zero it is ignored
static size_t suppathpos_to_list(const SupernodePathPos *snodepathpos,
                                 size_t start, size_t suplen,
                                 hkey_t *nlist, Orientation *olist,
                                 size_t maxlen)
{
  SupernodePath *path = snodepathpos->path;
  CallerSupernode *snode;
  hkey_t *nodes;
  Orientation *orients;

  size_t pos, i, j = 0, end = start+suplen, len;

  for(pos = start; pos < end && j < maxlen; pos++)
  {
    snode = path->supernodes[pos];
    nodes = snode_nodes(snode);
    orients = snode_orients(snode);

    if(path->superorients[pos] == FORWARD)
    {
      // Supernode is forwards
      len = MIN2(snode->num_of_nodes, maxlen-j);
      memcpy(nlist+j, nodes, len*sizeof(hkey_t));
      memcpy(olist+j, orients, len*sizeof(Orientation));
      j += len;
    }
    else
    {
      // Supernode is reverse
      for(i = snode->num_of_nodes - 1; j <= maxlen-1; i--) {
        nlist[j] = nodes[i];
        olist[j] = opposite_orientation(orients[i]);
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
                         size_t flank5pkmers, size_t threadid,
                         size_t max_allele_len, size_t max_flank_len)
{
  // tmp variables
  size_t max = MAX2(max_allele_len, max_flank_len);
  hkey_t tmp_e[max];
  Orientation tmp_o[max];
  size_t i, num_kmers;

  // 5p flank
  gzprintf(out, ">var_%zu.%zu_5p_flank length=%zu\n", threadid, bnum, flank5pkmers);
  print_branch(flank5pe, flank5po, flank5pkmers, true, db_graph, out);

  // 3p flank
  num_kmers = suppathpos_to_list(spp_arr[0], spp_arr[0]->pos, 1,
                                 tmp_e, tmp_o, max_flank_len);
  gzprintf(out, ">var_%zu.%zu_3p_flank length=%zu\n", threadid, bnum, num_kmers);
  print_branch(tmp_e, tmp_o, num_kmers, false, db_graph, out);

  // Print alleles
  for(i = 0; i < num_of_paths; i++)
  {
    num_kmers = suppathpos_to_list(spp_arr[i], 0, spp_arr[i]->pos,
                                   tmp_e, tmp_o, max_allele_len);
    gzprintf(out, ">var_%zu.%zu_branch_%zu length=%zu\n", threadid, bnum, i, num_kmers);
    print_branch(tmp_e, tmp_o, num_kmers, false, db_graph, out);
  }

  gzputc(out, '\n');
}

// Returns 1 if successfully walked across supernode, 0 otherwise
static inline void walk_supernode_end(GraphWalker *wlk, CallerSupernode *snode,
                                      SuperOrientation snorient,
                                      Nucleotide *lost_nuc)
{
  // Only need to traverse the first and last nodes of a supernode
  size_t last = snode->num_of_nodes-1, kmer_size = wlk->db_graph->kmer_size;
  hkey_t last_node;
  Orientation last_orient;
  BinaryKmer last_bkmer;

  if(last > 0) {
    if(snorient == FORWARD) {
      last_node = snode_nodes(snode)[last];
      last_orient = snode_orients(snode)[last];
    } else {
      last_node = snode_nodes(snode)[0];
      last_orient = rev_orient(snode_orients(snode)[0]);
    }

    last_bkmer = db_graph_oriented_bkmer(wlk->db_graph, last_node, last_orient);
    graph_traverse_force_jump(wlk, last_node, last_bkmer, false);
    // don't need counter paths here (we're at the end of a supernode)
    *lost_nuc = binary_kmer_first_nuc(last_bkmer, kmer_size);
  }
  else {
    // XOR snorient => negate/reverse if snorient is REVERSE
    // 0 0 1 1
    // 0 1 0 1
    last_bkmer = db_node_bkmer(wlk->db_graph, snode_nodes(snode)[0]);
    *lost_nuc = db_node_first_nuc(last_bkmer, snode_orients(snode)[0] ^ snorient,
                                  kmer_size);
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
                             CallerNodeBuf *nbuf,
                             CallerSupernode *snode_store,
                             SupernodePathPos *snodepos_store,
                             size_t *snode_count_ptr,
                             size_t *snodepos_count_ptr,
                             size_t max_allele_len)
{
  CallerSupernode *snode;
  hkey_t node2, *nodes; // node,node2 are start/end nodes of supernode
  int hashret;
  khiter_t k;
  boolean supernode_already_exists;
  SuperOrientation snorient;

  const dBGraph *db_graph = wlk->db_graph;
  const size_t kmer_size = db_graph->kmer_size;

  #ifdef DEBUG_CALLER
    char tmp[MAX_KMER_SIZE+1];
    printf("new allele\n");
  #endif

  size_t snode_count = *snode_count_ptr;
  size_t snodepos_count = *snodepos_count_ptr;

  size_t supindx, kmers_in_path = 0;

  for(supindx = 0; ; supindx++)
  {
    #ifdef DEBUG_CALLER
      binary_kmer_to_str(db_node_bkmer(db_graph,node), kmer_size, tmp);
      printf(" load_allele_path: %s:%i\n", tmp, or);
    #endif

    // Find or add supernode
    k = kh_put(supnode_hsh, snode_hash, (uint64_t)node, &hashret);
    supernode_already_exists = (hashret == 0);

    if(supernode_already_exists)
    {
      // already in supernode hash table
      snode = kh_value(snode_hash, k);

      #ifdef DEBUG_CALLER
        printf("  (supernode already exists)\n");
      #endif
    }
    else
    {
      // add to supernode hash table
      snode = &snode_store[snode_count++];
      snode->nbuf = nbuf;
      caller_supernode_create(node, or, snode, db_graph);
      nodes = snode_nodes(snode);

      kh_value(snode_hash, k) = snode;

      // Add end node to hash
      if(snode->num_of_nodes > 1)
      {
        node2 = nodes[node == nodes[0] ? snode->num_of_nodes-1 : 0];
        k = kh_put(supnode_hsh, snode_hash, (uint64_t)node2, &hashret);
        kh_value(snode_hash, k) = snode;
      }
    }

    snorient = supernode_get_orientation(snode, node, or);

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

    // Done adding new supernode
    // Check path length
    kmers_in_path += snode->num_of_nodes;
    if(kmers_in_path > max_allele_len) break;

    // Check if we've already traversed this supernode
    if(db_node_has_traversed(visited, node, or)) break;
    db_node_set_traversed(visited, node, or);

    // Find next node
    Nucleotide lost_nuc;
    walk_supernode_end(wlk, snode, snorient, &lost_nuc);

    size_t i, num_edges;
    hkey_t *next_nodes;
    Orientation *next_orients;
    BinaryKmer next_bkmers[4], next_bkmer;
    Nucleotide next_bases[4];
    int nxt_idx;
    boolean isfork;

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
    for(i = 0; i < num_edges; i++)
    {
      next_bkmers[i] = db_node_bkmer(db_graph, next_nodes[i]);
      next_bases[i] = db_node_last_nuc(next_bkmers[i], next_orients[i], kmer_size);
    }

    nxt_idx = graph_walker_choose(wlk, num_edges, next_nodes, next_bases, &isfork);
    if(nxt_idx == -1) break;

    node = next_nodes[nxt_idx];
    or = next_orients[nxt_idx];
    next_bkmer = db_node_oriented_bkmer(next_bkmers[nxt_idx], or, kmer_size);

    graph_traverse_force_jump(wlk, node, next_bkmer, isfork);
    graph_walker_node_add_counter_paths(wlk, lost_nuc);
  }
  // printf("DONE\n");

  *snode_count_ptr = snode_count;
  *snodepos_count_ptr = snodepos_count;
}

// spphash is a temporary hash table used to through out duplicate supernode
// path positions
// We can have duplicate paths if two or more colours find the same path through
// the graph during bubble discovery
static size_t remove_snpath_pos_dupes(SupernodePathPos **results, size_t arrlen,
                                      khash_t(snpps_hsh) *spphash)
{
  size_t i, j;
  int hashret;
  kh_clear(snpps_hsh, spphash);

  for(i = 0, j = 0; i < arrlen; i++) {
    kh_put(snpps_hsh, spphash, results[i], &hashret);
    if(hashret != 0) results[j++] = results[i];
  }
  return j;
}

static void get_prev_supernode_pos(const SupernodePathPos *spp,
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

static char is_bubble_flank(SupernodePathPos *const* spp_arr, size_t num)
{
  if(num == 0) return 0;

  size_t i;
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

// Returns 0 or 1
static boolean path_in_colour(const SupernodePathPos *pp, size_t col,
                              const dBGraph *db_graph)
{
  CallerNodeBuf *nbuf;
  size_t i, j;
  for(i = 0; i <= pp->pos; i++) {
    nbuf = pp->path->supernodes[i]->nbuf;
    for(j = 0; j < nbuf->len; j++)
      if(db_node_has_col(db_graph, nbuf->nodes[j], col))
        return false;
  }
  return true;
}

// Returns number of paths
static size_t remove_ref_paths(SupernodePathPos **spp_arr, size_t num_paths,
                               const size_t *ref_cols, size_t num_ref,
                               const dBGraph *db_graph)
{
  size_t r, p, seen[num_ref];
  memset(seen, 0, sizeof(size_t)*num_ref);

  for(p = 0; p < num_paths; )
  {
    for(r = 0; r < num_ref; r++)
    {
      if(path_in_colour(spp_arr[p], ref_cols[r], db_graph))
      {
        // Drop path if already seen
        if(seen[r]) break;
        seen[r] = 1;
      }
    }
    if(r < num_ref)
    {
      // Drop path
      SupernodePathPos *tmp_spp;
      SWAP(spp_arr[p], spp_arr[num_paths-1], tmp_spp);
      num_paths--;
    }
    else p++;
  }

  return num_paths;
}

// Potential bubble - filter ref and duplicate alleles
static void resolve_bubble(SupernodePathPos **snodepathposes, size_t num,
                           gzFile out, size_t *bnum, size_t threadid,
                           size_t max_allele_len, size_t max_flank_len,
                           hkey_t *flank5pe, Orientation *flank5po,
                           size_t *flank5pkmers, const dBGraph *db_graph,
                           const size_t *ref_cols, size_t num_ref,
                           khash_t(snpps_hsh) *spp_hash)
{
  if(num < 2 || !is_bubble_flank(snodepathposes, num)) return;
  num = remove_ref_paths(snodepathposes, num, ref_cols, num_ref, db_graph);
  if(num < 2) return;
  num = remove_snpath_pos_dupes(snodepathposes, num, spp_hash);
  if(num < 2) return;

  if(*flank5pkmers == 0)
  {
    // Haven't fetched 5p flank yet
    boolean iscycle;
    int len = supernode_extend(db_graph, &flank5pe, &flank5po, 0,
                               &iscycle, &max_flank_len, false);
    *flank5pkmers = (len == -1 ? max_flank_len : (unsigned)len);
    supernode_reverse(flank5pe, flank5po, *flank5pkmers);
  }

  print_bubble(out, *bnum, db_graph,
               snodepathposes, num,
               flank5pe, flank5po, *flank5pkmers, threadid,
               max_allele_len, max_flank_len);
  (*bnum)++;
}

static void find_bubbles(hkey_t fork_n, Orientation fork_o,
                         const dBGraph *db_graph, GraphWalker *wlk,
                         uint64_t *visited,
                         khash_t(supnode_hsh) *snode_hash,
                         khash_t(snpps_hsh) *spp_hash,
                         CallerNodeBuf *nbuf,
                         SupernodePath *paths,
                         CallerSupernode *snode_store,
                         SupernodePathPos *snodepos_store,
                         gzFile out, size_t *bnum, size_t threadid,
                         size_t max_allele_len, size_t max_flank_len,
                         const size_t *ref_cols, size_t num_ref)
{
  // Set number of supernodes, positions, nodes to zero
  size_t snode_count = 0, snodepos_count = 0;
  nbuf->len = 0;

  // Clear the hash table of supernodes
  kh_clear(supnode_hsh, snode_hash);

  hkey_t nodes[4];
  Orientation orients[4];
  Nucleotide bases[4];
  size_t i, num_next;

  num_next = db_graph_next_nodes(db_graph, db_node_bkmer(db_graph, fork_n),
                                 fork_o, db_node_edges(db_graph, 0, fork_n),
                                 nodes, orients, bases);

  #ifdef DEBUG_CALLER
    char tmpstr[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_bkmer(db_graph, fork_n), db_graph->kmer_size, tmpstr);
    printf("fork %s:%i out-degree:%i\n", tmpstr, (int)fork_o, (int)num_next);
  #endif

  // loop over alleles, then colours
  size_t supindx, num_of_paths = 0;
  Colour colour, colours_loaded = db_graph->num_of_cols_used;
  Nucleotide lost_nuc;
  SupernodePath *path;
  CallerSupernode *snode;
  int col_has_node[4], num_edges_in_col;

  for(colour = 0; colour < colours_loaded; colour++)
  {
    num_edges_in_col = 0;

    for(i = 0; i < num_next; i++) {
      col_has_node[i] = (db_node_has_col(db_graph, nodes[i], colour) > 0);
      num_edges_in_col += col_has_node[i];
    }

    for(i = 0; i < num_next; i++)
    {
      if(col_has_node[i])
      {
        path = paths + num_of_paths++;

        graph_walker_init(wlk, db_graph, colour, colour, fork_n, fork_o);
        lost_nuc = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

        graph_traverse_force(wlk, nodes[i], bases[i], num_edges_in_col > 1);
        graph_walker_node_add_counter_paths(wlk, lost_nuc);

        // Constructs a path of supernodes (SupernodePath)
        load_allele_path(nodes[i], orients[i], path, snode_hash, wlk, visited,
                         nbuf, snode_store, snodepos_store,
                         &snode_count, &snodepos_count, max_allele_len);

        graph_walker_finish(wlk);

        // Remove mark traversed
        for(supindx = 0; supindx < snode_count; supindx++)
        {
          snode = snode_store + supindx;
          size_t last = snode->num_of_nodes-1;
          db_node_fast_clear_traversed(visited, snode_nodes(snode)[0]);
          db_node_fast_clear_traversed(visited, snode_nodes(snode)[last]);
        }
      }
    }
  }

  hkey_t flank5_nstore[max_flank_len], *flank5pe = flank5_nstore;
  Orientation flank5_ostore[max_flank_len], *flank5po = flank5_ostore;
  size_t flank5pkmers = 0;
  flank5pe[0] = fork_n;
  flank5po[0] = opposite_orientation(fork_o);

  // Loop over supernodes checking if they are 3p flanks
  for(i = 0; i < snode_count; i++)
  {
    #ifdef DEBUG_CALLER
      char tmpsup[MAX_KMER_SIZE+1];
      BinaryKmer bkmer = db_node_bkmer(db_graph, snode_nodes(&snode_store[i])[0]);
      binary_kmer_to_str(bkmer, db_graph->kmer_size, tmpsup);
      printf("check supernode: %s\n", tmpsup);
    #endif

    SupernodePathPos *pp = snode_store[i].first_pathpos;
    // pp may be null if the node could not be traversed by any path
    // e.g. supernode is also 5p flank and no assembly information available
    // DEV: is this still true?
    if(pp != NULL && pp->next != NULL)
    {
      // possible 3p flank (i.e. bubble end)
      // there a 4 possible allele paths per colour
      // each path can hit a node at most NUM_OF_SHADES times
      size_t maxnodes = max_allele_len * db_graph->num_of_cols * 4;
      SupernodePathPos *spp_forward[maxnodes];
      SupernodePathPos *spp_reverse[maxnodes];
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

      // #ifdef DEBUG_CALLER
      // gzprintf(out, " num_forward: %zu; num_reverse: %zu\n", num_forward, num_reverse);
      // #endif

      resolve_bubble(spp_forward, num_forward, out, bnum, threadid,
                     max_allele_len, max_flank_len, flank5pe, flank5po,
                     &flank5pkmers, db_graph, ref_cols, num_ref, spp_hash);

      resolve_bubble(spp_reverse, num_reverse, out, bnum, threadid,
                     max_allele_len, max_flank_len, flank5pe, flank5po,
                     &flank5pkmers, db_graph, ref_cols, num_ref, spp_hash);

      /*
      if(!is_bubble_flank(spp_forward, num_forward)) num_forward = 0;
      if(!is_bubble_flank(spp_reverse, num_reverse)) num_reverse = 0;

      // Remove paths that are in ref if we have >1 of such paths
      if(num_forward > 1) {
        num_forward = remove_ref_paths(spp_forward, num_forward, ref_cols, num_ref, db_graph);
        if(num_forward > 1)
          num_forward = remove_snpath_pos_dupes(spp_forward, num_forward, spp_hash);
      }

      if(num_reverse > 1) {
        num_reverse = remove_ref_paths(spp_reverse, num_reverse, ref_cols, num_ref, db_graph);
        if(num_reverse > 1)
          num_reverse = remove_snpath_pos_dupes(spp_reverse, num_reverse, spp_hash);
      }

      if(num_forward > 1 || num_reverse > 1)
      {
        if(flank5pkmers == 0)
        {
          // Haven't fetched 5p flank yet
          boolean iscycle;
          int len = supernode_extend(db_graph, &flank5pe, &flank5po, 0,
                                     &iscycle, &max_flank_len, false);
          flank5pkmers = (len == -1 ? max_flank_len : (unsigned)len);
          supernode_reverse(flank5pe, flank5po, flank5pkmers);
        }

        // Check if actual bubbles in fw / rv
        if(num_forward > 1)
        {
          print_bubble(out, *bnum, db_graph,
                       spp_forward, num_forward,
                       flank5pe, flank5po, flank5pkmers, threadid,
                       max_allele_len, max_flank_len);
          (*bnum)++;
        }

        if(num_reverse > 1)
        {
          print_bubble(out, *bnum, db_graph,
                       spp_reverse, num_reverse,
                       flank5pe, flank5po, flank5pkmers, threadid,
                       max_allele_len, max_flank_len);
          (*bnum)++;
        }
      }
      */
    }
  }
}


struct caller_region_t
{
  const dBGraph *db_graph;
  const gzFile out;
  const uint64_t start_hkey, end_hkey;
  const size_t threadid;
  const size_t max_allele_len, max_flank_len, *ref_cols, num_ref;
  size_t num_of_bubbles;
};

void* bubble_caller(void *args)
{
  struct caller_region_t *tdata = (struct caller_region_t *)args;
  const dBGraph *db_graph = tdata->db_graph;
  const size_t max_allele_len = tdata->max_allele_len;
  const size_t max_flank_len = tdata->max_flank_len;
  gzFile out = tdata->out;

  // Arrays to re-use
  size_t i, npaths = 4 * db_graph->num_of_cols;

  size_t node_bits = round_bits_to_words64(db_graph->ht.capacity);
  uint64_t *visited = calloc2(2*node_bits, sizeof(uint64_t));

  size_t maxnodes = max_allele_len * db_graph->num_of_cols * 4;

  SupernodePath snode_paths[npaths];
  CallerSupernode **supernodes = malloc2(maxnodes * sizeof(CallerSupernode*));
  SuperOrientation *superorients = malloc2(maxnodes * sizeof(SuperOrientation));

  for(i = 0; i < npaths; i++) {
    snode_paths[i].supernodes = supernodes + i * max_allele_len;
    snode_paths[i].superorients = superorients + i * max_allele_len;
  }

  CallerSupernode *snode_store = malloc2(maxnodes * sizeof(CallerSupernode));
  SupernodePathPos *snodepos_store = malloc2(maxnodes * sizeof(SupernodePathPos));

  khash_t(supnode_hsh) *snode_hash = kh_init(supnode_hsh);
  khash_t(snpps_hsh) *spp_hash = kh_init(snpps_hsh);

  // Max usage is 4 * max_allele_len * cols
  // size_t nodeslen = 4 * max_allele_len * 2;
  hkey_t *nodes = malloc2(maxnodes * sizeof(hkey_t));
  Orientation *orients = malloc2(maxnodes * sizeof(Orientation));

  CallerNodeBuf nbuf = {.nodes = nodes, .orients = orients,
                        .len = 0, .cap = maxnodes};

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  // BinaryKmer tmpkmer;
  // hkey_t node;
  // tmpkmer = binary_kmer_from_str("AAGGTGACTGATCTGCCTGAGAAATGAACCT", db_graph->kmer_size);
  // node = hash_table_find(&db_graph->ht, tmpkmer);

  BinaryKmer *table = db_graph->ht.table;
  BinaryKmer *ptr = table + tdata->start_hkey;
  BinaryKmer *end = table + tdata->end_hkey;

  for(; ptr < end; ptr++) {
    if(HASH_ENTRY_ASSIGNED(*ptr)) {
      hkey_t node = ptr - table;
      Edges edges = db_node_edges(db_graph, 0, node);
      if(edges_get_outdegree(edges, FORWARD) > 1) {
        find_bubbles(node, FORWARD, db_graph, &wlk, visited,
                     snode_hash, spp_hash, &nbuf,
                     snode_paths, snode_store, snodepos_store,
                     out, &tdata->num_of_bubbles, tdata->threadid,
                     max_allele_len, max_flank_len,
                     tdata->ref_cols, tdata->num_ref);
      }
      if(edges_get_outdegree(edges, REVERSE) > 1) {
        find_bubbles(node, REVERSE, db_graph, &wlk, visited,
                     snode_hash, spp_hash, &nbuf,
                     snode_paths, snode_store, snodepos_store,
                     out, &tdata->num_of_bubbles, tdata->threadid,
                     max_allele_len, max_flank_len,
                     tdata->ref_cols, tdata->num_ref);
      }
    }
  }


  graph_walker_dealloc(&wlk);

  free(visited);
  free(snode_store);
  free(snodepos_store);
  free(supernodes);
  free(superorients);
  kh_destroy(supnode_hsh, snode_hash);
  kh_destroy(snpps_hsh, spp_hash);
  free(nbuf.nodes);
  free(nbuf.orients);

  return NULL;
}

// max_allele_len, max_flank_len in kmers
void invoke_bubble_caller(const dBGraph *db_graph, const char* out_file,
                          int num_threads, char **tmp_paths,
                          size_t max_allele_len, size_t max_flank_len,
                          const size_t *ref_cols, size_t num_ref,
                          const CmdArgs *args)
{
  assert(db_graph->num_edge_cols == 1);

  // Open output file
  gzFile out = gzopen(out_file, "w");

  if(out == NULL)
    die("Cannot open paths bubble caller output file: %s", out_file);

  // Print header
  print_calling_header(db_graph, out, out_file, args);

  size_t num_of_bubbles = 0;

  if(num_threads <= 1)
  {
    struct caller_region_t tdata
      = {.db_graph = db_graph, .out = out,
         .start_hkey = 0, .end_hkey = db_graph->ht.capacity, .threadid = 0,
         .max_allele_len = max_allele_len, .max_flank_len = max_flank_len,
         .ref_cols = ref_cols, .num_ref = num_ref, .num_of_bubbles = 0};
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

      struct caller_region_t tmptdata
        = {.db_graph = db_graph, .out = tmpout,
           .start_hkey = start, .end_hkey = end, .threadid = i,
           .max_allele_len = max_allele_len, .max_flank_len = max_flank_len,
           .ref_cols = ref_cols, .num_ref = num_ref, .num_of_bubbles = 0};

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
  status("%s bubbles called with Paths-Bubble-Caller\n", num_bubbles_str);
  status("  saved to: %s\n", out_file);

  gzclose(out);
}
