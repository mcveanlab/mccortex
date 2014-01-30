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
#include "cmd.h"
#include "seq_reader.h"
#include "path_store.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "caller_supernode.h"
#include "supernode.h"

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

void bubble_caller_print_header(const dBGraph *db_graph, gzFile out,
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

static void print_branch(dBNode *nodes, size_t len, boolean print_first_kmer,
                         const dBGraph *db_graph, gzFile out)
{
  size_t i = print_first_kmer, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;

  if(print_first_kmer) {
    char tmp[MAX_KMER_SIZE+1];
    bkmer = db_graph_oriented_bkmer(db_graph, nodes[0].key, nodes[0].orient);
    binary_kmer_to_str(bkmer, kmer_size, tmp);
    gzputs(out, tmp);
  }

  // i = 1 if print_first_kmer, otherwise 0
  for(; i < len; i++) {
    bkmer = db_node_get_bkmer(db_graph, nodes[i].key);
    nuc = db_node_last_nuc(bkmer, nodes[i].orient, kmer_size);
    gzputc(out, dna_nuc_to_char(nuc));
  }

  gzputc(out, '\n');
}

// Returns number of nodes
// if maxlen is zero it is ignored
static size_t suppathpos_to_list(const SupernodePathPos *snodepathpos,
                                 size_t start, size_t suplen,
                                 dBNode *nlist, size_t maxlen)
{
  SupernodePath *path = snodepathpos->path;
  CallerSupernode *snode;
  dBNode *nodes;

  size_t pos, i, j = 0, end = start+suplen, len;

  for(pos = start; pos < end && j < maxlen; pos++)
  {
    snode = path->supernodes[pos];
    nodes = snode_nodes(snode);

    if(path->superorients[pos] == FORWARD)
    {
      // Supernode is forwards
      len = MIN2(snode->num_of_nodes, maxlen-j);
      memcpy(nlist+j, nodes, len*sizeof(dBNode));
      j += len;
    }
    else
    {
      // Supernode is reverse
      for(i = snode->num_of_nodes - 1; j <= maxlen-1; i--) {
        nlist[j].key = nodes[i].key;
        nlist[j].orient = opposite_orientation(nodes[i].orient);
        j++;
        if(i == 0) break;
      }
    }
  }
  return j;
}

static void print_bubble(gzFile out, size_t bnum,
                         SupernodePathPos **spp_arr, size_t num_of_paths,
                         dBNode *flank5p, size_t flank5pkmers,
                         size_t max_allele_len, size_t max_flank_len,
                         size_t threadid, const dBGraph *db_graph)
{
  // tmp variables
  size_t max = MAX2(max_allele_len, max_flank_len);
  dBNode tmp_nodes[max];
  size_t i, num_kmers;

  // 5p flank
  gzprintf(out, ">var_%zu.%zu_5p_flank length=%zu\n", threadid, bnum, flank5pkmers);
  print_branch(flank5p, flank5pkmers, true, db_graph, out);

  // 3p flank
  num_kmers = suppathpos_to_list(spp_arr[0], spp_arr[0]->pos, 1,
                                 tmp_nodes, max_flank_len);
  gzprintf(out, ">var_%zu.%zu_3p_flank length=%zu\n", threadid, bnum, num_kmers);
  print_branch(tmp_nodes, num_kmers, false, db_graph, out);

  // Print alleles
  for(i = 0; i < num_of_paths; i++)
  {
    num_kmers = suppathpos_to_list(spp_arr[i], 0, spp_arr[i]->pos,
                                   tmp_nodes, max_allele_len);
    gzprintf(out, ">var_%zu.%zu_branch_%zu length=%zu\n", threadid, bnum, i, num_kmers);
    print_branch(tmp_nodes, num_kmers, false, db_graph, out);
  }

  gzputc(out, '\n');
}

// Returns 1 if successfully walked across supernode, 0 otherwise
static inline void walk_supernode_end(GraphWalker *wlk, CallerSupernode *snode,
                                      SuperOrientation snorient)
{
  // Only need to traverse the first and last nodes of a supernode
  const size_t last = snode->num_of_nodes-1;
  dBNode lastnode;
  BinaryKmer last_bkmer;

  if(last > 0) {
    if(snorient == FORWARD) {
      lastnode = snode_nodes(snode)[last];
    } else {
      lastnode = db_node_reverse(snode_nodes(snode)[0]);
    }
    last_bkmer = db_graph_oriented_bkmer(wlk->db_graph, lastnode.key, lastnode.orient);
    graph_walker_jump_snode_end(wlk, lastnode.key, last_bkmer);
  }
}

// Constructs a path of supernodes (SupernodePath)
// returns number of supernodes loaded
static void load_allele_path(dBNode node,
                             SupernodePath *path,
                             khash_t(supnode_hsh) *snode_hash,
                             GraphWalker *wlk, // walker set to go
                             RepeatWalker *rptwlk, // cleared and ready
                             // these 4 params are tmp memory
                             dBNodeBuffer *nbuf,
                             CallerSupernode *snode_store,
                             SupernodePathPos *snodepos_store,
                             size_t *snode_count_ptr,
                             size_t *snodepos_count_ptr,
                             size_t max_allele_len)
{
  CallerSupernode *snode;
  dBNode node2, *nodes; // node,node2 are start/end nodes of supernode
  int hashret;
  khiter_t k;
  boolean supernode_already_exists;
  SuperOrientation snorient;

  const dBGraph *db_graph = wlk->db_graph;

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
      binary_kmer_to_str(db_node_get_bkmer(db_graph,node.key), kmer_size, tmp);
      printf(" load_allele_path: %s:%i\n", tmp, node.orient);
    #endif

    // Find or add supernode beginning with given node
    k = kh_put(supnode_hsh, snode_hash, (uint64_t)node.key, &hashret);
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
      caller_supernode_create(node, snode, db_graph);
      nodes = snode_nodes(snode);

      kh_value(snode_hash, k) = snode;

      // Add end node to hash
      if(snode->num_of_nodes > 1)
      {
        node2 = nodes[node.key == nodes[0].key ? snode->num_of_nodes-1 : 0];
        k = kh_put(supnode_hsh, snode_hash, (uint64_t)node2.key, &hashret);
        kh_value(snode_hash, k) = snode;
      }
    }

    snorient = supernode_get_orientation(snode, node);

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

    // Traverse to the end of the supernode
    walk_supernode_end(wlk, snode, snorient);

    // Find next node
    uint8_t num_edges;
    const dBNode *next_nodes;
    const Nucleotide *next_bases;

    if(snorient == FORWARD) {
      num_edges = snode->num_next;
      next_nodes = snode->next_nodes;
      next_bases = snode->next_bases;
    }
    else {
      num_edges = snode->num_prev;
      next_nodes = snode->prev_nodes;
      next_bases = snode->prev_bases;
    }

    if(!graph_traverse_nodes(wlk, num_edges, next_nodes, next_bases) ||
       !rpt_walker_attempt_traverse(rptwlk, wlk)) break;

    node = wlk->node;
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
  dBNodeBuffer *nbuf;
  size_t i, j;
  for(i = 0; i <= pp->pos; i++) {
    nbuf = pp->path->supernodes[i]->nbuf;
    for(j = 0; j < nbuf->len; j++)
      if(db_node_has_col(db_graph, nbuf->data[j].key, col))
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
                           dBNodeBuffer *flank5p,
                           const dBGraph *db_graph,
                           const size_t *ref_cols, size_t num_ref,
                           khash_t(snpps_hsh) *spp_hash)
{
  if(num < 2 || !is_bubble_flank(snodepathposes, num)) return;
  num = remove_ref_paths(snodepathposes, num, ref_cols, num_ref, db_graph);
  if(num < 2) return;
  num = remove_snpath_pos_dupes(snodepathposes, num, spp_hash);
  if(num < 2) return;

  if(flank5p->len == 0)
  {
    // Haven't fetched 5p flank yet
    // flank5p[0] already contains the first node
    flank5p->len = 1;
    supernode_extend(flank5p, max_flank_len, db_graph);
    supernode_reverse(flank5p->data, flank5p->len);
  }

  print_bubble(out, *bnum,
               snodepathposes, num,
               flank5p->data, flank5p->len,
               max_allele_len, max_flank_len,
               threadid, db_graph);
  (*bnum)++;
}

static void find_bubbles(hkey_t fork_n, Orientation fork_o,
                         const dBGraph *db_graph,
                         GraphWalker *wlk, RepeatWalker *rptwlk,
                         khash_t(supnode_hsh) *snode_hash,
                         khash_t(snpps_hsh) *spp_hash,
                         dBNodeBuffer *nbuf,
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

  dBNode nodes[4];
  Nucleotide bases[4];
  size_t i, num_next;

  num_next = db_graph_next_nodes(db_graph, db_node_get_bkmer(db_graph, fork_n),
                                 fork_o, db_node_edges(db_graph, 0, fork_n),
                                 nodes, bases);

  #ifdef DEBUG_CALLER
    char tmpstr[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_get_bkmer(db_graph, fork_n), db_graph->kmer_size, tmpstr);
    printf("fork %s:%i out-degree:%i\n", tmpstr, (int)fork_o, (int)num_next);
  #endif

  // loop over alleles, then colours
  size_t supindx, num_of_paths = 0;
  Colour colour, colours_loaded = db_graph->num_of_cols_used;
  SupernodePath *path;
  CallerSupernode *snode;
  int node_has_col[4], num_edges_in_col;

  for(colour = 0; colour < colours_loaded; colour++)
  {
    if(!db_node_has_col(db_graph, fork_n, colour)) continue;

    // Determine if this fork is a fork in the current colour
    num_edges_in_col = 0;

    for(i = 0; i < num_next; i++) {
      node_has_col[i] = (db_node_has_col(db_graph, nodes[i].key, colour) > 0);
      num_edges_in_col += node_has_col[i];
    }

    for(i = 0; i < num_next; i++)
    {
      if(node_has_col[i])
      {
        path = paths + num_of_paths;
        num_of_paths++;

        dBNode node = {.key = fork_n, .orient = fork_o};
        graph_walker_init(wlk, db_graph, colour, colour, node);

        graph_traverse_force(wlk, nodes[i].key, bases[i], num_edges_in_col > 1);

        // Constructs a path of supernodes (SupernodePath)
        load_allele_path(nodes[i], path, snode_hash, wlk, rptwlk,
                         nbuf, snode_store, snodepos_store,
                         &snode_count, &snodepos_count, max_allele_len);

        graph_walker_finish(wlk);

        // Remove mark traversed
        rpt_walker_fast_clear(rptwlk, NULL, 0);

        for(supindx = 0; supindx < snode_count; supindx++)
        {
          snode = snode_store + supindx;
          size_t last = snode->num_of_nodes-1;
          rpt_walker_fast_clear2(rptwlk, snode_nodes(snode)[0]);
          rpt_walker_fast_clear2(rptwlk, snode_nodes(snode)[last]);
        }
      }
    }
  }

  dBNode flank5_store[max_flank_len];
  dBNodeBuffer flank5p = {.data = flank5_store, .len = 0,
                          .capacity = max_flank_len};

  flank5p.data[0].key = fork_n;
  flank5p.data[0].orient = opposite_orientation(fork_o);

  // Loop over supernodes checking if they are 3p flanks
  for(i = 0; i < snode_count; i++)
  {
    #ifdef DEBUG_CALLER
      char tmpsup[MAX_KMER_SIZE+1];
      BinaryKmer bkmer = db_node_get_bkmer(db_graph, snode_nodes(&snode_store[i])[0].key);
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
      // there are 4 possible allele paths per colour
      // each path can hit a node at most max_allele_len times
      // DEV: move this to the heap
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
      //   gzprintf(out, " num_forward: %zu; num_reverse: %zu\n",
      //            num_forward, num_reverse);
      // #endif

      resolve_bubble(spp_forward, num_forward, out, bnum, threadid,
                     max_allele_len, max_flank_len, &flank5p,
                     db_graph, ref_cols, num_ref, spp_hash);

      resolve_bubble(spp_reverse, num_reverse, out, bnum, threadid,
                     max_allele_len, max_flank_len, &flank5p,
                     db_graph, ref_cols, num_ref, spp_hash);
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

  RepeatWalker rptwlk;
  rpt_walker_alloc(&rptwlk, db_graph->ht.capacity, 22); // 4MB

  // Max usage is 4 * max_allele_len * cols
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

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, maxnodes);

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
      hkey_t node = (hkey_t)(ptr - table);
      Edges edges = db_node_get_edges(db_graph, 0, node);
      if(edges_get_outdegree(edges, FORWARD) > 1) {
        find_bubbles(node, FORWARD, db_graph, &wlk, &rptwlk,
                     snode_hash, spp_hash, &nbuf,
                     snode_paths, snode_store, snodepos_store,
                     out, &tdata->num_of_bubbles, tdata->threadid,
                     max_allele_len, max_flank_len,
                     tdata->ref_cols, tdata->num_ref);
      }
      if(edges_get_outdegree(edges, REVERSE) > 1) {
        find_bubbles(node, REVERSE, db_graph, &wlk, &rptwlk,
                     snode_hash, spp_hash, &nbuf,
                     snode_paths, snode_store, snodepos_store,
                     out, &tdata->num_of_bubbles, tdata->threadid,
                     max_allele_len, max_flank_len,
                     tdata->ref_cols, tdata->num_ref);
      }
    }
  }

  rpt_walker_dealloc(&rptwlk);
  graph_walker_dealloc(&wlk);

  free(snode_store);
  free(snodepos_store);
  free(supernodes);
  free(superorients);
  kh_destroy(supnode_hsh, snode_hash);
  kh_destroy(snpps_hsh, spp_hash);
  db_node_buf_dealloc(&nbuf);

  return NULL;
}

// max_allele_len, max_flank_len in kmers
void invoke_bubble_caller(const dBGraph *db_graph, gzFile gzout,
                          const size_t num_threads, char **tmp_paths,
                          size_t max_allele_len, size_t max_flank_len,
                          const size_t *ref_cols, size_t num_ref)
{
  assert(db_graph->num_edge_cols == 1);

  size_t num_of_bubbles = 0;

  if(num_threads <= 1)
  {
    struct caller_region_t tdata
      = {.db_graph = db_graph, .out = gzout,
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

    const uint64_t capacity = db_graph->ht.capacity;
    uint64_t start = 0, end;
    gzFile tmpout;
    size_t i;

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
      while((len = gzread(in, data, MEGABYTE)) > 0)
        gzwrite(gzout, data, (unsigned int)len);
      gzclose(in);

      num_of_bubbles += tdata[i].num_of_bubbles;
    }
  }

  char num_bubbles_str[100];
  ulong_to_str(num_of_bubbles, num_bubbles_str);
  status("%s bubbles called with Paths-Bubble-Caller\n", num_bubbles_str);
}
