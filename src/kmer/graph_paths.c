#include "global.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_paths.h"
#include "packed_path.h"
#include "binary_seq.h"
#include "path_format.h"
#include "sorted_path_set.h"

#include "bit_array/bit_macros.h"

// Similar to path_file_filter.c:path_file_load_check()
// Check kmer size matches and sample names match
void graphs_paths_compatible(const GraphFileReader *graphs, size_t num_graphs,
                             const PathFileReader *paths, size_t num_paths)
{
  size_t g, p, kmer_size;
  size_t ctx_max_cols = 0, ctp_max_cols = 0;
  size_t ctx_max_kmers = 0, ctp_max_kmers = 0;

  if(num_graphs) kmer_size = graphs[0].hdr.kmer_size;
  else if(num_paths) kmer_size = paths[0].hdr.kmer_size;
  else return;

  for(g = 0; g < num_graphs; g++) {
    if(graphs[g].hdr.kmer_size != kmer_size) {
      die("Kmer-size doesn't match between files [%zu vs %u]: %s",
          kmer_size, graphs[g].hdr.kmer_size, graphs[g].fltr.orig_path.buff);
    }
    ctx_max_cols = MAX2(ctx_max_cols, graph_file_usedcols(&graphs[g]));
    ctx_max_kmers = MAX2(ctx_max_kmers, graphs[g].num_of_kmers);
  }

  for(p = 0; p < num_paths; p++) {
    if(paths[p].hdr.kmer_size != kmer_size) {
      die("Kmer-size doesn't match between files [%zu vs %u]: %s",
          kmer_size, paths[p].hdr.kmer_size, paths[p].fltr.orig_path.buff);
    }
    ctp_max_cols = MAX2(ctp_max_cols, path_file_usedcols(&paths[p]));
    ctp_max_kmers = MAX2(ctp_max_kmers, paths[p].hdr.num_kmers_with_paths);
  }

  const FileFilter *fltr = (num_graphs ? &graphs[0].fltr : &paths[0].fltr);
  db_graph_check_kmer_size(kmer_size, fltr->orig_path.buff);

  // ctx_max_kmers may be zero if reading from a stream
  if(ctx_max_kmers > 0 && ctp_max_kmers > ctx_max_kmers)
    die("More kmers in path files than in graph files!");

  if(ctp_max_cols > ctx_max_cols)
    die("More colours in path files than in graph files!");

  // Check sample names
  size_t pinto, pfrom, ginto, gfrom, i, j;
  StrBuf *pname, *gname;

  // Ugly loops I'm afraid
  // number of paths / graphs and ultimately number of colours should be low
  for(p = 0; p < num_paths; p++) {
    for(i = 0; i < paths[p].fltr.ncols; i++) {
      pinto = path_file_intocol(&paths[p], i);
      pfrom = path_file_fromcol(&paths[p], i);
      pname = &paths[p].hdr.sample_names[pfrom];
      for(g = 0; g < num_graphs; g++) {
        for(j = 0; j < graphs[g].fltr.ncols; j++) {
          ginto = graph_file_intocol(&graphs[g], j);
          gfrom = graph_file_fromcol(&graphs[g], j);
          gname = &graphs[g].hdr.ginfo[gfrom].sample_name;
          if(pinto == ginto && strcmp(pname->buff, gname->buff) != 0) {
            die("Sample names don't match\n%s:%zu%s\n%s:%zu%s\n",
                graphs[g].fltr.orig_path.buff, graphs[g].fltr.cols[j], gname->buff,
                paths[p].fltr.orig_path.buff, paths[p].fltr.cols[i], pname->buff);
          }
        }
      }
    }
  }
}


//
// Thread safe wrapper for path_store [mt for multithreaded]
//

// Returns true if new to colour, false otherwise
// packed points to <PathLen><PackedSeq>
// Returns address of path in PathStore by setting newidx
bool graph_paths_find_or_add_mt(dBNode node, Colour ctpcol,
                                const uint8_t *packed, size_t plen,
                                PathStore *pstore, PathIndex *newidx)
{
  // path_nbytes is the number of bytes in <PackedSeq>
  size_t path_nbytes = (plen+3)/4;

  // debug
  // print_path(hkey, packed, pstore);

  // 1) Get lock for kmer
  // calling bitlock_yield_acquire instead of bitlock_acquire causes
  bitlock_yield_acquire(pstore->kmer_locks, node.key);

  // const PathIndex next = pstore_get_pindex(pstore, node.key);
  const PathIndex next = *(volatile PathIndex*)&pstore->kmer_paths_write[node.key];

  // 2) Search for path
  // PathIndex match = path_store_find(pstore, next, packed, path_nbytes);
  PathIndex match;
  uint8_t *colset;
  size_t phash_pos = 0;
  int pret = path_hash_find_or_insert_mt(&pstore->phash, node.key, packed,
                                         pstore->store, pstore->colset_bytes,
                                         &phash_pos);

  if(pret == 0)
  {
    // Path already exists -> add colour -> release lock
    match = path_hash_get_pindex(&pstore->phash, phash_pos);
    bool added = false;

    if(pstore->num_of_cols > 1)
    {
      // Add colour to bitset
      colset = packedpath_get_colset(pstore->store+match);
      added = !bitset_get(colset, ctpcol);
      if(added) __sync_add_and_fetch((volatile size_t*)&pstore->num_col_paths, 1);
      bitset_set(colset, ctpcol);
    }

    // Relase lock - we're done
    bitlock_release(pstore->kmer_locks, node.key);
    *newidx = match;
    return added;
  }

  // 3) shift store->next to make space
  // Note: if you update mem_bytes, check where it is used - pstore->num_of_bytes
  //       should not include extra mem
  size_t mem_bytes = packedpath_mem2(pstore->colset_bytes, path_nbytes);
  uint8_t *new_path;

  // atomic { new_path = pstore->next; pstore->next += mem_bytes; }
  // __sync_fetch_and_add(x,y) x and y need to be of same type
  new_path = __sync_fetch_and_add((uint8_t * volatile *)(void*)&pstore->next,
                                  (uint8_t *)mem_bytes);

  if(new_path + mem_bytes > pstore->end || pret == -1)
  {
    // DEV: be better, pause threads, dump graph, clear memory
    die("Out of path memory! [%s]", pret == -1 ? "PathHash" : "PathLinkedList");
  }

  PathIndex pindex = (uint64_t)(new_path - pstore->store);

  // 4) Copy new entry
  // Point new entry to existing linked list
  packedpath_set_prev(new_path, next);
  // Update PathHash
  path_hash_set_pindex(&pstore->phash, phash_pos, pindex);

  // bitset
  colset = packedpath_get_colset(new_path);
  memset(colset, 0, pstore->colset_bytes);
  bitset_set(colset, ctpcol);

  // Length + Path
  memcpy(colset+pstore->colset_bytes, packed, sizeof(PathLen) + path_nbytes);

  /* Note: below not true if we are writing to a diff set to walking */
  // path must be written before we move the kmer path pointer forward
  // although there is a write-lock (pstore->kmer_locks), threads currently
  // traversing the graph would fall over
  __sync_synchronize();

  // 5) update kmer pointer
  pstore_set_pindex(pstore, node.key, pindex);

  // Update number of kmers with paths if this the first path for this kmer
  if(next == PATH_NULL) {
    __sync_add_and_fetch((volatile size_t*)&pstore->num_kmers_with_paths, 1);
  }

  // Update number of paths
  __sync_add_and_fetch((volatile size_t*)&pstore->num_of_paths, 1);
  __sync_add_and_fetch((volatile size_t*)&pstore->num_col_paths, 1);
  __sync_add_and_fetch((volatile size_t*)&pstore->num_of_bytes, mem_bytes);

  // status("npaths: %zu nkmers: %zu", pstore->num_of_paths,
  //        pstore->num_kmers_with_paths);

  // 6) release kmer lock
  bitlock_release(pstore->kmer_locks, node.key);

  *newidx = pindex;
  return true;
}

//
// Remove all redundant paths
//

static inline void graph_paths_remove_redundant_node(hkey_t hkey,
                                                     SortedPathSet *set,
                                                     size_t *paths_removed_ptr,
                                                     size_t *bytes_removed_ptr,
                                                     dBGraph *db_graph)
{
  // Construct SortedPathSet
  PathStore *pstore = &db_graph->pstore;
  size_t i, num_orig_entries, num_orig_bytes;
  PathIndex pindex0, pindex1;

  sorted_path_set_init(set, pstore, hkey);

  num_orig_entries = set->members.len;
  num_orig_bytes = sorted_path_get_bytes_sum(set);

  sorted_path_set_slim(set);

  // Update PathStore if set has changed
  if(num_orig_entries != set->members.len)
  {
    ctx_assert(set->members.len > 0 && num_orig_entries > 0);

    pindex0 = pindex1 = set->members.data[0].pindex;
    pstore_set_pindex(pstore, hkey, pindex0);

    for(i = 0; i+1 < set->members.len; i++, pindex0 = pindex1) {
      pindex1 = set->members.data[i+1].pindex;
      packedpath_set_prev(pstore->store + pindex0, pindex1);
    }

    *paths_removed_ptr += num_orig_entries - set->members.len;
    *bytes_removed_ptr += num_orig_bytes - sorted_path_get_bytes_sum(set);
  }
}

typedef struct {
  uint32_t threadid, nthreads;
  dBGraph *db_graph;
  size_t paths_removed, bytes_removed;
} RemoveRedundantPathsJob;

static void graph_paths_remove_redundant_thread(void *arg)
{
  RemoveRedundantPathsJob *job = (RemoveRedundantPathsJob*)arg;
  dBGraph *db_graph = job->db_graph;
  SortedPathSet set;
  sorted_path_set_alloc(&set);
  HASH_ITERATE_PART(&db_graph->ht, job->threadid, job->nthreads,
                    graph_paths_remove_redundant_node,
                    &set, &job->paths_removed, &job->bytes_removed, db_graph);
  sorted_path_set_dealloc(&set);
}

void graph_paths_remove_redundant(dBGraph *db_graph, size_t num_threads)
{
  status("Removing redundant GraphPaths using %zu threads", num_threads);

  size_t i;
  RemoveRedundantPathsJob *jobs;
  jobs = ctx_malloc(num_threads * sizeof(RemoveRedundantPathsJob));

  for(i = 0; i < num_threads; i++) {
    jobs[i] = (RemoveRedundantPathsJob){.threadid = i,
                                        .nthreads = num_threads,
                                        .db_graph = db_graph,
                                        .paths_removed = 0,
                                        .bytes_removed = 0};
  }

  util_run_threads(jobs, num_threads, sizeof(jobs[0]),
                   num_threads, graph_paths_remove_redundant_thread);

  // Update header
  PathStore *pstore = &db_graph->pstore;
  size_t paths_removed = 0, bytes_removed = 0;

  for(i = 0; i < num_threads; i++) {
    paths_removed += jobs[i].paths_removed;
    bytes_removed += jobs[i].bytes_removed;
  }

  // Print
  char num_paths_str[100], old_paths_str[100], new_paths_str[100];
  char num_bytes_str[100], old_bytes_str[100], new_bytes_str[100];

  ulong_to_str(paths_removed, num_paths_str);
  ulong_to_str(pstore->num_of_paths, old_paths_str);
  ulong_to_str(pstore->num_of_paths-paths_removed, new_paths_str);
  ulong_to_str(bytes_removed, num_bytes_str);
  ulong_to_str(pstore->num_of_bytes, old_bytes_str);
  ulong_to_str(pstore->num_of_bytes-bytes_removed, new_bytes_str);

  status("Removed %s paths [%s -> %s] and %s bytes [%s -> %s]",
         num_paths_str, old_paths_str, new_paths_str,
         num_bytes_str, old_bytes_str, new_bytes_str);

  ctx_assert(pstore->num_of_paths >= paths_removed);
  ctx_assert(pstore->num_of_bytes >= bytes_removed);

  pstore->num_of_paths -= paths_removed;
  pstore->num_of_bytes -= bytes_removed;

  ctx_free(jobs);
}

//
// Integrity checks on graph+paths
//

// 1) check node after node has indegree >1 in sample ctxcol
// 2) follow path, check each junction matches up with a node with outdegree >1
// col is graph colour
// packed is just <PackedBases>
bool graph_paths_check_valid(dBNode node, size_t ctxcol, const uint8_t *packed,
                            size_t nbases, const dBGraph *db_graph)
{
  check_ret(db_graph->num_edge_cols == db_graph->num_of_cols ||
            db_graph->node_in_cols != NULL);

  BinaryKmer bkmer;
  Edges edges;
  dBNode nodes[4];
  Nucleotide nucs[4];
  size_t i, j, n, edgecol = db_graph->num_edge_cols > 1 ? ctxcol : 0;
  // length is kmers and junctions
  size_t klen, plen;

  for(klen = 0, plen = 0; plen < nbases; klen++)
  {
    bkmer = db_node_get_bkmer(db_graph, node.key);
    edges = db_node_get_edges(db_graph, node.key, edgecol);

    // Check this node is in this colour
    if(db_graph->node_in_cols != NULL) {
      check_ret(db_node_has_col(db_graph, node.key, ctxcol));
    } else if(db_graph->col_covgs != NULL) {
      check_ret(db_node_get_covg(db_graph, node.key, ctxcol) > 0);
    }

    #ifdef CTXVERBOSE
      char bkmerstr[MAX_KMER_SIZE+1];
      binary_kmer_to_str(bkmer, db_graph->kmer_size, bkmerstr);
      fprintf(stderr, "klen: %zu plen: %zu %zu:%i %s",
             klen, plen, (size_t)node.key, node.orient, bkmerstr);
    #endif

    if(klen == 1) {
      dBNode rnode = db_node_reverse(node);
      Edges backedges = db_node_edges_in_col(rnode, ctxcol, db_graph);
      int outdegree = edges_get_outdegree(backedges, rnode.orient);
      if(outdegree <= 1) {
        fprintf(stderr, "outdegree: %i col: %zu", (int)outdegree, ctxcol);
      }
      check_ret(outdegree > 1);
    }

    n = db_graph_next_nodes(db_graph, bkmer, node.orient,
                            edges, nodes, nucs);

    check_ret(n > 0);

    // Reduce to nodes in our colour if edges limited
    if(db_graph->num_edge_cols == 1 && db_graph->node_in_cols != NULL) {
      for(i = 0, j = 0; i < n; i++) {
        if(db_node_has_col(db_graph, nodes[i].key, ctxcol)) {
          nodes[j] = nodes[i];
          nucs[j] = nucs[i];
          j++;
        }
      }
      n = j; // update number of next nodes
      check_ret(n > 0);
    }

    // If fork check nucleotide
    if(n > 1) {
      Nucleotide expbase = binary_seq_get(packed, plen);

      for(i = 0; i < n && nucs[i] != expbase; i++);
      if(i == n) {
        fprintf(stderr, "plen: %zu expected: %c\n", plen, dna_nuc_to_char(expbase));
        fprintf(stderr, "Got: ");
        for(i = 0; i < n; i++) fprintf(stderr, " %c", dna_nuc_to_char(nucs[i]));
        fprintf(stderr, "\n");
      }
      check_ret(i < n && nucs[i] == expbase);
      node = nodes[i];
      plen++;
    }
    else {
      node = nodes[0];
    }
  }

  return true;
}

static bool packed_path_check(hkey_t hkey, const uint8_t *packed,
                              const GraphPathPairing *gp,
                              const dBGraph *db_graph)
{
  const PathStore *pstore = &db_graph->pstore;
  PathLen len_bases;
  Orientation orient;
  const uint8_t *colset, *seq;
  size_t i;

  colset = packedpath_get_colset(packed);

  seq = packedpath_seq(packed, pstore->colset_bytes);
  packedpath_get_len_orient(packed, pstore->colset_bytes, &len_bases, &orient);

  dBNode node = {.key = hkey, .orient = orient};

  // Check length
  size_t nbytes = sizeof(PathIndex) + pstore->colset_bytes +
                  sizeof(PathLen) + packedpath_len_nbytes(len_bases);

  check_ret(packed + nbytes <= pstore->end);

  // Check at least one colour is set
  uint8_t colset_or = 0;
  for(i = 0; i < pstore->colset_bytes; i++) colset_or |= colset[i];
  check_ret(colset_or != 0);

  // print path
  // print_path(node.key, packed+sizeof(PathIndex)+pstore->colset_bytes, pstore);

  // Check for each colour the path has
  for(i = 0; i < gp->n; i++) {
    if(bitset_get(colset, gp->ctpcols[i])) {
      check_ret(graph_paths_check_valid(node, gp->ctxcols[i], seq, len_bases,
                                        db_graph));
    }
  }

  return true;
}

static bool kmer_check_paths(hkey_t hkey, const GraphPathPairing *gp,
                             const dBGraph *db_graph,
                             size_t *npaths_ptr, size_t *nkmers_ptr)
{
  const PathStore *pstore = &db_graph->pstore;
  PathIndex pindex = pstore_get_pindex(pstore, hkey);
  uint8_t *packed;
  size_t num_paths = 0;

  while(pindex != PATH_NULL)
  {
    packed = pstore->store+pindex;
    check_ret(packed_path_check(hkey, packed, gp, db_graph));
    pindex = packedpath_get_prev(packed);
    num_paths++;
  }

  *npaths_ptr += num_paths;
  *nkmers_ptr += (num_paths > 0);
  return true;
}

// Returns false on first error
bool graph_paths_check_all_paths(const GraphPathPairing *gp,
                                 const dBGraph *db_graph)
{
  size_t num_paths = 0, num_kmers = 0, act_num_paths, act_num_kmers;

  HASH_ITERATE(&db_graph->ht, kmer_check_paths, gp, db_graph,
               &num_paths, &num_kmers);

  act_num_paths = db_graph->pstore.num_of_paths;
  act_num_kmers = db_graph->pstore.num_kmers_with_paths;

  check_ret2(num_paths == act_num_paths, "%zu vs %zu", num_paths, act_num_paths);
  check_ret2(num_kmers == act_num_kmers, "%zu vs %zu", num_kmers, act_num_kmers);
  return true;
}

bool graph_paths_check_path(hkey_t node, PathIndex pindex,
                              const GraphPathPairing *gp,
                              const dBGraph *db_graph)
{
  uint8_t *packed = db_graph->pstore.store+pindex;
  return packed_path_check(node, packed, gp, db_graph);
}

// For debugging
static void _pstore_update_counts(hkey_t hkey, const dBGraph *db_graph,
                                  size_t *nkmers, size_t *npaths,
                                  size_t *nvisited)
{
  const PathStore *pstore = &db_graph->pstore;
  PathIndex pindex = pstore_get_pindex(pstore, hkey);
  size_t n;

  (*nvisited)++;

  if(pindex == PATH_NULL) return;
  (*nkmers)++;

  for(n = 1; (pindex = packedpath_get_prev(pstore->store+pindex)) != PATH_NULL; n++);
  (*npaths) += n;
}

void graph_paths_check_counts(const dBGraph *db_graph)
{
  const PathStore *pstore = &db_graph->pstore;
  size_t nkmers = 0, npaths = 0, nvisited = 0;

  HASH_ITERATE(&db_graph->ht, _pstore_update_counts,
               db_graph, &nkmers, &npaths, &nvisited);

  ctx_assert(nvisited == db_graph->ht.num_kmers);
  ctx_assert(nkmers == pstore->num_kmers_with_paths);
  ctx_assert(npaths == pstore->num_of_paths);
}
