#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"

#ifdef CTXVERBOSE
#define DEBUG_WALKER 1
#endif

#define USE_COUNTER_PATHS 1

// For GraphWalker to work we assume all edges are merged into one colour
// (i.e. graph->num_edge_cols == 1)
// If only one colour loaded we assume all edges belong to this colour

static void print_path_list(FollowPath **arr, size_t num)
{
  size_t i, j;
  for(i = 0; i < num; i++) {
    FollowPath *path = arr[i];
    printf("   %p %p ", path, path->bases);
    for(j = 0; j < path->len; j++)
      putc(dna_nuc_to_char(path->bases[j]), stdout);
    printf(" [%i/%i]\n", (int)path->pos, (int)path->len);
  }
}

void graph_walker_print_state(const GraphWalker *wlk)
{
  BinaryKmer bkmer = db_node_bkmer(wlk->db_graph, wlk->node);
  char bkmerstr[MAX_KMER_SIZE+1], bkeystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, bkmerstr);
  binary_kmer_to_str(bkmer, wlk->db_graph->kmer_size, bkeystr);
  printf(" GWState:%s (%s:%i)\n", bkmerstr, bkeystr, wlk->orient);
  printf("  num_curr: %zu\n", wlk->num_curr);
  print_path_list(wlk->curr_paths, wlk->num_curr);
  printf("  num_new: %zu\n", wlk->num_new);
  print_path_list(wlk->curr_paths + wlk->num_curr, wlk->num_new);
  printf("  num_counter: %zu\n", wlk->num_counter);
  print_path_list(wlk->counter_paths, wlk->num_counter);
  printf("--\n");
}

static inline void resize_paths(GraphWalker *wlk, PathLen new_len)
{
  size_t i;
  size_t prev_num_paths = wlk->max_num_paths, prev_path_len = wlk->max_path_len;

  // printf("allpaths: %p-%p\n", wlk->allpaths, wlk->allpaths + wlk->max_num_paths);
  // printf("data: %p-%p\n", wlk->data, wlk->data + wlk->max_num_paths * wlk->max_path_len);
  // for(i = 0; i < wlk->num_unused; i++)
  //   printf("unused %p %p\n", wlk->unused_paths[i], wlk->unused_paths[i]->bases);
  // for(i = 0; i < wlk->num_curr+wlk->num_new; i++)
  //   printf("currnt %p %p\n", wlk->curr_paths[i], wlk->curr_paths[i]->bases);
  // for(i = 0; i < wlk->num_counter; i++)
  //   printf("countr %p %p\n", wlk->counter_paths[i], wlk->counter_paths[i]->bases);

  if(wlk->num_unused == 0)
    wlk->max_num_paths *= 2;

  if(new_len > wlk->max_path_len)
    wlk->max_path_len = roundup2pow(new_len);

  // message("RESIZE\n");
  // message(" maxnum %zu -> %zu\n", prev_num_paths, wlk->max_num_paths);
  // message(" maxlen %zu -> %zu\n", prev_path_len, wlk->max_path_len);

  size_t data_mem = wlk->max_num_paths * wlk->max_path_len * sizeof(Nucleotide);
  wlk->data = realloc2(wlk->data, data_mem);

  if(prev_path_len < wlk->max_path_len)
  {
    // Shift data up if needed
    for(i = prev_num_paths-1; i > 0; i--)
    {
      memmove(wlk->data + i * wlk->max_path_len,
              wlk->data + i * prev_path_len,
              prev_path_len * sizeof(Nucleotide));
    }
  }

  for(i = 0; i < prev_num_paths; i++)
    wlk->allpaths[i].bases = wlk->data + i * wlk->max_path_len;

  // for(i = 0; i < prev_num_paths; i++)
  //   printf(" %zu) %p.bases -> %p\n", i, wlk->allpaths+i, wlk->allpaths[i].bases);

  if(prev_num_paths < wlk->max_num_paths)
  {
    size_t paths_mem = wlk->max_num_paths * sizeof(FollowPath);
    size_t ptr_mem = wlk->max_num_paths * sizeof(FollowPath*);

    FollowPath *oldpaths = wlk->allpaths;
    wlk->allpaths = realloc2(wlk->allpaths, paths_mem);
    wlk->unused_paths = realloc2(wlk->unused_paths, ptr_mem);
    wlk->curr_paths = realloc2(wlk->curr_paths, ptr_mem);
    wlk->counter_paths = realloc2(wlk->counter_paths, ptr_mem);

    // Update pointers
    for(i = 0; i < wlk->num_unused; i++)
      wlk->unused_paths[i] = wlk->allpaths + (wlk->unused_paths[i] - oldpaths);

    for(i = 0; i < wlk->num_curr+wlk->num_new; i++)
      wlk->curr_paths[i] = wlk->allpaths + (wlk->curr_paths[i] - oldpaths);

    for(i = 0; i < wlk->num_counter; i++)
      wlk->counter_paths[i] = wlk->allpaths + (wlk->counter_paths[i] - oldpaths);

    // Add new paths to pool of unused paths
    for(i = prev_num_paths; i < wlk->max_num_paths; i++)
    {
      wlk->allpaths[i].bases = wlk->data + i * wlk->max_path_len;
      wlk->unused_paths[wlk->num_unused++] = wlk->allpaths + i;
    }
  }

  // printf("allpaths: %p-%p\n", wlk->allpaths, wlk->allpaths + wlk->max_num_paths);
  // printf("data: %p-%p\n", wlk->data, wlk->data + wlk->max_num_paths * wlk->max_path_len);
  // for(i = 0; i < wlk->num_unused; i++)
  //   printf("unused %p %p\n", wlk->unused_paths[i], wlk->unused_paths[i]->bases);
  // for(i = 0; i < wlk->num_curr+wlk->num_new; i++)
  //   printf("currnt %p %p\n", wlk->curr_paths[i], wlk->curr_paths[i]->bases);
  // for(i = 0; i < wlk->num_counter; i++)
  //   printf("countr %p %p\n", wlk->counter_paths[i], wlk->counter_paths[i]->bases);
  // printf("\n\n");
}

size_t graph_walker_est_mem()
{
  const size_t init_max_path_len = 8;
  const size_t init_max_num_paths = 16;

  size_t mem = init_max_num_paths * init_max_path_len * sizeof(Nucleotide) +
               init_max_num_paths * sizeof(FollowPath) +
               init_max_num_paths * sizeof(FollowPath*) +
               init_max_num_paths * sizeof(FollowPath*) +
               init_max_num_paths * sizeof(FollowPath*);
  return mem;
}

void graph_walker_alloc(GraphWalker *wlk)
{
  // Setup capacity path length and number of paths
  wlk->max_path_len = 8;
  wlk->num_unused = wlk->max_num_paths = 16;

  wlk->data = malloc2(wlk->max_num_paths * wlk->max_path_len * sizeof(Nucleotide));
  wlk->allpaths = malloc2(wlk->max_num_paths * sizeof(FollowPath));
  wlk->unused_paths = malloc2(wlk->max_num_paths * sizeof(FollowPath*));
  wlk->curr_paths = malloc2(wlk->max_num_paths * sizeof(FollowPath*));
  wlk->counter_paths = malloc2(wlk->max_num_paths * sizeof(FollowPath*));

  // message("graph_walker_alloc num_unused: %zu %p\n", wlk->num_unused, wlk->data);

  size_t i;
  for(i = 0; i < wlk->num_unused; i++) {
    wlk->allpaths[i].bases = wlk->data + i * wlk->max_path_len;
    wlk->unused_paths[i] = wlk->allpaths + i;
  }

  wlk->num_curr = wlk->num_new = wlk->num_counter = 0;
  wlk->fork_count = 0;

  // GraphWalker gw = {.db_graph = NULL, .ctxcol = 0, .ctpcol = 0,
  //                   .node = HASH_NOT_FOUND, .orient = FORWARD,
  //                   .bkmer = BINARY_KMER_ZERO_MACRO,
  //                   .data = wlk->data, .allpaths = wlk->allpaths,
  //                   .max_path_len = wlk->max_path_len,
  //                   .max_num_paths = wlk->max_num_paths,
  //                   .unused_paths = wlk->unused_paths,
  //                   .curr_paths = wlk->curr_paths,
  //                   .counter_paths = wlk->counter_paths,
  //                   .num_unused = wlk->num_unused,
  //                   .num_curr = 0, .num_new = 0, .num_counter = 0,
  //                   .fork_count = 0};
  // memcpy(wlk, &gw, sizeof(GraphWalker));
}

void graph_walker_dealloc(GraphWalker *wlk)
{
  free(wlk->data);
  free(wlk->allpaths);
  free(wlk->unused_paths);
  free(wlk->curr_paths);
  free(wlk->counter_paths);
}

// Returns number of paths picked up
static inline size_t pickup_paths(const PathStore *paths, GraphWalker *wlk,
                                  PathIndex index, Orientation orient,
                                  boolean counter)
{
  // printf("pickup %s paths from: %zu:%i\n", counter ? "counter" : "curr",
  //        (size_t)index, orient);

  Orientation porient;
  PathLen len;
  size_t *num = counter ? &wlk->num_counter : &wlk->num_curr;
  size_t start_pos = *num;

  while(index != PATH_NULL)
  {
    packedpack_len_orient(paths->store+index, paths, &len, &porient);
    if(packedpath_has_col(paths->store+index, path, wlk->ctpcol) && orient == porient)
    {
      if(len > wlk->max_path_len || wlk->num_unused == 0)
        resize_paths(wlk, len);

      FollowPath **arr = counter ? wlk->counter_paths : wlk->curr_paths;
      // Take from unused heap
      arr[*num] = wlk->unused_paths[--(wlk->num_unused)];
      path_store_fetch_bases(paths, index, arr[*num]->bases, len);

      // size_t j;
      // for(j = 0; j < len; j++)
      //   putc(dna_nuc_to_char(arr[*num]->bases[j]), stdout);
      // printf(" [len:%zu] %p\n", len, arr[*num]->bases);

      arr[*num]->pos = 0;
      arr[*num]->len = len;
      (*num)++;
    }

    index = packedpath_prev(paths->store+index);
  }

  return *num - start_pos;
}

void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol,
                       hkey_t node, Orientation orient)
{
  // #ifdef DEBUG_WALKER
  //   char str[MAX_KMER_SIZE+1];
  //   binary_kmer_to_str(db_node_bkmer(graph, node), graph->kmer_size, str);
  //   printf("INIT %s:%i\n", str, orient);
  // #endif

  // Check that this walker is fresh.
  // If previously used, ensure graph_walker_finish() was called.
  assert(wlk->num_curr == 0);
  assert(wlk->num_counter == 0);

  // Check that the graph is loaded properly (all edges merged into one colour)
  assert(graph->num_edge_cols == 1);

  GraphWalker gw = {.db_graph = graph, .ctxcol = ctxcol, .ctpcol = ctpcol,
                    .node = node, .orient = orient,
                    .data = wlk->data, .allpaths = wlk->allpaths,
                    .unused_paths = wlk->unused_paths,
                    .curr_paths = wlk->curr_paths,
                    .counter_paths = wlk->counter_paths,
                    .max_path_len = wlk->max_path_len,
                    .max_num_paths = wlk->max_num_paths,
                    .num_unused = wlk->max_num_paths,
                    .num_curr = 0, .num_new = 0, .num_counter = 0,
                    .fork_count = 0};

  memcpy(wlk, &gw, sizeof(GraphWalker));

  // Get bkmer oriented correctly (not bkey)
  wlk->bkmer = db_graph_oriented_bkmer(graph, node, orient);

  // Pick up new paths
  const PathStore *paths = &wlk->db_graph->pdata;
  PathIndex index = db_node_paths(wlk->db_graph, wlk->node);
  wlk->num_new = pickup_paths(paths, wlk, index, wlk->orient, false);
  wlk->num_curr -= wlk->num_new;
}

void graph_walker_finish(GraphWalker *wlk)
{
  // Put all paths back on unused
  size_t i;
  for(i = 0; i < wlk->num_curr+wlk->num_new; i++)
    wlk->unused_paths[wlk->num_unused++] = wlk->curr_paths[i];
  for(i = 0; i < wlk->num_counter; i++)
    wlk->unused_paths[wlk->num_unused++] = wlk->counter_paths[i];
  wlk->num_curr = wlk->num_new = wlk->num_counter = 0;
}

// Hash a path using its length, sequence and current offset/position
static inline uint32_t follow_path_fasthash(const FollowPath *path)
{
  // Hash upto last 32 bases
  uint32_t i, hash = path->len, max = MIN2(path->len, 16);
  for(i = path->len-max; i < path->len; i++) {
    hash ^= path->bases[i];
    hash <<= 2;
  }
  return hash ^ path->pos;
}

// Hash a binary kmer + GraphWalker paths with offsets
uint32_t graph_walker_fasthash(const GraphWalker *wlk, const BinaryKmer bkmer)
{
  uint64_t i, hash64 = bkmer.b[0];

  for(i = 1; i < NUM_BKMER_WORDS; i++)
    hash64 ^= bkmer.b[i];

  // Fold in half, use only bottom 32bits
  uint_fast32_t hash32 = (hash64 ^ (hash64>>32)) & 0xffffffff;

  for(i = 0; i < wlk->num_curr+wlk->num_new; i++)
    hash32 ^= follow_path_fasthash(wlk->curr_paths[i]);

  hash32 = rot32(hash32, 15);

  for(i = 0; i < wlk->num_counter; i++)
    hash32 ^= follow_path_fasthash(wlk->counter_paths[i]);

  hash32 = rot32(hash32, 15);

  return hash32;
}

#define return_step(i,s) { \
  GraphStep _stp = {.idx = (int8_t)(i), .status = (s)}; \
  return _stp; \
}

// Returns index of choice or -1
// Sets is_fork_in_col true if there is a fork in the given colour
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                              const hkey_t next_nodes[4],
                              const Nucleotide next_bases[4])
{
  // #ifdef DEBUG_WALKER
  //   printf("CHOOSE\n");
  //   print_state(wlk);
  // #endif

  if(num_next == 0) return_step(-1, GRPHWLK_NOCOVG);
  if(num_next == 1) return_step( 0, GRPHWLK_FORWARD);

  int8_t indices[4] = {0,1,2,3};
  hkey_t nodes_store[4];
  Nucleotide bases_store[4];
  const hkey_t *nodes = nodes_store;
  const Nucleotide* bases = bases_store;
  size_t i, j;

  // Reduce next nodes that are in this colour
  if(wlk->db_graph->num_of_cols > 1)
  {
    for(i = 0, j = 0; i < num_next; i++)
    {
      if(db_node_has_col(wlk->db_graph, next_nodes[i], wlk->ctxcol)) {
        nodes_store[j] = next_nodes[i];
        bases_store[j] = next_bases[i];
        indices[j] = (int8_t)i;
        j++;
      }
    }

    num_next = j;

    if(num_next == 1) return_step(indices[0], GRPHWLK_COLFWD);
    if(num_next == 0) return_step(-1,         GRPHWLK_NOCOLCOVG);
  }
  else {
    nodes = next_nodes;
    bases = next_bases;
  }

  // We have hit a fork
  // abandon if no path info
  if(wlk->num_curr == 0) return_step(-1, GRPHWLK_NOPATHS);

  // Do all the oldest paths pick a consistent next node?
  FollowPath *oldest_path = wlk->curr_paths[0];
  PathLen greatest_age = oldest_path->pos;
  Nucleotide greatest_nuc = oldest_path->bases[oldest_path->pos];

  FollowPath *path;

  for(i = 1; i < wlk->num_curr; i++) {
    path = wlk->curr_paths[i];
    if(path->pos < greatest_age) break;
    if(path->bases[path->pos] != greatest_nuc)
      return_step(-1, GRPHWLK_SPLIT_PATHS);
  }

  // Does every next node have a path?
  #ifdef USE_COUNTER_PATHS
  size_t c[4] = {0};

  for(i = 0; i < wlk->num_curr && c[0]+c[1]+c[2]+c[3] < num_next; i++)
  {
    path = wlk->curr_paths[i];
    c[path->bases[path->pos]] = 1;
  }

  for(i = 0; i < wlk->num_counter && c[0]+c[1]+c[2]+c[3] < num_next; i++)
  {
    path = wlk->counter_paths[i];
    c[path->bases[path->pos]] = 1;
  }

  // Fail if missing assembly info
  if(c[0]+c[1]+c[2]+c[3] < num_next) return_step(-1, GRPHWLK_MISSING_PATHS);
  if(c[0]+c[1]+c[2]+c[3] > num_next) die("Counter path corruption");
  #endif

  // There is unique next node
  // Find the correct next node chosen by the paths
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      return_step(indices[i], GRPHWLK_USEPATH);

  // If we reach here something has gone wrong
  // print some debug information then exit
  char str[MAX_KMER_SIZE+1];
  binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, str);
  message("Fork: %s\n", str);

  for(i = 0; i < num_next; i++) {
    BinaryKmer bkmer = db_node_bkmer(wlk->db_graph, nodes[i]);
    binary_kmer_to_str(bkmer, wlk->db_graph->kmer_size, str);
    message("  %s [%c]\n", str, dna_nuc_to_char(bases[i]));
  }

  message("curr_paths:\n");
  for(i = 0; i < wlk->num_curr; i++) {
    path = wlk->curr_paths[i];
    message(" %c [%u/%u]\n", dna_nuc_to_char(path->bases[path->pos]),
            path->pos, path->len);
  }

  message("counter_paths:\n");
  for(i = 0; i < wlk->num_counter; i++) {
    path = wlk->counter_paths[i];
    message(" %c [%u/%u]\n", dna_nuc_to_char(path->bases[path->pos]),
            path->pos, path->len);
  }

  message("walker: ctx %zu ctp %zu\n", wlk->ctxcol, wlk->ctpcol);
  message("[path corruption] {%zu:%c}", num_next, dna_nuc_to_char(greatest_nuc));

  die("Did you build this .ctp against THIS EXACT .ctx? (REALLY?)");
}

#undef return_step

// Force a traversal
// If fork is true, node is the result of taking a fork -> slim down paths
// prev is array of nodes with edges to the node we are moving to
void graph_traverse_force_jump(GraphWalker *wlk, hkey_t node, BinaryKmer bkmer,
                               boolean fork)
{
  assert(node != HASH_NOT_FOUND);

  // #ifdef DEBUG_WALKER
  //   char str[MAX_KMER_SIZE+1];
  //   binary_kmer_to_str(bkmer, wlk->db_graph->kmer_size, str);
  //   printf("FORCE JUMP %s (fork:%s)\n", str, fork ? "yes" : "no");
  // #endif

  if(fork)
  {
    // We passed a fork - take all paths that agree with said nucleotide and
    // haven't ended, also update ages
    Nucleotide base = binary_kmer_last_nuc(bkmer);
    FollowPath *path;
    size_t i, j;

    // Check curr paths
    for(i = 0, j = 0; i < wlk->num_curr + wlk->num_new; i++)
    {
      path = wlk->curr_paths[i];
      if(path->bases[path->pos] == base && path->pos+1 < path->len)
      {
        wlk->curr_paths[j++] = path;
        path->pos++;
      }
      else {
        // Put path back in unused pool
        wlk->unused_paths[wlk->num_unused++] = path;
      }
    }
    wlk->num_curr = j;
    wlk->num_new = 0;

    // Check counter paths
    for(i = 0, j = 0; i < wlk->num_counter; i++)
    {
      path = wlk->counter_paths[i];
      if(path->bases[path->pos] == base && path->pos+1 < path->len)
      {
        wlk->counter_paths[j++] = path;
        path->pos++;
      }
      else {
        // Put path back in unused pool
        wlk->unused_paths[wlk->num_unused++] = path;
      }
    }
    wlk->num_counter = j;

    // Statistics
    wlk->fork_count++;
  }

  const dBGraph *db_graph = wlk->db_graph;

  wlk->node = node;
  wlk->bkmer = bkmer;
  wlk->orient = db_node_get_orientation(wlk->bkmer, db_node_bkmer(db_graph, node));

  // Take `new' paths
  wlk->num_curr += wlk->num_new;

  // Pick up new paths
  const PathStore *paths = &db_graph->pdata;
  PathIndex index = db_node_paths(db_graph, wlk->node);
  wlk->num_new = pickup_paths(paths, wlk, index, wlk->orient, false);
  wlk->num_curr -= wlk->num_new;
}

void graph_traverse_force(GraphWalker *wlk, hkey_t node, Nucleotide base,
                          boolean fork)
{
  assert(node != HASH_NOT_FOUND);
  BinaryKmer bkmer;
  bkmer = binary_kmer_left_shift_add(wlk->bkmer, wlk->db_graph->kmer_size, base);
  graph_traverse_force_jump(wlk, node, bkmer, fork);
}

void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    hkey_t prev_nodes[4],
                                    Orientation prev_orients[4],
                                    size_t num_prev)
{
  const PathStore *paths = &wlk->db_graph->pdata;
  PathIndex index;

  // Check that counter paths take this nuc
  Nucleotide next_base = binary_kmer_last_nuc(wlk->bkmer);
  // message("{{ graph_walker_add_counter_paths: %zu }}\n", num_prev);

  size_t i, j, k, num_paths;
  Edges edges;
  FollowPath **new_paths;

  for(i = 0; i < num_prev; i++)
  {
    index = db_node_paths(wlk->db_graph, prev_nodes[i]);
    num_paths = pickup_paths(paths, wlk, index, prev_orients[i], true);
    wlk->num_counter -= num_paths;

    edges = db_node_edges(wlk->db_graph, wlk->ctxcol, prev_nodes[i]);
    if(edges_get_outdegree(edges, prev_orients[i]) > 1) {
      new_paths = wlk->counter_paths + wlk->num_counter;
      for(j = 0, k = 0; j < num_paths; j++) {
        if(new_paths[j]->len > 1 && new_paths[j]->bases[0] == next_base) {
          new_paths[j]->pos++;
          new_paths[k++] = new_paths[j];
        }
        else {
          // Put back in pool of unused paths
          wlk->unused_paths[wlk->num_unused++] = new_paths[j];
        }
      }
      num_paths = k;
    }

    wlk->num_counter += num_paths;
  }
}

void graph_walker_node_add_counter_paths(GraphWalker *wlk, Nucleotide prev_nuc)
{
  const dBGraph *db_graph = wlk->db_graph;
  hkey_t prev_nodes[4];
  Orientation prev_orients[4];
  Nucleotide prev_bases[4];
  size_t i, num_prev_nodes;

  // Use local version if possible
  BinaryKmer bkmer = (wlk->orient == FORWARD ? wlk->bkmer
                                             : db_node_bkmer(db_graph, wlk->node));
  Orientation orient = opposite_orientation(wlk->orient);

  Edges edges = db_node_edges(db_graph, wlk->ctxcol, wlk->node) &
                ~nuc_orient_to_edge(dna_nuc_complement(prev_nuc), orient);

  num_prev_nodes = db_graph_next_nodes(db_graph, bkmer, orient, edges,
                                       prev_nodes, prev_orients, prev_bases);

  // Reverse orientation
  for(i = 0; i < 4; i++) prev_orients[i] = !prev_orients[i];
  graph_walker_add_counter_paths(wlk, prev_nodes, prev_orients, num_prev_nodes);
}

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk)
{
  const dBGraph *db_graph = wlk->db_graph;
  Edges edges = db_node_edges(db_graph, wlk->ctxcol, wlk->node);
  edges = edges_with_orientation(edges, wlk->orient);

  hkey_t nodes[4];
  Orientation orients[4];
  Nucleotide bases[4];
  size_t num_next;

  // char tmpbkmer[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(wlk->bkmer, db_graph->kmer_size, tmpbkmer);
  // printf("bkmer: %s\n", tmpbkmer);

  num_next = db_graph_next_nodes(db_graph, wlk->bkmer, FORWARD, edges,
                                 nodes, orients, bases);

  return graph_traverse_nodes(wlk, num_next, nodes, bases);
}

boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             const hkey_t nodes[4], const Nucleotide bases[4])
{
  wlk->last_step = graph_walker_choose(wlk, num_next, nodes, bases);
  int idx = wlk->last_step.idx;
  if(idx == -1) return false;
  graph_traverse_force(wlk, nodes[idx], bases[idx],
                       graphstep_is_fork(wlk->last_step));
  return true;
}

// Fast traverse - avoid a bkmer_revcmp
static inline void graph_walker_fast(GraphWalker *wlk, const dBNode prev_node,
                                     const dBNode next_node, boolean fork)
{
  const size_t kmer_size = wlk->db_graph->kmer_size;
  BinaryKmer bkmer = db_node_bkmer(wlk->db_graph, next_node.key);
  Nucleotide nuc;

  // Only one path between two nodes
  if(wlk->node == prev_node.key && wlk->orient == prev_node.orient) {
    nuc = db_node_last_nuc(bkmer, next_node.orient, kmer_size);
    graph_traverse_force(wlk, next_node.key, nuc, fork);
  }
  else {
    bkmer = db_node_oriented_bkmer(bkmer, next_node.orient, kmer_size);
    graph_traverse_force_jump(wlk, next_node.key, bkmer, fork);
  }

  // char tmpbkmer[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, tmpbkmer);
  // printf("  forced: %s\n", tmpbkmer);
}

// Fast traversal of a list of nodes using the supplied GraphWalker
// Only visits nodes deemed informative + last node
void graph_walker_fast_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                boolean forward)
{
  if(n == 0) return;

  size_t i;
  boolean infork[3] = {false}, outfork[3] = {false};
  Edges edges;
  dBNode nodes[3];

  edges = db_node_edges(wlk->db_graph, wlk->ctxcol, wlk->node);
  outfork[0] = edges_get_outdegree(edges, wlk->orient) > 1;

  dBNode first_node = {.key = wlk->node, .orient = wlk->orient};
  nodes[0] = first_node;
  nodes[1] = forward ? arr[0] : db_node_reverse(arr[n-1]);

  edges = db_node_edges(wlk->db_graph, wlk->ctxcol, nodes[1].key);
  outfork[1] = edges_get_outdegree(edges, nodes[1].orient) > 1;
  infork[1] = edges_get_indegree(edges, nodes[1].orient) > 1;

  for(i = 0; i+1 < n; i++)
  {
    // printf("i: %zu %zu:%i\n", i, (size_t)nodes[1].key, (int)nodes[1].orient);
    nodes[2] = forward ? arr[i+1] : db_node_reverse(arr[n-i-2]);

    edges = db_node_edges(wlk->db_graph, wlk->ctxcol, nodes[2].key);
    outfork[2] = edges_get_outdegree(edges, nodes[2].orient) > 1;
    infork[2] = edges_get_indegree(edges, nodes[2].orient) > 1;

    // Traverse nodes[i] if:
    // - previous node had out-degree > 1 OR
    // - next node has in-degree > 1
    if(outfork[0] || infork[2]) {
      graph_walker_fast(wlk, nodes[0], nodes[1], outfork[0]);
    }

    // Rotate edges, nodes
    infork[0] = infork[1]; infork[1] = infork[2];
    outfork[0] = outfork[1]; outfork[1] = outfork[2];
    nodes[0] = nodes[1]; nodes[1] = nodes[2];
  }

  // Traverse last node
  graph_walker_fast(wlk, nodes[0], nodes[1], outfork[0]);
}

// Traversal of every node in a list of nodes using the supplied GraphWalker
// Visits each node specifed
void graph_walker_slow_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                boolean forward)
{
  Edges edges;
  boolean prev_fork;
  BinaryKmer bkmer;
  dBNode next;
  size_t i;
  const dBGraph *db_graph = wlk->db_graph;

  for(i = 0; i < n; i++) {
    edges = db_node_edges(db_graph, wlk->ctxcol, wlk->node);
    prev_fork = edges_get_outdegree(edges, wlk->orient) > 1;
    next = forward ? arr[i] : db_node_reverse(arr[n-1-i]);
    bkmer = db_node_bkmer(db_graph, next.key);
    graph_traverse_force_jump(wlk, next.key, bkmer, prev_fork);
  }
}
