#include "global.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"

#ifdef DEBUG
#define DEBUG_WALKER 1
#endif

#ifdef DEBUG_WALKER
static void print_state(const GraphWalker *wlk)
{
  char tmp[100], tmp2[100];
  binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, tmp);
  ConstBinaryKmerPtr bkmerptr = db_node_bkmer(wlk->db_graph, wlk->node);
  binary_kmer_to_str(bkmerptr, wlk->db_graph->kmer_size, tmp2);
  printf("gw:%s (%s:%i)\n", tmp, tmp2, wlk->orient);
  printf(" num_curr_paths: %zu; path_cap: %zu; num_new_paths: %zu;\n",
         wlk->num_curr_paths, wlk->paths_cap, wlk->num_new_paths);
  size_t i;
  for(i = 0; i < wlk->num_curr_paths; i++) {
    binary_paths_dump_path(wlk->curr_paths[i]);
  }
  printf("-\n");
  for(i = 0; i < wlk->num_new_paths; i++) {
    binary_paths_dump_path(wlk->curr_paths[wlk->num_curr_paths+i]);
  }
  printf("-\n");
}
#endif

static void resize_curr_paths(GraphWalker *wlk)
{
  #ifdef DEBUG_WALKER
    printf("%s:%i RESIZE %zu\n", __FILE__, __LINE__, wlk->paths_cap);
    print_state(wlk);
  #endif

  size_t i, old_size = wlk->paths_cap;
  path_t *old_paths_data = wlk->paths_data;
  wlk->paths_cap *= 2;
  wlk->curr_paths = realloc(wlk->curr_paths, wlk->paths_cap * sizeof(path_t*));
  wlk->paths_data = realloc(wlk->paths_data, wlk->paths_cap * sizeof(path_t));

  if(wlk->curr_paths == NULL || wlk->paths_data == NULL) die("Out of memory");

  // Update old addresses
  for(i = 0; i < old_size; i++)
    wlk->curr_paths[i] = wlk->paths_data + (wlk->curr_paths[i] - old_paths_data);

  // alloc new
  for(i = old_size; i < wlk->paths_cap; i++)
  {
    wlk->curr_paths[i] = wlk->paths_data + i;
    path_alloc(wlk->curr_paths[i], wlk->db_graph->pdata.num_of_cols);
  }
}

static void pickup_paths(GraphWalker *wlk)
{
  wlk->num_curr_paths += wlk->num_new_paths;
  wlk->num_new_paths = 0;
  size_t min_len = 0, end_pos = wlk->num_curr_paths;

  if(wlk->num_curr_paths > 0)
    min_len = wlk->curr_paths[0]->len - wlk->curr_paths[0]->pos;

  if(end_pos+1 >= wlk->paths_cap)
    resize_curr_paths(wlk);

  // Pick up new paths
  uint64_t prev_index, index = db_node_paths(wlk->db_graph, wlk->node, wlk->orient);
  const binary_paths_t *paths = &wlk->db_graph->pdata;

  while(index != PATH_NULL)
  {
    path_t *nxtpath = wlk->curr_paths[end_pos];
    binary_paths_fetch(paths, index, nxtpath);
    prev_index = nxtpath->prev;

    if(path_has_col(nxtpath, wlk->colour) && nxtpath->len > min_len) {
      end_pos++;
      if(end_pos+1 >= wlk->paths_cap) resize_curr_paths(wlk);
    }

    index = prev_index;
  }

  wlk->num_new_paths = end_pos - wlk->num_curr_paths;
}

void graph_walker_alloc(GraphWalker *wlk, uint32_t num_of_cols)
{
  size_t i, curr_cap = 256, counter_cap = 256;
  size_t paths_cap = curr_cap + counter_cap;

  path_t **curr_paths = malloc(curr_cap * sizeof(path_t*));
  path_t **counter_paths = malloc(counter_cap * sizeof(path_t*));
  path_t *paths_data = malloc(paths_cap * sizeof(path_t));

  for(i = 0; i < curr_cap; i++) {
    curr_paths[i] = paths_data + i;
    path_alloc(curr_paths[i], num_of_cols);
  }

  for(i = curr_cap; i < paths_cap; i++) {
    counter_paths[i] = paths_data + i;
    path_alloc(counter_paths[i], num_of_cols);
  }

  wlk->curr_paths = curr_paths;
  wlk->counter_paths = counter_paths;
  wlk->paths_data = paths_data;
  wlk->paths_cap = paths_cap;
}

void graph_walker_dealloc(GraphWalker *wlk)
{
  size_t i;
  path_t *paths = wlk->paths_data;
  for(i = 0; i < wlk->paths_cap; i++) path_dealloc(paths+i);
  free(wlk->curr_paths);
  free(wlk->paths_data);
}

void graph_walker_init(GraphWalker *wlk, const dBGraph *graph, Colour colour,
                       hkey_t node, Orientation orient)
{
  #ifdef DEBUG_WALKER
    char str[100];
    binary_kmer_to_str(db_node_bkmer(graph, node), graph->kmer_size, str);
    printf("INIT %s:%i\n", str, orient);
  #endif

  GraphWalker gw = {.db_graph = graph, .colour = colour,
                    .node = node, .orient = orient,
                    .curr_paths = wlk->curr_paths,
                    .paths_data = wlk->paths_data,
                    .num_curr_paths = 0, .num_counter_paths = 0,
                    .num_new_paths = 0,
                    .paths_cap = wlk->paths_cap};

  memcpy(wlk, &gw, sizeof(GraphWalker));

  // Get bkmer oriented correctly (not bkey)
  db_graph_oriented_bkmer(graph, node, orient, wlk->bkmer);

  // Load paths from first kmer
  pickup_paths(wlk);
}

// Gets context up to (but not including) the node you pass
// graph_walker_init_context(wlk, graph, colour, node, orient)
// graph_traverse_force_jump(node, orient)
// Remember to call finish when done with wlk
void graph_walker_init_context(GraphWalker *wlk, const dBGraph *db_graph,
                               uint64_t *visited, Colour colour,
                               hkey_t node, Orientation orient)
{
  #ifdef DEBUG_WALKER
    char str[100];
    binary_kmer_to_str(db_node_bkmer(db_graph, node), db_graph->kmer_size, str);
    printf("INIT CONTEXT %s:%i\n", str, orient);
  #endif

  hkey_t nodes[MAX_WALK_BACK_NODES+1];
  Orientation orients[MAX_WALK_BACK_NODES+1];
  Nucleotide bases[MAX_WALK_BACK_NODES+1];

  Orientation dir = rev_orient(orient);

  graph_walker_init(wlk, db_graph, colour, node, dir);
  db_node_set_traversed(visited, node, dir);

  nodes[0] = node;
  orients[0] = dir;
  bases[0] = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);

  int i, num_nodes = 1;

  while(num_nodes <= MAX_WALK_BACK_NODES && graph_traverse(wlk) &&
        !db_node_has_traversed(visited, wlk->node, wlk->orient))
  {
    db_node_set_traversed(visited, wlk->node, wlk->orient);
    nodes[num_nodes] = wlk->node;
    orients[num_nodes] = wlk->orient;
    bases[num_nodes] = binary_kmer_first_nuc(wlk->bkmer, db_graph->kmer_size);
    num_nodes++;
  }

  graph_walker_finish(wlk);

  // Remove marks on all kmers
  for(i = 0; i < num_nodes; i++)
    db_node_fast_clear_traversed(visited, nodes[i]);

  graph_walker_init(wlk, db_graph, colour,
                    nodes[num_nodes-1], rev_orient(orients[num_nodes-1]));

  // Walk back over the kmers
  boolean was_fork = false;
  Edges edges;
  Nucleotide base;
  for(i = num_nodes-1; i > 0; ) {
    edges = db_node_edges(db_graph, nodes[i]);
    was_fork = edges_get_indegree(edges, orients[i]) > 1;
    i--;
    base = binary_nuc_complement(bases[i]);
    #ifdef DEBUG_WALKER
      printf(" take base %c\n", binary_nuc_to_char(base));
    #endif
    graph_traverse_force(wlk, nodes[i], base, was_fork);
  }

  #ifdef DEBUG_WALKER
    printf("GOT CONTEXT\n\n");
    print_state(wlk);
  #endif
}

void graph_walker_finish(GraphWalker *wlk)
{
  (void)wlk;
  // uint64_t *pp_indices = wlk->prev_paths;
  // binary_paths_t *paths = &wlk->db_graph->pdata;
  // size_t i, num = wlk->num_pp;
  // for(i = 0; i < num; i++)
  //   binary_paths_set_col(paths, pp_indices[i], wlk->colour);
}

// Returns index of choice or -1
int graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                        const hkey_t next_nodes[4],
                        const Nucleotide next_bases[4])
{
  // #ifdef DEBUG_WALKER
  //   printf("CHOOSE\n");
  //   print_state(wlk);
  // #endif

  if(num_next == 0) return -1;
  if(num_next == 1) return 0;

  // Reduce next nodes that are
  int indices[4];
  hkey_t nodes[4];
  Nucleotide bases[4];
  size_t i, j;

  for(i = 0, j = 0; i < num_next; i++)
  {
    if(db_node_has_col(wlk->db_graph, next_nodes[i], wlk->colour)) {
      nodes[j] = next_nodes[i];
      bases[j] = next_bases[i];
      indices[j] = i;
      j++;
    }
  }

  if(j == 0)
  {
    // Take all options
    memcpy(nodes, next_nodes, 4 * sizeof(hkey_t));
    memcpy(bases, next_bases, 4 * sizeof(Nucleotide));
    for(i = 0; i < 4; i++) indices[i] = i;
  }
  else {
    num_next = j;
  }

  if(num_next > 1 && wlk->num_curr_paths == 0) return -1;
  if(num_next == 1) return indices[0];

  // Do all the oldest shades pick a consistent next node?
  path_t *oldest_path = wlk->curr_paths[0];
  size_t greatest_age = oldest_path->pos;
  Nucleotide greatest_nuc = oldest_path->bases[oldest_path->pos];

  path_t *path;

  for(i = 1; i < wlk->num_curr_paths; i++) {
    path = wlk->curr_paths[i];
    if(path->pos < greatest_age) break;
    if(path->bases[path->pos] != greatest_nuc) return -1;
  }

  // There is unique next node
  // Find the correct next node chosen by the paths
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      return indices[i];

  die("Something went wrong. [path corruption] {%zu:%c}",
      num_next, binary_nuc_to_char(greatest_nuc));
}

// Force a traversal
// If fork is true, node is the result of taking a fork -> slim down paths
// prev is array of nodes with edges to the node we are moving to
void graph_traverse_force_jump(GraphWalker *wlk, hkey_t node, BinaryKmer bkmer,
                               boolean fork)
{
  #ifdef DEBUG_WALKER
    char str[100];
    binary_kmer_to_str(bkmer, wlk->db_graph->kmer_size, str);
    printf("FORCE JUMP %s (fork:%s)\n", str, fork ? "yes" : "no");
  #endif

  if(fork)
  {
    Nucleotide base = binary_kmer_last_nuc(bkmer);

    // take all paths that agree with said nucleotide and haven't ended
    // Update ages
    path_t **curr_paths = wlk->curr_paths, *path, *tmp_path;
    size_t i, j;
    // int mismatched_path_age = -1;

    for(i = 0, j = 0; i < wlk->num_curr_paths + wlk->num_new_paths; i++)
    {
      path = curr_paths[i];
      if(path->bases[path->pos] == base && path->pos+1 < path->len)
      {
        path->pos++;
        SWAP(curr_paths[j], curr_paths[i], tmp_path);
        j++;
      }
    }

    wlk->num_curr_paths = j;
    wlk->num_new_paths = 0;
  }

  const dBGraph *db_graph = wlk->db_graph;

  wlk->node = node;
  binary_kmer_assign(wlk->bkmer, bkmer);
  wlk->orient = db_node_get_orientation(wlk->bkmer, db_node_bkmer(db_graph, node));

  pickup_paths(wlk);
}

void graph_traverse_force(GraphWalker *wlk, hkey_t node, Nucleotide base,
                          boolean fork)
{
  BinaryKmer bkmer;
  binary_kmer_assign(bkmer, wlk->bkmer);
  binary_kmer_left_shift_add(bkmer, wlk->db_graph->kmer_size, base);
  graph_traverse_force_jump(wlk, node, bkmer, fork);
}

// void graph_walker_add_counter_paths(GraphWalker *wlk, hkey_t node,
//                                     hkey_t prev_nodes[4],
//                                     Orientation prev_orients[4],
//                                     size_t num_prev)
// {
//   size_t i;
//   for(i = 0; i < num_prev; i++)
//   {
//     // Get paths in orientation that go to node
//     // Add to counter_paths
//   }
// }

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk)
{
  const dBGraph *db_graph = wlk->db_graph;
  Edges edges = edges_with_orientation(db_graph->edges[wlk->node], wlk->orient);

  hkey_t nodes[4];
  BinaryKmer bkmers[4];
  Nucleotide bases[4];
  size_t i, num_next;

  num_next = db_graph_next_nodes(wlk->db_graph, wlk->bkmer, edges, nodes, bkmers);

  for(i = 0; i < num_next; i++)
    bases[i] = binary_kmer_last_nuc(bkmers[i]);

  return graph_traverse_nodes(wlk, num_next, nodes, bases);
}

boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             const hkey_t nodes[4], const Nucleotide bases[4])
{
  int nxt_indx = graph_walker_choose(wlk, num_next, nodes, bases);
  // printf("  nxt_indx: %i\n", nxt_indx);
  if(nxt_indx == -1) return false;
  graph_traverse_force(wlk, nodes[nxt_indx], bases[nxt_indx], num_next > 1);
  return true;
}
