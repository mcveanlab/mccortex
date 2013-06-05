#include "global.h"
#include "db_graph.h"
#include "graph_walker.h"

static void resize_curr_paths(GraphWalker *wlk)
{
  size_t i, old_size = wlk->paths_cap;
  wlk->paths_cap *= 2;
  wlk->curr_paths = realloc(wlk->curr_paths, wlk->paths_cap * sizeof(path_t*));
  wlk->paths_data = realloc(wlk->paths_data, wlk->paths_cap * sizeof(path_t));

  if(wlk->curr_paths == NULL || wlk->paths_data == NULL) die("Out of memory");

  for(i = old_size; i < wlk->paths_cap; i++)
  {
    wlk->curr_paths[i] = wlk->paths_data + i*sizeof(path_t);
    path_alloc(wlk->curr_paths[i]);
  }
}

static void pickup_paths(GraphWalker *wlk)
{
  uint64_t index = db_graph_kmer_path(wlk->db_graph, wlk->node, wlk->orient);

  if(index != PATH_NULL)
  {
    printf("Checking paths list [%zu]\n", (size_t)index);

    const binary_paths_t *paths = &wlk->db_graph->pdata;
  
    // Fetch first path
    size_t i = wlk->num_paths;
    binary_paths_fetch(paths, index, wlk->curr_paths[i]);

    do
    {
      printf("  new possible path\n");
      if(path_has_col(wlk->curr_paths[i], wlk->colour)) {
        printf("  -> picked up path\n");

        // Unmark this paths in the paths database
        if(wlk->num_pp == wlk->pp_cap) {
          wlk->pp_cap *= 2;
          wlk->prev_paths = realloc(wlk->prev_paths, wlk->pp_cap*sizeof(uint64_t));
        }
        wlk->prev_paths[wlk->num_pp++] = index;
        binary_paths_del_col(paths, index, wlk->colour);

        i++;
        if(i+1 >= wlk->paths_cap) resize_curr_paths(wlk);
      }
    }
    while(binary_paths_prev(paths, wlk->curr_paths[i], wlk->curr_paths[i+1]));

    wlk->num_paths = i;
  }
}

void graph_walker_alloc(GraphWalker *wlk)
{
  size_t i, paths_cap = 256;
  path_t **curr_paths = malloc(paths_cap * sizeof(path_t*));
  path_t *paths_data = malloc(paths_cap * sizeof(path_t));

  for(i = 0; i < paths_cap; i++) {
    curr_paths[i] = paths_data + i*sizeof(path_t);
    path_alloc(curr_paths[i]);
  }

  wlk->curr_paths = curr_paths;
  wlk->paths_data = paths_data;
  wlk->paths_cap = paths_cap;

  wlk->pp_cap = 1024;
  wlk->prev_paths = malloc(wlk->pp_cap * sizeof(uint64_t));
}

void graph_walker_dealloc(GraphWalker *wlk)
{
  free(wlk->prev_paths);
  free(wlk->curr_paths);
  free(wlk->paths_data);
}

void graph_walker_init(GraphWalker *wlk, dBGraph *graph, Colour colour,
                       hkey_t node, Orientation orient)
{
  GraphWalker gw = {.db_graph = graph, .colour = colour,
                    .node = node, .orient = orient,
                    .curr_paths = wlk->curr_paths,
                    .paths_data = wlk->paths_data,
                    .num_paths = 0, .paths_cap = wlk->paths_cap,
                    .prev_paths = wlk->prev_paths,
                    .num_pp = 0, .pp_cap = wlk->pp_cap};

  memcpy(wlk, &gw, sizeof(GraphWalker));

  // Get bkmer oriented correctly (not bkey)
  db_graph_oriented_bkmer(graph, node, orient, wlk->bkmer);

  // Load paths from first kmer
  pickup_paths(wlk);
}

void graph_walker_finish(GraphWalker *wlk)
{
  uint64_t *pp_indices = wlk->prev_paths;
  binary_paths_t *paths = &wlk->db_graph->pdata;
  size_t i, num = wlk->num_pp;
  for(i = 0; i < num; i++)
    binary_paths_set_col(paths, pp_indices[i], wlk->colour);
}

// Returns index of choice or -1
int graph_walker_choose(GraphWalker *wlk, size_t num_next, Nucleotide bases[4])
{
  if(num_next == 0 || (num_next > 1 && wlk->num_paths == 0)) return -1;
  if(num_next == 1) return 0;

  // Do all the oldest shades pick a consistent next node?
  path_t *oldest_path = wlk->curr_paths[0];
  size_t greatest_age = oldest_path->pos;
  Nucleotide greatest_nuc = oldest_path->bases[oldest_path->pos];

  size_t i;
  path_t *path;

  for(i = 1; i < wlk->num_paths; i++) {
    path = wlk->curr_paths[i];
    if(path->pos < greatest_age) break;
    if(path->bases[path->pos] != greatest_nuc) return -1;
  }

  // There is unique next node
  // Find the correct next node chosen by the paths
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      return i;

  die("Something went wrong. [path corruption]");
}

// Force a traversal
// If fork is true, node is the result of taking a fork -> slim down paths
void graph_traverse_force(GraphWalker *wlk, hkey_t node, BinaryKmer bkmer,
                          Orientation orient, boolean fork)
{
  if(fork)
  {
    Nucleotide next_base = binary_kmer_last_nuc(bkmer);

    // take all paths that agree with said nucleotide and haven't ended
    // Update ages
    path_t **curr_paths = wlk->curr_paths, *path, *tmp_path;
    size_t i, j;
    int mismatched_path_age = -1;

    for(i = 0, j = 0; i < wlk->num_paths; i++)
    {
      path = curr_paths[i];
      if(path->bases[path->pos] == next_base) {
        path->pos++;
        if(path->pos < path->core.len) {
          SWAP(curr_paths[j], curr_paths[i], tmp_path);
          j++;
        }
      }
      else {
        mismatched_path_age = path->pos;
        break;
      }
    }

    if(mismatched_path_age != -1) {
      while(j > 0 && curr_paths[j]->pos == (unsigned)mismatched_path_age) {
        // drop path j
        j--;
      }
    }

    wlk->num_paths = j;
  }

  wlk->node = node;
  wlk->orient = orient;
  memcpy(wlk->bkmer, bkmer, sizeof(BinaryKmer));

  pickup_paths(wlk);
}

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk)
{
  const dBGraph *db_graph = wlk->db_graph;
  Edges edges = edges_with_orientation(db_graph->edges[wlk->node], wlk->orient);

  hkey_t nodes[4];
  BinaryKmer bkmers[4];
  Orientation orients[4];

  uint8_t num_next = db_graph_next_nodes(wlk->db_graph, wlk->bkmer, edges,
                                         nodes, bkmers, orients);

  return graph_traverse_nodes(wlk, num_next, nodes, bkmers, orients);
}

boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             hkey_t nodes[4], BinaryKmer bkmers[4],
                             Orientation orients[4])
{
  Nucleotide next_bases[4];
  size_t i;

  for(i = 0; i < num_next; i++)
    next_bases[i] = binary_kmer_last_nuc(bkmers[i]);

  int nxt_indx = graph_walker_choose(wlk, num_next, next_bases);
  if(nxt_indx == -1) return false;

  graph_traverse_force(wlk, nodes[nxt_indx], bkmers[nxt_indx], orients[nxt_indx],
                       num_next > 1);

  return true;
}
