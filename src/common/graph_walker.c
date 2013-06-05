#include "global.h"
#include "db_graph.h"
#include "graph_walker.h"

static void resize_curr_paths(GraphWalker *wlk)
{
  size_t i, old_size = wlk->paths_capacity;
  wlk->paths_capacity *= 2;
  wlk->curr_paths = realloc(wlk->curr_paths, wlk->paths_capacity * sizeof(path_t*));
  wlk->paths_data = realloc(wlk->paths_data, wlk->paths_capacity * sizeof(path_t));

  if(wlk->curr_paths == NULL || wlk->paths_data == NULL) die("Out of memory");

  for(i = old_size; i < wlk->paths_capacity; i++)
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
      if(bitset_has(wlk->curr_paths[i]->core.colours, wlk->colour)) {
        printf("  -> picked up path\n");
        i++;
        if(i+1 >= wlk->paths_capacity) resize_curr_paths(wlk);
      }
    }
    while(binary_paths_prev(paths, wlk->curr_paths[i], wlk->curr_paths[i+1]));

    wlk->num_paths = i;
  }
}

void graph_walker_alloc(GraphWalker *wlk)
{
  size_t i, paths_capacity = 256;
  path_t **curr_paths = malloc(paths_capacity * sizeof(path_t*));
  path_t *paths_data = malloc(paths_capacity * sizeof(path_t));

  for(i = 0; i < paths_capacity; i++) {
    curr_paths[i] = paths_data + i*sizeof(path_t);
    path_alloc(curr_paths[i]);
  }

  wlk->curr_paths = curr_paths;
  wlk->paths_data = paths_data;
  wlk->paths_capacity = paths_capacity;
}

void graph_walker_init(GraphWalker *wlk, dBGraph *graph, Colour colour,
                       hkey_t node, Orientation orient)
{
  GraphWalker gw = {.db_graph = graph, .colour = colour,
                    .node = node, .orient = orient,
                    .curr_paths = wlk->curr_paths,
                    .paths_data = wlk->paths_data,
                    .paths_capacity = wlk->paths_capacity,
                    .num_paths = 0};

  memcpy(wlk, &gw, sizeof(GraphWalker));

  // Get bkmer oriented correctly (not bkey)
  db_graph_oriented_bkmer(graph, node, orient, wlk->bkmer);

  // Load paths from first kmer
  pickup_paths(wlk);
}

void graph_walker_dealloc(GraphWalker *wlk)
{
  free(wlk->curr_paths);
  free(wlk->paths_data);
}

// Index of choice or -1
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
  size_t i, j;
  int nxt_indx;

  for(i = 0; i < num_next; i++)
    next_bases[i] = binary_kmer_last_nuc(bkmers[i]);

  nxt_indx = graph_walker_choose(wlk, num_next, next_bases);
  if(nxt_indx == -1) return false;

  if(num_next > 1)
  {
    // num_next > 1, num_paths > 0
    path_t **curr_paths = wlk->curr_paths;

    // Do all the oldest shades pick a consistent next node?
    path_t *oldest_path = curr_paths[0];
    size_t greatest_age = oldest_path->pos;
    Nucleotide greatest_nuc = oldest_path->bases[oldest_path->pos];

    path_t *path, *tmp_path;

    for(i = 1; i < wlk->num_paths && curr_paths[i]->pos == greatest_age; i++);

    // Find index of first path that disagrees
    for(j = i; j < wlk->num_paths; j++)
    {
      path = wlk->curr_paths[j];
      if(path->bases[path->pos] != greatest_nuc)
      {
        // Work backwards to remove paths of matching age
        j--;
        while(j >= i && wlk->curr_paths[j]->pos == path->pos) j--;
        wlk->num_paths = j;
        break;
      }
    }

    // take all paths that agree with said nucleotide and haven't ended
    // Update ages
    for(i = 0, j = 0; i < wlk->num_paths; i++)
    {
      path = wlk->curr_paths[i];
      path->pos++;
      if(path->pos < path->core.len)
        SWAP(wlk->curr_paths[j], wlk->curr_paths[i], tmp_path);
      j++;
    }
    wlk->num_paths = j;
  }

  wlk->node = nodes[nxt_indx];
  wlk->orient = orients[nxt_indx];
  memcpy(wlk->bkmer, bkmers[nxt_indx], sizeof(BinaryKmer));

  // Pick up new paths
  pickup_paths(wlk);

  return true;
}
