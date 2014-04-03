#include "global.h"
#include "all_tests.h"
#include "build_graph.h"
#include "graph_walker.h"
#include "repeat_walker.h"

static void test_walk(GraphWalker *gwlk, RepeatWalker *rptwlk,
                      dBNode node0, dBNodeBuffer *nbuf,
                      const dBGraph *graph,
                      size_t expnkmers, const char *ans)
{
  db_node_buf_reset(nbuf);
  graph_walker_init(gwlk, graph, 0, 0, node0);

  do {
    db_node_buf_add(nbuf, gwlk->node);
  }
  while(graph_traverse(gwlk) && rpt_walker_attempt_traverse(rptwlk, gwlk));

  // db_nodes_print(nbuf->data, nbuf->len, graph, stdout);
  // printf("\n");
  // printf("%s\n", graph_step_str[gwlk->last_step.status]);

  TASSERT2(nbuf->len == expnkmers, "%zu / %zu", nbuf->len, expnkmers);

  char tmp[nbuf->len+MAX_KMER_SIZE];
  db_nodes_to_str(nbuf->data, nbuf->len, graph, tmp);
  TASSERT2(strcmp(tmp,ans) == 0, "%s vs %s", tmp, ans);

  graph_walker_finish(gwlk);
  rpt_walker_fast_clear(rptwlk, nbuf->data, nbuf->len);
}

static void test_repeat_loop()
{
  TASSERT(sizeof(FollowPath) == 20);

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = 0, .ins_gap_max = 0,
                            .one_way_gap_traverse = true, .use_end_check = true,
                            .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  // Sequence with repeat
  char seq[] = "ATTTGGAACTCCGGA"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "GATAGGGCCAGT"
               "CGTCAGGAGCTAACT";

  char p0[] = "ATTTGGAACTCCGGA""GATAGGGCCAGT""GATAGGGCCAGT";
  char p1[] = "GATAGGGCCAGT""GATAGGGCCAGT""CGTCAGGAGCTAACT";

  // Allocate graph, but don't add any sequence
  _construct_graph_with_paths(&graph, kmer_size, ncols, NULL, 0, params);

  GenPathWorker *gen_path_wrkr = gen_paths_workers_alloc(1, &graph, NULL);

  GraphWalker gwlk;
  RepeatWalker rptwlk;
  graph_walker_alloc(&gwlk);
  rpt_walker_alloc(&rptwlk, graph.ht.capacity, 12);

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 1024);

  // Construct graph but no paths
  build_graph_from_str_mt(&graph, 0, seq, strlen(seq));
  TASSERT2(graph.ht.num_kmers == 15+12+15, "%zu", (size_t)graph.ht.num_kmers);

  // Find first node in sequence
  dBNode node0 = db_graph_find_str(&graph, seq);
  TASSERT(node0.key != HASH_NOT_FOUND);

  // 1) With no paths
  char ans0[] = "ATTTGGAACTCCGGA""GATAGGGCCAGT";
  test_walk(&gwlk, &rptwlk, node0, &nbuf, &graph, 15+2, ans0);

  // 2) Add small paths - produces collapsed down seq with two copy repeat
  gen_paths_from_str_mt(gen_path_wrkr, p0, params);
  gen_paths_from_str_mt(gen_path_wrkr, p1, params);
  char ans1[] = "ATTTGGAACTCCGGA""GATAGGGCCAGT""GATAGGGCCAGT""CGTCAGGAGCTAACT";
  test_walk(&gwlk, &rptwlk, node0, &nbuf, &graph, 15+12+12+5, ans1);

  // 3) Add long paths
  gen_paths_from_str_mt(gen_path_wrkr, seq, params);
  test_walk(&gwlk, &rptwlk, node0, &nbuf, &graph, strlen(seq)+1-kmer_size, seq);

  graph_walker_dealloc(&gwlk);
  rpt_walker_dealloc(&rptwlk);
  db_node_buf_dealloc(&nbuf);
  gen_paths_workers_dealloc(gen_path_wrkr, 1);
  _deconstruct_graph_with_paths(&graph);
}

void test_repeat_walker()
{
  test_status("Testing repeat_walker.h");
  test_repeat_loop();
}
