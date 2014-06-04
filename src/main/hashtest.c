#include "global.h"

#include "cmd.h"
#include "cmd_mem.h"
#include "util.h"
#include "db_graph.h"
#include "binary_kmer.h"

static const char usage[] =
"usage: hashtest [options] <num_ops>\n"
"  Test hash table speed.  Assume kmer size of "QUOTE_VALUE(MAX_KMER_SIZE)" if none given\n";

int main(int argc, char **argv)
{
  cortex_init();
  ctx_msg_out = stderr;

  CmdArgs args;
  cmd_alloc(&args, argc, argv);

  argv = args.argv;
  argc = args.argc;

  size_t kmer_size = MAX_KMER_SIZE;

  if(argc > 0 && (!strcmp(argv[0],"--kmer_size") || !strcmp(argv[0],"-k"))) {
    if(argc == 1) die("%s <k> requires an argument", argv[0]);
    if(!parse_entire_size(argv[1], &kmer_size) || !(kmer_size&1)) {
      die("Invalid kmer-size (%s %s): requires odd number %i <= k <= %i",
          argv[0], argv[1], MIN_KMER_SIZE, MAX_KMER_SIZE);
    }
    if(kmer_size < MIN_KMER_SIZE || kmer_size > MAX_KMER_SIZE) {
      die("Please recompile with correct kmer size (%zu)", kmer_size);
    }
    argc -= 2;
    argv += 2;
  }

  if(args.argc != 1) print_usage(usage, NULL);

  unsigned long i, num_ops;
  if(!parse_entire_ulong(argv[0], &num_ops))
    print_usage(usage, "Invalid <num_ops>");

  // Decide on memory
  size_t kmers_in_hash, graph_mem;
  kmers_in_hash = cmd_get_kmers_in_hash(&args, 0, num_ops, num_ops, true, &graph_mem);
  cmd_check_mem_limit(args.mem_to_use, graph_mem);

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 0, kmers_in_hash);
  hash_table_print_stats(&db_graph.ht);

  BinaryKmer bkmer;
  bool found;

  for(i = 0; i < num_ops; i++)
  {
    bkmer = binary_kmer_random(kmer_size);
    hash_table_find_or_insert(&db_graph.ht, bkmer, &found);
  }

  hash_table_print_stats(&db_graph.ht);
  db_graph_dealloc(&db_graph);

  cmd_free(&args);
  cortex_destroy();
  return EXIT_SUCCESS;
}
