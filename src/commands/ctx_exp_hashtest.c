#include "global.h"
#include "commands.h"
#include "util.h"
#include "db_graph.h"
#include "binary_kmer.h"

const char exp_hashtest_usage[] =
"usage: "CMD" hashtest [options] <num_ops>\n"
"\n"
"  Test hash table speed.\n"
"\n"
"  -h, --help        This help message\n"
"  -m, --memory <M>  Memory to use\n"
"  -n, --nkmers <N>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -k, --kmer <K>    Kmer size must be odd ("QUOTE_VALUE(MAX_KMER_SIZE)" >= k >= "QUOTE_VALUE(MIN_KMER_SIZE)")\n"
"  -F, --func-only   Only use the hash function, do not store kmers\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"func-only",    no_argument,       NULL, 'F'},
  {NULL, 0, NULL, 0}
};

int ctx_exp_hashtest(int argc, char **argv)
{
  size_t kmer_size = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool store_kmers = true;

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'k': cmd_check(!kmer_size,cmd); kmer_size = cmd_uint32_nonzero(cmd, optarg); break;
      case 'F': cmd_check(store_kmers,cmd); store_kmers = false; break;
    }
  }

  if(!kmer_size) die("kmer size not set with -k <K>");
  if(kmer_size < MIN_KMER_SIZE || kmer_size > MAX_KMER_SIZE)
    die("Please recompile with correct kmer size (%zu)", kmer_size);
  if(!(kmer_size&1)) {
    die("Invalid kmer-size (%zu): requires odd number %i <= k <= %i",
        kmer_size, MIN_KMER_SIZE, MAX_KMER_SIZE);
  }

  if(optind+1 != argc) cmd_print_usage(NULL);

  size_t i, num_ops;
  if(!parse_entire_size(argv[optind], &num_ops))
    cmd_print_usage("Invalid <num_ops>");

  // Decide on memory
  size_t kmers_in_hash = 0, graph_mem = 0, bits_per_kmer = sizeof(BinaryKmer)*8;
  dBGraph db_graph;

  if(store_kmers)
  {
    // Min and max number of kmers both `num_ops`, since each iterations adds a
    // (probably unique) kmer to the graph
    kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                          memargs.mem_to_use_set,
                                          memargs.num_kmers,
                                          memargs.num_kmers_set,
                                          bits_per_kmer,
                                          num_ops, num_ops,
                                          true, &graph_mem);

    cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

    db_graph_alloc(&db_graph, kmer_size, 1, 0, kmers_in_hash, 0);
    hash_table_print_stats(&db_graph.ht);
  }

  BinaryKmer bkmer;
  bool found;
  uint32_t hash = 0;

  if(store_kmers) {
    for(i = 0; i < num_ops; i++) {
      bkmer = binary_kmer_random(kmer_size);
      hash_table_find_or_insert(&db_graph.ht, bkmer, &found);
    }
  } else {
    for(i = 0; i < num_ops; i++) {
      bkmer = binary_kmer_random(kmer_size);
      hash ^= binary_kmer_hash(bkmer, 0);
    }
  }

  if(store_kmers) {
    hash_table_print_stats(&db_graph.ht);
    db_graph_dealloc(&db_graph);
  } else {
    status("Output hash: %u", hash);
  }

  return EXIT_SUCCESS;
}
