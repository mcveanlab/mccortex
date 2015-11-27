#include "global.h"
#include "commands.h"
#include "util.h"
#include "db_graph.h"
#include "binary_kmer.h"

const char exp_hashtest_usage[] =
"usage: "CMD" hashtest [options] <num_ops>\n"
"\n"
"  Test hash table speed. If threads is set to 0, use single-threaded code.\n"
"\n"
"  -h, --help        This help message\n"
"  -m, --memory <M>  Memory to use\n"
"  -n, --nkmers <N>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T> Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -k, --kmer <K>    Kmer size must be odd ("QUOTE_VALUE(MAX_KMER_SIZE)" >= k >= "QUOTE_VALUE(MIN_KMER_SIZE)")\n"
"  -F, --func-only   Only use the hash function, do not store kmers\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"func-only",    no_argument,       NULL, 'F'},
  {NULL, 0, NULL, 0}
};

struct HashLoopJob {
  dBGraph *db_graph;
  bool single_threaded;
  size_t start, end;
  size_t hash; // return value
};

static inline void hash_loop(void *arg, size_t threadid)
{
  (void)threadid;
  struct HashLoopJob *jptr = (struct HashLoopJob*)arg, j = *jptr;

  size_t i;
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  bool found;
  uint32_t hash = 0;

  if(j.db_graph && j.single_threaded) {
    for(i = j.start; i < j.end; i++) {
      bkmer.b[0] = i;
      hash_table_find_or_insert(&j.db_graph->ht, bkmer, &found);
    }
  } else if(j.db_graph) {
    for(i = j.start; i < j.end; i++) {
      bkmer.b[0] = i;
      hash_table_find_or_insert_mt(&j.db_graph->ht, bkmer, &found,
                                   j.db_graph->bktlocks);
    }
  } else {
    for(i = j.start; i < j.end; i++) {
      bkmer.b[0] = i;
      hash ^= binary_kmer_hash(bkmer, 0);
    }
  }

  jptr->hash = hash;
}

int ctx_exp_hashtest(int argc, char **argv)
{
  size_t nthreads = 0, kmer_size = 0;
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
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'k': cmd_check(!kmer_size,cmd); kmer_size = cmd_uint32_nonzero(cmd, optarg); break;
      case 'F': cmd_check(store_kmers,cmd); store_kmers = false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" hashtest -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  bool single_threaded = false;
  if(nthreads == 0) { single_threaded = true; nthreads = 1; }

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

    db_graph_alloc(&db_graph, kmer_size, 1, 0, kmers_in_hash, DBG_ALLOC_BKTLOCKS);
    hash_table_print_stats(&db_graph.ht);
  }

  status("[threads] using %zu thread%s (%s-threaded code)",
         nthreads, util_plural_str(nthreads),
         single_threaded ? "single" : "multi");

  struct HashLoopJob jobs[nthreads];
  size_t hash = 0;

  for(i = 0; i < nthreads; i++) {
    size_t start = i * (num_ops / nthreads);
    size_t end = (i+1 == nthreads ? num_ops : start + (num_ops / nthreads));
    jobs[i] = (struct HashLoopJob){.db_graph = store_kmers ? &db_graph : NULL,
                                   .single_threaded = single_threaded,
                                   .start = start, .end = end, .hash = 0};
  }

  util_run_threads(jobs, nthreads, sizeof(jobs[0]), nthreads, hash_loop);

  for(i = 0; i < nthreads; i++) hash += jobs[i].hash;

  if(store_kmers) {
    hash_table_print_stats(&db_graph.ht);
    db_graph_dealloc(&db_graph);
  }

  status("Output hash: %zu", hash);

  return EXIT_SUCCESS;
}
