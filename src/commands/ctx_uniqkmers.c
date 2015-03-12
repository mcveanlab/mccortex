#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_checks.h"
#include "build_graph.h"
#include "seqout.h"
#include "seq_reader.h"

const char uniqkmers_usage[] =
"usage: "CMD" uniqkmers [options] <N>\n"
"\n"
"  Generate <N> random unique kmers not in the graph.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -o, --out <bub.fa>      Output file [default: STDOUT]\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
//
"  -k, --kmer <K>          Kmer size (required if only giving --seq input)\n"
"  -g, --graph <in.ctx>    Load kmers from the graph file\n"
"  -1, --seq <in.fa>       Load kmers from a sequence file [SAM/BAM/FASTQ etc.]\n"
"  -F, --flank <in.fa>     Add flanking kmers to <in.fa>\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"graph",        required_argument, NULL, 'g'},
  {"seq",          required_argument, NULL, '1'},
  {"flank",        required_argument, NULL, 'F'},
  {NULL, 0, NULL, 0}
};

static inline dBNode db_graph_add_random_node(dBGraph *db_graph,
                                              BinaryKmer *bkmer_ptr)
{
  size_t i;
  bool found;
  for(i = 0; i < 1000; i++) {
    BinaryKmer bkmer = binary_kmer_random(db_graph->kmer_size);
    dBNode node = db_graph_find_or_add_node(db_graph, bkmer, &found);
    if(!found) {
      *bkmer_ptr = bkmer;
      return node;
    }
  }
  die("Ran 1000 times but couldn't find a unique binary kmer");
}

/**
 * Generate set of kmers from appending bkmer to sequence on given side (left/right).
 * If new kmer keys do not match flank kmer key, add all kmer keys to the graph.
 * Otherwise don't add any kmers and return false.
 * @param left_side  If true add binary kmer to left side of `r`, otherwise add
 *                   to the right side
 * @return           true if addd all kmer keys, false if duplicate kmer key
 */
static inline bool _is_valid_flank(BinaryKmer bkmer, const read_t *r,
                                   bool left_side, dBGraph *db_graph)
{
  if(r->seq.end == 0) return true;

  BinaryKmer bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  BinaryKmer tmp_bkmer = bkmer, tmp_bkey = bkey;
  BinaryKmer bkeys[MAX_KMER_SIZE];
  const size_t kmer_size = db_graph->kmer_size;
  const size_t nkmers = MIN2(kmer_size, r->seq.end);
  Nucleotide nuc;
  size_t i;

  bkeys[0] = tmp_bkey;

  if(left_side) {
    for(i = 1; i < nkmers; i++) {
      nuc = dna_char_to_nuc(r->seq.b[i-1]);
      tmp_bkmer = binary_kmer_left_shift_add(tmp_bkmer, kmer_size, nuc);
      tmp_bkey = binary_kmer_get_key(tmp_bkmer, kmer_size);
      if(binary_kmers_are_equal(bkey, tmp_bkey)) return false;
      bkeys[i] = tmp_bkey;
    }
  }
  else {
    for(i = 1; i < nkmers; i++) {
      nuc = dna_char_to_nuc(r->seq.b[r->seq.end-i]);
      tmp_bkmer = binary_kmer_right_shift_add(tmp_bkmer, kmer_size, nuc);
      tmp_bkey = binary_kmer_get_key(tmp_bkmer, kmer_size);
      if(binary_kmers_are_equal(bkey, tmp_bkey)) return false;
      bkeys[i] = tmp_bkey;
    }
  }

  // Add new kmers
  bool found;
  for(i = 1; i < nkmers; i++)
    hash_table_find_or_insert(&db_graph->ht, bkeys[i], &found);

  return true;
}

static inline void _add_uniq_flanks(read_t *r, const char *path,
                                    FILE *fout, seq_format fmt,
                                    dBGraph *db_graph)
{
  // printf("\n\nADDING FLANKS!\n\n");

  const char *ptr;
  for(ptr = r->seq.b; *ptr; ptr++)
    if(!char_is_acgt(*ptr))
      die("Invalid character in read: '%c' [%s]", *ptr, path);

  size_t i, kmer_size = db_graph->kmer_size;

  int side;
  char bkmerstr[2][MAX_KMER_SIZE+1];

  // side: 0 => left, 1 => right
  // Add to right side, then do left side,
  // in case sequence is shorter than kmer_size
  for(side = 1; side >= 0; side--)
  {
    for(i = 0; i < 100; i++) {
      BinaryKmer bkmer;
      dBNode node = db_graph_add_random_node(db_graph, &bkmer);

      if(_is_valid_flank(bkmer, r, side == 0, db_graph)) {
        binary_kmer_to_str(bkmer, kmer_size, bkmerstr[side]);
        if(side == 1) strm_buf_append_str(&r->seq, bkmerstr[1]);
        break;
      }
      else hash_table_delete(&db_graph->ht, node.key);
    }

    if(i == 100) die("couldn't find unique padding sequence");
  }

  // printf("bkmer0: %s bkmer1: %s\n", bkmerstr[0], bkmerstr[1]);

  if(fmt == SEQ_FMT_PLAIN) {
    fputs(bkmerstr[0], fout);
    fputs(r->seq.b, fout);
  }
  else if(fmt == SEQ_FMT_FASTA) {
    fprintf(fout, ">%s\n%s%s\n", r->name.b, bkmerstr[0], r->seq.b);
  }
  else if(fmt == SEQ_FMT_FASTQ) {
    fprintf(fout, "@%s\n%s%s\n+\n", r->name.b, bkmerstr[0], r->seq.b);
    for(i = 0; i < kmer_size; i++) fputc('.', fout);
    fputs(r->qual.b, fout);
    for(i = r->qual.end; i < r->seq.end; i++) fputc('.', fout);
  }
  else
    die("Unknown format: %i", (int)fmt);
}

int ctx_uniqkmers(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t i, kmer_size = 0;

  GraphFileReader tmp_gfile;
  GraphFileBuffer gfilebuf;
  gfile_buf_alloc(&gfilebuf, 8);

  seq_file_t *tmp_sfile;
  SeqFilePtrBuffer sfilebuf, flankbuf;
  seq_file_ptr_buf_alloc(&sfilebuf, 16);
  seq_file_ptr_buf_alloc(&flankbuf, 16);

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'k': cmd_check(!kmer_size,cmd); kmer_size = cmd_kmer_size(cmd, optarg); break;
      case 'g':
        graph_file_reset(&tmp_gfile);
        graph_file_open2(&tmp_gfile, optarg, "r", 0);
        file_filter_flatten(&tmp_gfile.fltr, 0);
        gfile_buf_push(&gfilebuf, &tmp_gfile, 1);
        break;
      case '1':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&sfilebuf, tmp_sfile);
        break;
      case 'F':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --flank file %s", optarg);
        seq_file_ptr_buf_add(&flankbuf, tmp_sfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" uniqkmers -h` for help. Bad option: %s", argv[optind-1]);
      default: die("Bad option: %s", cmd);
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(optind+1 != argc) cmd_print_usage("Expected <N> for number of kmers only");

  size_t num_uniqkmers = 0;

  if(!parse_entire_size(argv[argc-1], &num_uniqkmers))
    cmd_print_usage("Invalid number of unique kmers: %s", argv[argc-1]);

  if(gfilebuf.len == 0 && !kmer_size)
    die("kmer size not set with -k <K>");
  else if(!kmer_size)
    kmer_size = gfilebuf.b[0].hdr.kmer_size;
  else if(gfilebuf.len > 0 && gfilebuf.b[0].hdr.kmer_size != kmer_size) {
    die("Kmer size mismatches: %zu vs %zu (you do not have to specify -k)",
        (size_t)gfilebuf.b[0].hdr.kmer_size, kmer_size);
  }

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfilebuf.b, gfilebuf.len, NULL, 0, -1);

  bool ctx_uses_stdin = false;
  size_t ctx_max_kmers = num_uniqkmers, ctx_sum_kmers = num_uniqkmers;

  for(i = 0; i < gfilebuf.len; i++) {
    ctx_max_kmers = MAX2(ctx_max_kmers, graph_file_nkmers(&gfilebuf.b[i]));
    ctx_sum_kmers += graph_file_nkmers(&gfilebuf.b[i]);
    ctx_uses_stdin |= file_filter_isstdin(&gfilebuf.b[i].fltr);
  }

  if(ctx_uses_stdin) ctx_sum_kmers = SIZE_MAX;

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  FILE *fout = futil_fopen_create(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 1, kmers_in_hash, DBG_ALLOC_BKTLOCKS);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < gfilebuf.len; i++) {
    graph_load(&gfilebuf.b[i], gprefs, &stats);
    graph_file_close(&gfilebuf.b[i]);
    gprefs.empty_colours = false;
  }

  hash_table_print_stats(&db_graph.ht);

  // Load sequence files into colour zero
  build_graph_from_seq(&db_graph, sfilebuf.b, sfilebuf.len, nthreads, 0);
  build_graph_from_seq(&db_graph, flankbuf.b, flankbuf.len, nthreads, 0);

  if(sfilebuf.len > 0 || flankbuf.len > 0)
    hash_table_print_stats(&db_graph.ht);

  // Print random kmers
  BinaryKmer bkmer;
  char bkmerstr[MAX_KMER_SIZE+1];

  seq_format fmt = SEQ_FMT_FASTA;

  // Add random kmers to flank input sequences
  if(flankbuf.len > 0)
  {
    seq_file_t *sf;
    read_t r;
    seq_read_alloc(&r);

    for(i = 0; i < flankbuf.len; i++)
    {
      // Copy path in case we cannot reopen sequence file
      char *path = strdup(flankbuf.b[i]->path);
      sf = flankbuf.b[i] = seq_reopen(flankbuf.b[i]);
      if(sf == NULL) die("Couldn't reopen file: %s", path);
      free(path);

      while(seq_read(sf, &r) > 0) {
        _add_uniq_flanks(&r, sf->path, fout, fmt, &db_graph);
      }
      seq_close(sf);
    }
    seq_read_dealloc(&r);
  }

  // Generate random kmers not in the graph
  char sname[100];
  for(i = 0; i < num_uniqkmers; i++)
  {
    db_graph_add_random_node(&db_graph, &bkmer);
    binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
    sprintf(sname, ">kmer%zu", i);
    seqout_print_strs(sname, bkmerstr, kmer_size, NULL, 0, fmt, fout);
  }

  char num_kmers_str[100];
  ulong_to_str(num_uniqkmers, num_kmers_str);
  status("  wrote %s kmers to: %s\n", num_kmers_str, futil_outpath_str(out_path));
  fclose(fout);

  gfile_buf_dealloc(&gfilebuf);
  seq_file_ptr_buf_dealloc(&sfilebuf);
  seq_file_ptr_buf_dealloc(&flankbuf);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
