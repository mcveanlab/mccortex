#include "global.h"

#include "seq_file.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "seq_reader.h"
#include "graph_format.h"

const char reads_usage[] =
"usage: "CMD" reads [options] <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Filters reads based on which have a kmer in the graph. \n"
"\n"
"  -m, --memory <mem>          Memory to use\n"
"  -n, --nkmers <kmers>        Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -f, --fasta                 Output as gzipped FASTA\n"
"  -q, --fastq                 Output as gzipped FASTQ [default]\n"
"  -v, --invert                Print reads/read pairs with no kmer in graph\n"
"  -1, --seq <in> <O>          Writes output to <O>.fq.gz\n"
"  -2, --seq2 <in1> <in2> <O>  Writes output to <O>.{1,2}.fq.gz\n"
"  -i, --seqi <in> <O>         Writes output to <O>.{1,2}.fq.gz\n"
"\n"
"  Can specify --seq/--seq2 multiple times.\n"
"\n";

typedef struct {
  dBGraph *const db_graph;
  LoadingStats *stats;
  char *in1, *in2;
  gzFile out1, out2;
  size_t num_of_reads_printed;
  void (*print)(const read_t *r, gzFile gz, size_t linewrap);
  bool invert;
} AlignReadsData;

static void get_out_path(char *path, size_t len, bool use_fq, int pe_num)
{
  if(pe_num == 0) {
    memcpy(path+len, use_fq ? ".fq.gz" : ".fa.gz", 6);
    path[len+6] = '\0';
  }
  else
    sprintf(path+len, use_fq ? ".%i.fq.gz" : ".%i.fa.gz", pe_num);
}

static void check_outfile_exists(char *outbase, bool is_pe, bool use_fq)
{
  size_t len = strlen(outbase);
  char path[len+1+8]; // .1.fq.gz
  memcpy(path, outbase, len);
  path[len] = '\0';

  if(is_pe) {
    get_out_path(path, len, use_fq, 1);
    if(futil_file_exists(path)) die("Output file already exists: %s", path);
    if(!futil_is_file_writable(path)) cmd_print_usage("Cannot write: %s", path);
    get_out_path(path, len, use_fq, 2);
    if(futil_file_exists(path)) die("Output file already exists: %s", path);
    if(!futil_is_file_writable(path)) cmd_print_usage("Cannot write: %s", path);
  }
  else {
    get_out_path(path, len, use_fq, 0);
    if(futil_file_exists(path)) die("Output file already exists: %s", path);
    if(!futil_is_file_writable(path)) cmd_print_usage("Cannot write: %s", path);
  }
}

static hkey_t find_node(BinaryKmer bkmer, const dBGraph *db_graph)
{
  BinaryKmer bkey = bkmer_get_key(bkmer, db_graph->kmer_size);
  return hash_table_find(&db_graph->ht, bkey);
}

static bool read_touches_graph(const read_t *r, const dBGraph *db_graph,
                                  LoadingStats *stats)
{
  bool found = false;
  size_t kmer_size = db_graph->kmer_size, num_contigs = 0, num_kmers_loaded = 0;

  if(r->seq.end >= kmer_size)
  {
    size_t search_pos = 0, start, end = 0, i;
    BinaryKmer bkmer; Nucleotide nuc;

    while((start = seq_contig_start(r, search_pos, kmer_size, 0,0)) < r->seq.end &&
          !found)
    {
      end = seq_contig_end(r, start, kmer_size, 0, 0, &search_pos);
      stats->total_bases_loaded += end - start;
      num_contigs++;

      bkmer = binary_kmer_from_str(r->seq.b + start, kmer_size);
      num_kmers_loaded++;
      if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }

      for(i = start+kmer_size; i < end; i++)
      {
        nuc = dna_char_to_nuc(r->seq.b[i]);
        bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
        num_kmers_loaded++;
        if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }
      }
    }
  }

  // Update stats
  stats->total_bases_read += r->seq.end;
  stats->num_kmers_loaded += num_kmers_loaded;
  stats->num_kmers_novel += num_kmers_loaded - found;
  if(num_contigs > 0) stats->num_good_reads++;
  else stats->num_bad_reads++;

  return found;
}

void filter_reads(read_t *r1, read_t *r2,
                  uint8_t qoffset1, uint8_t qoffset2, void *ptr)
{
  (void)qoffset1; (void)qoffset2;

  AlignReadsData *data = (AlignReadsData*)ptr;
  const dBGraph *db_graph = data->db_graph;
  LoadingStats *stats = data->stats;

  bool touches_graph = read_touches_graph(r1, db_graph, stats) ||
                          (r2 != NULL && read_touches_graph(r2, db_graph, stats));

  if(touches_graph != data->invert)
  {
    if(r2 != NULL) {
      // Print paired-end
      gzFile out2 = data->out2 != NULL ? data->out2 : data->out1;
      data->print(r1, data->out1, 0);
      data->print(r2, out2, 0);
    }
    else {
      // Print single-ended
      data->print(r1, data->out1, 0);
    }
    data->num_of_reads_printed++;
  }
}

int ctx_reads(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that there are at least 4 arguments

  // Test other input args
  // Check filelists are readable
  // Check output is writable

  bool use_fq = false, use_fa = false, invert = false;
  seq_file_t *seqfiles[argc];
  size_t num_sf = 0, sf = 0;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++)
  {
    if(!strcmp(argv[argi], "--fastq") || !strcmp(argv[argi],"-q")) use_fq = true;
    else if(!strcmp(argv[argi], "--fasta") || !strcmp(argv[argi],"-f")) use_fa = true;
    else if(!strcmp(argv[argi], "--invert") || !strcmp(argv[argi],"-v")) invert = true;
    else if(!strcmp(argv[argi], "--seq") || !strcmp(argv[argi],"-1"))
    {
      if(argi + 2 >= argc) cmd_print_usage("Missing arguments");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file: %s", argv[argi+1]);
      argi += 2;
    }
    else if(!strcmp(argv[argi], "--seq2") || !strcmp(argv[argi],"-2"))
    {
      if(argi + 4 >= argc) cmd_print_usage("Missing arguments");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read first --seq2 file: %s", argv[argi+1]);
      if((seqfiles[num_sf++] = seq_open(argv[argi+2])) == NULL)
        die("Cannot read second --seq2 file: %s", argv[argi+2]);
      argi += 3;
    }
    else if(!strcmp(argv[argi], "--seqi") || !strcmp(argv[argi],"-i"))
    {
      if(argi + 2 >= argc) cmd_print_usage("Missing arguments");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file: %s", argv[argi+1]);
      argi += 2;
    }
    else cmd_print_usage("Unexpected argument: %s", argv[argi]);
  }

  if(use_fq && use_fa)
    cmd_print_usage("Cannot print in both --fasta and --fastq");
  else if(!use_fa && !use_fq) use_fq = true;

  if(argi >= argc)
    cmd_print_usage("Please specify input graph files");

  int argend = argi;

  const size_t num_gfiles = (size_t)(argc - argend);
  size_t i;
  char **graph_paths = argv + argend;


  //
  // Open input graphs
  //
  GraphFileReader gfiles[num_gfiles];
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  graph_files_open(graph_paths, gfiles, num_gfiles,
                   &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Calculate memory use
  //
  size_t kmers_in_hash, graph_mem;

  kmers_in_hash = cmd_get_kmers_in_hash(args, 0, ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);
  cmd_check_mem_limit(args, graph_mem);

  //
  // Test output files
  //
  for(argi = 0; argi < argend; argi++)
  {
    if(!strcmp(argv[argi],"--seq") || !strcmp(argv[argi],"-1")) {
      check_outfile_exists(argv[argi+2], false, use_fq);
      argi += 2;
    }
    else if(!strcmp(argv[argi],"--seq2") || !strcmp(argv[argi],"-2")) {
      check_outfile_exists(argv[argi+3], true, use_fq);
      argi += 3;
    }
    else if(!strcmp(argv[argi],"--seqi") || !strcmp(argv[argi],"-i")) {
      check_outfile_exists(argv[argi+2], true, use_fq);
      argi += 2;
    }
  }

  //
  // Set up graph
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 0, kmers_in_hash);

  // Load graphs
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .must_exist_in_graph = false,
                              .empty_colours = true,
                              .boolean_covgs = false};

  for(i = 0; i < num_gfiles; i++) {
    gfiles[i].fltr.flatten = true;
    file_filter_update_intocol(&gfiles[i].fltr, 0);
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
  }

  status("Printing reads that do %stouch the graph\n", invert ? "not " : "");

  //
  // Filtering reads
  //
  read_t r1, r2;
  size_t total_reads_printed = 0;

  if(seq_read_alloc(&r1) == NULL || seq_read_alloc(&r2) == NULL)
    die("Out of memory");

  for(argi = 0; argi < argend; argi++)
  {
    char is_se = (!strcmp(argv[argi],"--seq")  || !strcmp(argv[argi],"-1"));
    char is_pe = (!strcmp(argv[argi],"--seq2") || !strcmp(argv[argi],"-2"));
    bool is_interleaved = (!strcmp(argv[argi],"--seqi")  || !strcmp(argv[argi],"-i"));

    if(is_se || is_pe || is_interleaved)
    {
      char *in1 = NULL, *in2 = NULL, *out = NULL;
      size_t init_reads, reads_loaded;

      AlignReadsData data = {&db_graph, &stats, in1, in2, NULL, NULL, 0,
                             use_fq ? seq_gzprint_fastq : seq_gzprint_fasta,
                             invert};

      if(is_se || is_interleaved) {
        in1 = argv[argi+1];
        out = argv[argi+2];
      }
      else if(is_pe) {
        in1 = argv[argi+1];
        in2 = argv[argi+2];
        out = argv[argi+3];
      }

      size_t pathlen = strlen(out);
      char path1[pathlen+8+1], path2[pathlen+8+1];

      memcpy(path1, out, pathlen);
      get_out_path(path1, pathlen, use_fq, is_pe ? 1 : 0);
      if((data.out1 = gzopen(path1, "w")) == NULL)
        die("Cannot write to: %s", path1);

      if(is_pe || is_interleaved) {
        memcpy(path2, out, pathlen);
        get_out_path(path2, pathlen, use_fq, 2);
        if((data.out2 = gzopen(path2, "w")) == NULL)
          die("Cannot write to: %s", path2);
      }

      init_reads = stats.num_se_reads + stats.num_pe_reads;

      if(is_pe) {
        status("reading: %s %s\n", in1, in2);
        status("writing: %s %s\n", path1, path2);
        seq_parse_pe_sf(seqfiles[sf], seqfiles[sf+1], 0, &r1, &r2,
                        filter_reads, &data);
        seq_close(seqfiles[sf]);
        seq_close(seqfiles[sf+1]);
        sf += 2;
      }
      else if(is_interleaved) {
        status("reading: %s (interleaved)\n", in1);
        status("writing: %s %s\n", path1, path2);
        seq_parse_interleaved_sf(seqfiles[sf], 0, &r1, &r2,
                                 filter_reads, &data);
        seq_close(seqfiles[sf]);
        sf++;
      }
      else {
        status("reading: %s\n", in1);
        status("writing: %s\n", path1);
        seq_parse_se_sf(seqfiles[sf], 0, &r1,
                        filter_reads, &data);
        seq_close(seqfiles[sf]);
        sf++;
      }

      gzclose(data.out1);
      if(is_pe || is_interleaved) gzclose(data.out2);

      total_reads_printed += data.num_of_reads_printed;
      reads_loaded = stats.num_se_reads + stats.num_pe_reads - init_reads;

      status("  Printed %zu / %zu inputs\n",
             data.num_of_reads_printed, reads_loaded);

      argi += ((is_se || is_interleaved) ? 2 : 3);
    }
  }

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  size_t total_reads = stats.num_se_reads + stats.num_pe_reads;

  status("Total printed %zu / %zu reads\n", total_reads_printed, total_reads);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
