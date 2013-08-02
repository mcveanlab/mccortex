#include "global.h"

#include <zlib.h>
#include "string_buffer.h"
#include "seq_file.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "seq_reader.h"
#include "binary_format.h"

static const char usage[] =
"usage: "CMD" reads [options] <in.ctx>[:cols] [in2.ctx ...]\n"
"  Filters reads based on which have a kmer in the graph. \n"
"  Options:\n"
"    -m <mem> Memory limit\n"
"    -h <mem> Number of kmers in hash table\n"
"\n"
"    --fasta   print output as gzipped FASTA\n"
"    --fastq   print output as gzipped FASTQ [default]\n"
"    --invert  print reads/read pairs with no kmer in graph\n"
"\n"
"    --se <in> <out>         Writes output to <out>.fq.gz\n"
"    --pe <in1> <in2> <out>  Writes output to <out>.1.fq.gz, <out>.2.fq.gz\n"
"\n"
"      Can specify --se/--pe multiple times.\n";

typedef struct {
  char *in1, *in2;
  gzFile out1, out2;
  size_t num_of_reads_printed;
  void (*print)(const read_t *r, gzFile gz, int linewrap);
  boolean invert;
} AlignReadsData;

static void get_out_path(char *path, size_t len, boolean use_fq, int pe_num)
{
  if(pe_num == 0) {
    memcpy(path+len, use_fq ? ".fq.gz" : ".fa.gz", 6);
    path[len+6] = '\0';
  }
  else
    sprintf(path+len, use_fq ? ".%i.fq.gz" : ".%i.fa.gz", pe_num);
}

static void check_outfile_exists(char *outbase, boolean is_pe, boolean use_fq)
{
  size_t len = strlen(outbase);
  char path[len+1+8]; // .1.fq.gz
  memcpy(path, outbase, len);
  path[len] = '\0';

  if(is_pe) {
    get_out_path(path, len, use_fq, 1);
    if(file_exists(path)) die("Output file already exists: %s", path);
    if(!test_file_writable(path)) print_usage(usage, "Cannot write: %s", path);
    get_out_path(path, len, use_fq, 2);
    if(file_exists(path)) die("Output file already exists: %s", path);
    if(!test_file_writable(path)) print_usage(usage, "Cannot write: %s", path);
  }
  else {
    get_out_path(path, len, use_fq, 0);
    if(file_exists(path)) die("Output file already exists: %s", path);
    if(!test_file_writable(path)) print_usage(usage, "Cannot write: %s", path);
  }
}

static hkey_t find_node(BinaryKmer bkmer, const dBGraph *db_graph)
{
  BinaryKmer tmpkey;
  db_node_get_key(bkmer, db_graph->kmer_size, tmpkey);
  return hash_table_find(&db_graph->ht, tmpkey);
}

static boolean read_touches_graph(const read_t *r, const dBGraph *db_graph,
                                  SeqLoadingStats *stats)
{
  boolean found = false;
  uint32_t kmer_size = db_graph->kmer_size, num_contigs = 0, kmers_loaded = 0;

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

      binary_kmer_from_str(r->seq.b + start, kmer_size, bkmer);
      kmers_loaded++;
      if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }

      for(i = start+kmer_size; i < end; i++)
      {
        nuc = binary_nuc_from_char(r->seq.b[i]);
        binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
        kmers_loaded++;
        if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }
      }
    }
  }

  // Update stats
  stats->total_bases_read += r->seq.end;
  stats->kmers_loaded += kmers_loaded;
  stats->unique_kmers += kmers_loaded - found;
  if(num_contigs > 0) stats->total_good_reads++;
  else stats->total_bad_reads++;

  return found;
}

void filter_reads(read_t *r1, read_t *r2,
                  int qoffset1, int qoffset2,
                  SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  (void)prefs; (void)stats;

  AlignReadsData *data = (AlignReadsData*)ptr;
  const dBGraph *db_graph = prefs->db_graph;

  boolean touches_graph = read_touches_graph(r1, db_graph, stats) ||
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
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 4) print_usage(usage, NULL);

  // Test other input args
  // Check filelists are readable
  // Check output is writable

  boolean use_fq = false, use_fa = false, invert = false;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcasecmp(argv[argi], "--fastq") == 0) use_fq = true;
    else if(strcasecmp(argv[argi], "--fasta") == 0) use_fa = true;
    else if(strcasecmp(argv[argi], "--invert") == 0) invert = true;
    else if(strcasecmp(argv[argi], "--se") == 0)
    {
      if(argi + 2 >= argc) print_usage(usage, "Missing arguments");
      char *in = argv[argi+1];
      if(!test_file_readable(in)) print_usage(usage, "Cannot read: %s", in);
      argi += 2;
    }
    else if(strcasecmp(argv[argi], "--pe") == 0)
    {
      if(argi + 4 >= argc) print_usage(usage, "Missing arguments");
      char *in1 = argv[argi+1], *in2 = argv[argi+2];
      if(!test_file_readable(in1)) print_usage(usage, "Cannot read: %s", in1);
      if(!test_file_readable(in2)) print_usage(usage, "Cannot read: %s", in2);
      argi += 3;
    }
    else print_usage(usage, "Unexpected argument: %s", argv[argi]);
  }

  if(use_fq && use_fa)
    print_usage(usage, "Cannot print in both --fasta and --fastq");
  else if(!use_fa && !use_fq) use_fq = true;

  int argend = argi;

  int i, num_binaries = argc - argend;
  char *binary_paths[num_binaries];

  if(num_binaries == 0)
    print_usage(usage, "Please specify input graph files");

  //
  // Probe binaries to get kmer-size
  //
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2, num_of_cols, max_col;
  uint64_t num_kmers, max_num_kmers = 0;

  for(i = 0; i < num_binaries; i++)
  {
    binary_paths[i] = argv[argend+i];

    if(!binary_probe(binary_paths[i], &is_binary, &kmer_size2,
                     &num_of_cols, &max_col, &num_kmers)) {
      print_usage(usage, "Cannot read binary file: %s", binary_paths[i]);
    } else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0) kmer_size = kmer_size2;
    else if(kmer_size != kmer_size2) {
      die("Graph kmer-sizes do not match [%u vs %u; %s; %s]\n",
          kmer_size, kmer_size2, binary_paths[i-1], binary_paths[i]);
    }

    max_num_kmers = MAX2(num_kmers, max_num_kmers);
  }

  //
  // Calculate memory use
  //
  size_t kmers_in_hash = cmd_get_kmers_in_hash(args, 0);
  size_t ideal_capacity = max_num_kmers*(1.0/IDEAL_OCCUPANCY);

  char num_kmers_str[100];
  ulong_to_str(max_num_kmers, num_kmers_str);

  if(kmers_in_hash < max_num_kmers) {
    print_usage(usage, "Not enough kmers in the hash, require: %s "
                       "(set bigger -h <kmers> or -m <mem>)", num_kmers_str);
  }
  else if(kmers_in_hash < ideal_capacity)
    warn("Low memory for binary size (require: %s)", num_kmers_str);

  //
  // Test output files
  //
  for(argi = 0; argi < argend; argi++)
  {
    if(strcasecmp(argv[argi],"--se") == 0) {
      check_outfile_exists(argv[argi+2], false, use_fq);
      argi += 2;
    }
    else if(strcasecmp(argv[argi],"--pe") == 0) {
      check_outfile_exists(argv[argi+3], true, use_fq);
      argi += 3;
    }
  }

  //
  // Set up graph
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, num_of_cols, kmers_in_hash);

  // Load binaries
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                           .boolean_covgs = false,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = true,
                           .update_ginfo = false,
                           .db_graph = &db_graph};

  for(i = 0; i < num_binaries; i++)
    binary_load(binary_paths[i], &prefs, stats, NULL);

  prefs.load_seq = true;
  prefs.load_binaries = false;

  if(invert) message("Printing reads that do not touch the graph\n");
  else message("Printing reads that touch the graph\n");

  //
  // Filtering reads
  //
  read_t r1, r2;
  size_t total_reads_printed = 0;

  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

  for(argi = 0; argi < argend; argi++)
  {
    char is_se = (strcasecmp(argv[argi], "--se") == 0);
    char is_pe = (strcasecmp(argv[argi], "--pe") == 0);

    if(is_se || is_pe)
    {
      char *in1 = NULL, *in2 = NULL, *out = NULL;
      size_t init_reads, reads_loaded;

      AlignReadsData data = {in1, in2, NULL, NULL, 0,
                             use_fq ? seq_gzprint_fastq : seq_gzprint_fasta,
                             invert};

      if(is_se) {
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

      if(is_pe) {
        memcpy(path2, out, pathlen);
        get_out_path(path2, pathlen, use_fq, 2);
        if((data.out2 = gzopen(path2, "w")) == NULL)
          die("Cannot write to: %s", path2);
        message("writing: %s\n", path2);
      }

      init_reads = stats->total_good_reads + stats->total_bad_reads +
                   stats->total_dup_reads;

      if(is_pe) {
        message("reading: %s %s\n", in1, in2);
        message("writing: %s %s\n", path1, path2);
        seq_parse_pe(in1, in2, &r1, &r2, &prefs, stats, filter_reads, &data);
      } else {
        message("reading: %s\n", in1);
        message("writing: %s\n", path1);
        seq_parse_se(in1, &r1, &r2, &prefs, stats, filter_reads, &data);
      }

      gzclose(data.out1);
      if(is_pe) gzclose(data.out2);

      total_reads_printed += data.num_of_reads_printed;
      reads_loaded = stats->total_good_reads + stats->total_bad_reads +
                     stats->total_dup_reads - init_reads;

      message("  Printed %zu / %zu inputs\n",
              data.num_of_reads_printed, reads_loaded);

      if(is_se) argi += 2;
      else argi += 3;
    }
  }

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  size_t total_reads = stats->total_good_reads + stats->total_bad_reads +
                       stats->total_dup_reads;

  message("Total printed %zu / %zu reads\n", total_reads_printed, total_reads);
  message("Done\n");

  // free
  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  exit(EXIT_SUCCESS);
}
