#include "global.h"

// third party libraries
#include <string_buffer.h>
#include "seq_file.h"

// cortex_var headers
#include "file_reader.h"
#include "file_util.h"
#include "util.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "graph_info.h"
#include "seq_reader.h"
#include "binary_format.h"

//
// Create/free/sum SeqLoadingStats
//

SeqLoadingStats* seq_loading_stats_create(unsigned long readlen_arrsize)
{
  SeqLoadingStats *stats = (SeqLoadingStats*)malloc(sizeof(SeqLoadingStats));
  if(stats == NULL) die("Out of memory");

  unsigned long *readlen_arr = NULL;

  if(readlen_arrsize > 0) {
    readlen_arr = (unsigned long*)calloc(readlen_arrsize, sizeof(unsigned long));
    if(readlen_arr == NULL) die("Out of memory");
  }

  stats->readlen_count_array = readlen_arr;
  stats->readlen_count_array_size = readlen_arrsize;

  stats->se_colourlists_loaded = stats->pe_colourlist_pairs_loaded = 0;
  stats->se_filelists_loaded = stats->pe_filelist_pairs_loaded = 0;
  stats->num_files_loaded = stats->binaries_loaded = 0;
  stats->num_se_reads = stats->num_pe_reads = 0;
  stats->total_good_reads = stats->total_bad_reads = stats->total_dup_reads = 0;
  stats->total_bases_read = stats->total_bases_loaded = 0;
  // unique_kmers is the number of new kmers added to the graph
  stats->kmers_loaded = stats->unique_kmers = 0;

  // Used for binaries and colourlists
  stats->num_of_colours_loaded = 0;

  return stats;
}

void seq_loading_stats_sum(SeqLoadingStats* dst, SeqLoadingStats* src)
{
  dst->se_colourlists_loaded += src->se_colourlists_loaded;
  dst->pe_colourlist_pairs_loaded += src->pe_colourlist_pairs_loaded;
  dst->se_filelists_loaded += src->se_filelists_loaded;
  dst->pe_filelist_pairs_loaded += src->pe_filelist_pairs_loaded;
  dst->num_files_loaded += src->num_files_loaded;
  dst->binaries_loaded += src->binaries_loaded;

  dst->num_se_reads += src->num_se_reads;
  dst->num_pe_reads += src->num_pe_reads;
  dst->total_good_reads += src->total_good_reads;
  dst->total_bad_reads += src->total_bad_reads;
  dst->total_dup_reads += src->total_dup_reads;
  dst->total_bases_read += src->total_bases_read;
  dst->total_bases_loaded += src->total_bases_loaded;
  dst->kmers_loaded += src->kmers_loaded;
  dst->unique_kmers += src->unique_kmers;

  // Used for binaries and colourlists
  dst->num_of_colours_loaded += src->num_of_colours_loaded;

  uint64_t i;
  uint64_t limit = MIN2(dst->readlen_count_array_size,
                        src->readlen_count_array_size);

  for(i = 0; i < limit; i++)
  {
    dst->readlen_count_array[i] += src->readlen_count_array[i];
  }

  if(dst->readlen_count_array_size < src->readlen_count_array_size)
  {
    uint64_t sum = 0;

    for(i = dst->readlen_count_array_size; i < src->readlen_count_array_size; i++)
    {
      sum += src->readlen_count_array[i];
    }

    dst->readlen_count_array[dst->readlen_count_array_size-1] += sum;
  }
}

void seq_loading_stats_free(SeqLoadingStats* stats)
{
  if(stats->readlen_count_array != NULL) free(stats->readlen_count_array);
  free(stats);
}

//
// Parse filelists and colourlists
//

// list_path1 and/or list_path2 may be null
void parse_filelists(const char *list_path1, const char *list_path2,
                     uint8_t are_colour_lists,
                     SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                     void (*read_func)(read_t *r1, read_t *r2,
                                       int qoffset1, int qoffset2,
                                       SeqLoadingPrefs *prefs,
                                       SeqLoadingStats *stats, void *ptr),
                     void *reader_ptr)
{
  if(list_path1 == NULL && list_path2 == NULL) return;

  boolean is_pe = (list_path1 != NULL && list_path2 != NULL);

  boolean has_sample_names = false;
  StrBuf **sample_names = NULL;
  uint32_t num_files1 = 0, num_files2 = 0, num_files = 0, i;
  size_t initial_kmers = prefs->db_graph->ht.unique_kmers;

  if(is_pe) {
    message("Parse paired-end %s: %s; %s\n",
            are_colour_lists ? "sequence colour lists" : "files",
            list_path1, list_path2);
  } else {
    message("Parse single-ended %s: %s\n",
            are_colour_lists ? "sequence colour list" : "files",
            list_path1 != NULL ? list_path1 : list_path2);
  }

  if(list_path1 != NULL)
  {
    num_files1 = load_paths_from_filelist(list_path1, NULL, are_colour_lists,
                                          NULL, &has_sample_names);
    num_files = num_files1;
  }

  if(list_path2 != NULL)
  {
    num_files2 = load_paths_from_filelist(list_path2, NULL, are_colour_lists,
                                          NULL, &has_sample_names);
    num_files = num_files2;
  }

  if(is_pe && num_files1 != num_files2)
  {
    die("list files have different number of entries\nf1: %s [%u]\nf2:%s [%u]",
        list_path1, num_files1, list_path2, num_files2);
  }

  char **filelist1 = (char**)malloc(num_files * sizeof(char*));
  char **filelist2 = (char**)malloc(num_files * sizeof(char*));

  // Set all entries to NULL
  for(i = 0; i < num_files; i++) { filelist1[i] = filelist2[i] = NULL; }

  // Check if we should load sample names
  if(are_colour_lists && prefs->update_ginfo) {
    sample_names = prefs->db_graph->ginfo.sample_names + prefs->into_colour;
  }

  if(list_path1 != NULL) {
    load_paths_from_filelist(list_path1, filelist1, are_colour_lists,
                             sample_names, &has_sample_names);
  }

  if(list_path2 != NULL) {
    load_paths_from_filelist(list_path2, filelist2, are_colour_lists,
                             sample_names, &has_sample_names);
  }

  SeqLoadingStats *new_stats
    = seq_loading_stats_create(stats->readlen_count_array_size);

  boolean is_cortex_binary1 = false, is_cortex_binary2 = false;
  uint32_t binary_kmer_size1, colours_in_binary1;
  uint32_t binary_kmer_size2, colours_in_binary2;
  uint64_t kmers_in_binary1, kmers_in_binary2;

  read_t *r1 = NULL, *r2 = NULL;
  
  if(!are_colour_lists) {
    r1 = seq_read_new();
    r2 = seq_read_new();
  }

  for(i = 0; i < num_files; i++)
  {
    char *p1 = filelist1[i], *p2 = filelist2[i];
    is_cortex_binary1 = is_cortex_binary2 = false;

    if(p1 == NULL && p2 == NULL) continue;

    if(p1 != NULL && !binary_probe(p1, &is_cortex_binary1,
                                   &binary_kmer_size1, &colours_in_binary1,
                                   &kmers_in_binary1))
    {
      die("Cannot read file: %s \n[listed in: %s]", p1, list_path1);
    }

    if(p2 != NULL && !binary_probe(p2, &is_cortex_binary2,
                                   &binary_kmer_size2, &colours_in_binary2,
                                   &kmers_in_binary2))
    {
      die("Cannot read file: %s \n[listed in: %s]", p2, list_path2);
    }

    if(!prefs->load_binaries)
    {
      if(is_cortex_binary1) {
        warn("Skipping binary [binary: %s; colourlist: %s]", p1, list_path1);
        p1 = NULL; is_cortex_binary1 = false;
      }
      if(is_cortex_binary2) {
        warn("Skipping binary [binary: %s; colourlist: %s]", p2, list_path2);
        p2 = NULL; is_cortex_binary2 = false;
      }
    }

    // Can only load single colour binaries in filelists
    if(!are_colour_lists && is_cortex_binary1 && colours_in_binary1 != 1)
    {
      warn("Cannot load multicolour binary along with sequence -- "
           "skipping [binary: %s; cols: %i; filelist: %s]",
           filelist1[i], colours_in_binary1, list_path1);
      p1 = NULL; is_cortex_binary1 = false;
    }
    if(!are_colour_lists && is_cortex_binary2 && colours_in_binary2 != 1)
    {
      warn("Cannot load multicolour binary along with sequence -- "
           "skipping [binary: %s; cols: %i; filelist: %s]",
           filelist2[i], colours_in_binary2, list_path2);
      p2 = NULL; is_cortex_binary2 = false;
    }

    // This is only used if we are reading a colour list
    uint32_t cols_read_in = 0;

    // Load binaries if passed
    if(is_cortex_binary1 && is_cortex_binary2)
    {
      if(colours_in_binary1 != colours_in_binary2) {
        warn("Number of colours in binaries doesn't match [%u != %u; %s; %s]",
             colours_in_binary1, colours_in_binary2, p1, p2);
      }

      // binary_load updates stats->num_of_colours_loaded
      binary_load(p1, prefs->db_graph,
                  prefs->into_colour, prefs->must_exist_in_colour,
                  prefs->empty_colours, prefs->load_as_union,
                  new_stats);
      // (undo auto-update of stats)
      new_stats->num_of_colours_loaded -= colours_in_binary1;

      binary_load(p2, prefs->db_graph,
                  prefs->into_colour, prefs->must_exist_in_colour,
                  prefs->empty_colours, prefs->load_as_union,
                  new_stats);
      // (undo auto-update of stats)
      new_stats->num_of_colours_loaded -= colours_in_binary2;

      cols_read_in = MAX2(colours_in_binary1, colours_in_binary2);

      // Mark null so we know we're dealt with this files
      p1 = p2 = NULL;
    }
    else if(is_cortex_binary1) {
      binary_load(p1, prefs->db_graph,
                  prefs->into_colour, prefs->must_exist_in_colour,
                  prefs->empty_colours, prefs->load_as_union,
                  new_stats);
      cols_read_in = colours_in_binary1;
      p1 = NULL;
    }
    else if(is_cortex_binary2) {
      binary_load(p2, prefs->db_graph,
                  prefs->into_colour, prefs->must_exist_in_colour,
                  prefs->empty_colours, prefs->load_as_union,
                  new_stats);
      cols_read_in = colours_in_binary2;
      p2 = NULL;
    }

    // Load non-binaries
    if(are_colour_lists)
    {
      if(p1 != NULL || p2 != NULL)
      {
        parse_filelists(p1, p2, READ_FALIST, prefs, new_stats,
                        read_func, reader_ptr);

        cols_read_in = MAX2(cols_read_in, 1);

        if(is_pe) new_stats->pe_colourlist_pairs_loaded++;
        else new_stats->se_colourlists_loaded++;
      }

      prefs->into_colour += cols_read_in;
      new_stats->num_of_colours_loaded += cols_read_in;
    }
    else
    {
      if(prefs->load_seq)
      {
        if(p1 != NULL && p2 != NULL)
        {
          seq_parse_pe(p1, p2, r1, r2, prefs, new_stats,
                       read_func, reader_ptr);
        }
        else if(p1 != NULL)
        {
          seq_parse_se(p1, r1, r2, prefs, new_stats,
                       read_func, reader_ptr);
        }
        else if(p2 != NULL)
        {
          seq_parse_se(p2, r1, r2, prefs, new_stats,
                       read_func, reader_ptr);
        }
      }
      else
      {
        if(p1 != NULL) {
          warn("Skipping sequence file [%s; list: %s]", p1, list_path1);
        }
        if(p2 != NULL) {
          warn("Skipping sequence file [%s; list: %s]", p2, list_path2);
        }
      }
    }
  }

  if(!are_colour_lists)
  {
    seq_read_free(r1);
    seq_read_free(r2);

    if(prefs->update_ginfo)
    {
      // Update the graph info object
      uint32_t mean_contig_length
        = calculate_mean_ulong(new_stats->readlen_count_array,
                               new_stats->readlen_count_array_size);

      graph_info_update_mean_readlen_and_total_seq(&prefs->db_graph->ginfo,
                                                   prefs->into_colour,
                                                   mean_contig_length,
                                                   new_stats->total_bases_loaded);
    }

    if(is_pe) new_stats->pe_filelist_pairs_loaded++;
    else new_stats->se_filelists_loaded++;

    new_stats->unique_kmers += prefs->db_graph->ht.unique_kmers - initial_kmers;

    // Print stats for this set of files
    message("\n");
    message("Number of sequence files loaded: %lu\n", new_stats->num_files_loaded);

    message("  number of SE reads: %lu\n", new_stats->num_se_reads);
    message("  number of PE reads: %lu\n", new_stats->num_pe_reads);
    message("  number of reads (good,bad,dupe): (%lu, %lu, %lu)\n",
            new_stats->total_good_reads, new_stats->total_bad_reads,
            new_stats->total_dup_reads);
    message("  sequence parsed: %lu\n", new_stats->total_bases_read);
    message("  sequence passed filters: %lu\n", new_stats->total_bases_loaded);
    message("  kmers parsed (of which novel): %lu (%lu)\n",
            new_stats->kmers_loaded, new_stats->unique_kmers);
    message("  total kmers in graph: %zu\n\n",
            (size_t)prefs->db_graph->ht.unique_kmers);
  }

  for(i = 0; i < num_files; i++) {
    if(filelist1[i] != NULL) free(filelist1[i]);
    if(filelist2[i] != NULL) free(filelist2[i]);
  }

  free(filelist1);
  free(filelist2);

  // Update cumulative stats
  seq_loading_stats_sum(stats, new_stats);
  seq_loading_stats_free(new_stats);
}

// Given a 'filelist' file, check all files pointed to exist.
// If path_array is not NULL, populate it with the paths.
// Exits with error if a file doesn't exist or isn't readable.
// Return the number of files pointed to.
uint32_t load_paths_from_filelist(const char *filelist_path, char **path_array,
                                  boolean sample_names_permitted,
                                  StrBuf **sample_names,
                                  boolean *has_sample_names)
{
  // Get absolute path
  char absolute_path[PATH_MAX + 1];
  char *filelist_abs_path = realpath(filelist_path, absolute_path);

  if(filelist_abs_path == NULL)
  {
    die("Cannot get absolute path to filelist: %s\n", filelist_path);
  }

  FILE *filelist_handle = fopen(filelist_abs_path, "r");

  if(filelist_handle == NULL)
  {
    die("Cannot open filelist: %s\n", filelist_abs_path);
  }

  // Get directory path
  StrBuf *dir = strbuf_new();
  file_reader_get_strbuf_of_dir_path(filelist_abs_path, dir);

  // Read filelist
  uint32_t lineno = 0;
  boolean seen_sample_names = false;

  StrBuf *line = strbuf_new();

  for(lineno = 0; strbuf_reset_readline(line, filelist_handle); lineno++)
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/')
      {
        strbuf_insert(line, 0, dir->buff, strbuf_len(dir));
      }

      // Replace the first '\t' with '\0'
      strtok(line->buff, "\t");

      // Get absolute paths
      char *path_ptr = realpath(line->buff, absolute_path);

      if(path_ptr == NULL) die("Cannot find file: %s\n", line->buff);

      if(path_array != NULL)
      {
        path_array[lineno] = strdup(line->buff);
        if(path_array[lineno] == NULL) {
          die("Out of memory loading filelist: %s\n", filelist_path);
        }
      }

      // Get sample name
      char *sample_name = strtok(NULL, "\t");

      if(lineno == 0)
        seen_sample_names = (sample_name != NULL);
      else if((sample_name != NULL) != seen_sample_names)
        die("Only some lines of list file have sample names: %s", filelist_path);

      if(seen_sample_names) {
        if(!sample_names_permitted)
          die("sample names not permitted: %s", filelist_path);

        strbuf_set(sample_names[lineno], sample_name);
      }
      else if(sample_names != NULL) {
        strbuf_set(sample_names[lineno], "undefined");
      }
    }
    else if(path_array != NULL) path_array[lineno] = NULL;
  }

  strbuf_free(line);
  strbuf_free(dir);

  fclose(filelist_handle);

  if(sample_names_permitted) *has_sample_names = seen_sample_names;

  return lineno;
}

// filename is a list of files, one for each colour (with optional second column
// of sample-ids). Check they all exists, there are not too many, and that each
// of them contains a list of valid binaries.

// Returns number for files in file passed
// Check that ctxlists contain only valid ctx files
// Check that colourlists don't contain any ctx files
uint32_t check_colour_or_ctx_list(const char *list_path, char is_colourlist,
                                  boolean binaries_allowed, boolean seq_allowed,
                                  uint32_t kmer_size)
{
  boolean has_sample_names;
  uint32_t num_files_in_list, i;

  num_files_in_list = load_paths_from_filelist(list_path, NULL, is_colourlist,
                                               NULL, &has_sample_names);

  // Empty file list is valid
  if(num_files_in_list == 0) return 0;

  if(is_colourlist && num_files_in_list > NUM_OF_COLOURS)
  {
    die("Too many colours in colourlist [cols: %i; compiled for: %i; list: %s]",
        num_files_in_list, NUM_OF_COLOURS, list_path);
  }

  char **file_paths = malloc(sizeof(char*) * num_files_in_list);
  assert(file_paths != NULL);

  load_paths_from_filelist(list_path, file_paths, is_colourlist,
                           NULL, &has_sample_names);

  for(i = 0; i < num_files_in_list; i++)
  {
    if(file_paths[i] != NULL)
    {
      boolean is_cortex_binary;
      uint32_t binary_kmer_size, colours_in_binary;
      uint64_t kmers_in_binary;

      if(!binary_probe(file_paths[i], &is_cortex_binary,
                       &binary_kmer_size, &colours_in_binary,
                       &kmers_in_binary))
      {
        die("Cannot read file: %s\nlisted in file: %s", file_paths[i], list_path);
      }

      if(!binaries_allowed && is_cortex_binary)
      {
        die("Didn't expect binary file in fasta/fastq/bam list\n"
            "  Found binary file: %s\n"
            "  in list: %s\n", file_paths[i], list_path);
      }
      else if(!seq_allowed && !is_cortex_binary)
      {
        die("Didn't expect seq file in ctx list\n"
            "  Found binary file: %s\n"
            "  in list: %s\n", file_paths[i], list_path);
      }
      else if(is_cortex_binary && binary_kmer_size != kmer_size)
      {
        // Kmer size mismatch
        die("Binary file has different kmer size [file: %s; kmer: %i]",
            list_path, kmer_size);
      }

      if(is_colourlist)
      {
        check_colour_or_ctx_list(file_paths[i], READ_FALIST,
                                 binaries_allowed, seq_allowed,
                                 kmer_size);
      }
    }
  }

  // Cleanup
  for(i = 0; i < num_files_in_list; i++) free(file_paths[i]);
  free(file_paths);

  return num_files_in_list;
}


// takes a filename (a colour_list) this file contains a list of filenames
// (one per colour), each of these represents an individual (and contains a list
// of single-colour binaries for that individual). Also takes two colour numbers.
// clean_colour is the clean colour, and into_colour is the colour ALL of these
// individuals are loaded into, thus:
// First take person 0's list of binaries and load them all into colour into_colour,
// BUT only load nodes that are already in the hash, and only load those edges
// that are in the clean colour. Then dump a single-colour binary of colour
// into_colour, with filename = colour name PLUS a suffix added on the end.
void dump_successive_cleaned_binaries(const char *filename, uint32_t into_colour,
                                      uint32_t clean_colour, const char *suffix,
                                      dBGraph *db_graph,
                                      const char *cleaned_against_name)
{
  assert(into_colour != clean_colour);

  // Get absolute path
  char absolute_path[PATH_MAX + 1];
  char *filename_abs_path = realpath(filename, absolute_path);

  if(filename_abs_path == NULL) {
    die("Cannot get absolute path to .colours file: %s\n", filename);
  }

  FILE *fp = fopen(filename, "r");
  if(fp == NULL) die("Cannot open file: %s\n", filename);

  // Set up loading settings
  SeqLoadingPrefs seq_loading_prefs
    = {into_colour, // load into colour
       false, 0, 0, 0, false, false,
       true, // load binaries
       clean_colour, // must exist in this colour
       true, // empty colours
       false, // don't load as union
       true, // update ginfo
       db_graph};

  // Get directory path
  StrBuf *dir = strbuf_new();
  StrBuf *line = strbuf_new();
  StrBuf *output_path_str = strbuf_new();

  file_reader_get_strbuf_of_dir_path(filename_abs_path, dir);

  GraphInfo *ginfo = &db_graph->ginfo;

  while(strbuf_reset_readline(line, fp))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/') {
        strbuf_insert(line, 0, dir->buff, strbuf_len(dir));
      }

      // Get absolute paths
      char *path_ptr = realpath(line->buff, absolute_path);
      if(path_ptr == NULL) die("Cannot find ctxlist: %s\n", line->buff);

      // Reset ginfo for colour we're loading into
      graph_info_reset_one_colour(ginfo, into_colour);

      SeqLoadingStats *loading_stats = seq_loading_stats_create(0);

      parse_filelists(path_ptr, NULL, READ_FALIST,
                      &seq_loading_prefs, loading_stats,
                      &seq_load_into_db_graph, NULL);

      seq_loading_stats_free(loading_stats);

      strbuf_reset(output_path_str);
      strbuf_append_buff(output_path_str, line);
      if(suffix != NULL) {
        strbuf_append_char(output_path_str, '_');
        strbuf_append_str(output_path_str, suffix);
      }
      strbuf_append_str(output_path_str, ".ctx");

      // Set ginfo for cleaned against graph
      StrBuf *name = ginfo->cleaning[into_colour].cleaned_against_graph_name;
      strbuf_set(name, cleaned_against_name);
      ginfo->cleaning[into_colour].cleaned_against_another_graph = true;

      binary_dump_graph(output_path_str->buff, db_graph,
                        CURR_CTX_VERSION, NULL, into_colour,
                        1, ginfo->num_of_shades_loaded);

      // reset that colour in the graph
      db_graph_wipe_colour(db_graph, into_colour);
    }
  }

  strbuf_free(output_path_str);
  strbuf_free(line);
  strbuf_free(dir);

  fclose(fp);
}

