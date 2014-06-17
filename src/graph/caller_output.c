#include "global.h"
#include "caller_output.h"
#include "file_util.h"
#include "graph_info.h"
#include "cmd.h" // define cmd_get_cmdline() and cmd_get_cwd()

void caller_gzprint_ginfo(const GraphInfo *ginfo, size_t ncols, gzFile gzout)
{
  size_t col;
  StrBuf *sample_name = strbuf_new();

  // Print sample names
  for(col = 0; col < ncols; col++, ginfo++)
  {
    const ErrorCleaning *ec = &ginfo->cleaning;

    // Find and replace double quotes with single quotes
    char *sname = ginfo->sample_name.buff;

    if(strcmp(sname, "undefined") == 0 || strchr(sname, '\t') != NULL ||
       strchr(sname, ' ') != NULL || strchr(sname, '\r') != NULL ||
       strchr(sname, '\n') != NULL)
    {
      strbuf_reset(sample_name);
      strbuf_sprintf(sample_name, "sample%zu", col);
    }
    else {
      strbuf_set(sample_name, ginfo->sample_name.buff);
    }

    gzprintf(gzout, "##colour=<ID=%s,name=\"%s\",colour=%i,"
                    "meanreadlen=%zu,totalseqloaded=%zu,"
                    "seqerror=%Lf,tipclipped=%s,removelowcovgsupernodes=%u,"
                    "removelowcovgkmer=%u,cleanedagainstgraph=%s>\n",
             sample_name->buff, ginfo->sample_name.buff, col,
             (size_t)ginfo->mean_read_length, (size_t)ginfo->total_sequence,
             ginfo->seq_err,
             ec->cleaned_tips ? "yes" : "no", ec->clean_snodes_thresh,
             ec->clean_kmers_thresh,
             ec->intersection_name.buff);
  }

  strbuf_free(sample_name);
}

// Print header with absolute path to a file
void caller_gzprint_path_hdr(gzFile gzout, const char *name, const char *path)
{
  char absolute_path[PATH_MAX + 1];

  if(realpath(path, absolute_path) == NULL)
    warn("Cannot get absolute path: %s\n", path);
  else
    path = absolute_path;

  gzprintf(gzout, "##%s=%s\n", name, path);
}

void caller_gzprint_header(gzFile gzout, const char* out_file,
                           const char *format_str, const dBGraph *db_graph)
{
  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  gzprintf(gzout, "##fileFormat=%s\n", format_str);
  gzprintf(gzout, "##ctxCmd=\"%s\"\n", cmd_get_cmdline());
  gzprintf(gzout, "##ctxCwd=%s\n", cmd_get_cwd());
  gzprintf(gzout, "##ctxDate=%s\n", datestr);
  gzprintf(gzout, "##ctxVersion=<version=%s,MAXK=%i>\n",
           CTX_VERSION, MAX_KMER_SIZE);
  gzprintf(gzout, "##ctxKmerSize=%u\n", db_graph->kmer_size);
  caller_gzprint_path_hdr(gzout, "outPath", out_file);
  gzprintf(gzout, "##ctxNumColoursUsedInCalling=%i\n", db_graph->num_of_cols);

  caller_gzprint_ginfo(db_graph->ginfo, db_graph->num_of_cols, gzout);
}
