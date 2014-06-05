#include "global.h"
#include "graph_info.h"

static void error_cleaning_init(ErrorCleaning *ec)
{
  ec->cleaned_tips = ec->cleaned_snodes = ec->cleaned_kmers = false;
  ec->clean_snodes_thresh = ec->clean_kmers_thresh = 0;
  ec->is_graph_intersection = false;
  strbuf_set(&ec->intersection_name, "undefined");
}

static void error_cleaning_alloc(ErrorCleaning *ec)
{
  strbuf_alloc(&ec->intersection_name, 256);
  error_cleaning_init(ec);
}

static void error_cleaning_dealloc(ErrorCleaning *ec)
{
  strbuf_dealloc(&ec->intersection_name);
}

static void error_cleaning_cpy(ErrorCleaning *dst, const ErrorCleaning *src)
{
  dst->cleaned_tips = src->cleaned_tips;
  dst->cleaned_snodes = src->cleaned_snodes;
  dst->cleaned_kmers = src->cleaned_kmers;
  dst->clean_snodes_thresh = src->clean_snodes_thresh;
  dst->clean_kmers_thresh = src->clean_kmers_thresh;
  dst->is_graph_intersection = src->is_graph_intersection;
  strbuf_set(&dst->intersection_name, src->intersection_name.buff);
}

static void error_cleaning_merge(ErrorCleaning *dst, const ErrorCleaning *src,
                                 uint64_t dst_totalseq, uint64_t src_totalseq)
{
  dst->cleaned_tips |= src->cleaned_tips;
  dst->cleaned_snodes |= src->cleaned_snodes;
  dst->cleaned_kmers |= src->cleaned_kmers;

  if(src->clean_snodes_thresh > 0 &&
     (dst->clean_snodes_thresh == 0 ||
      src->clean_snodes_thresh < dst->clean_snodes_thresh))
  {
    dst->clean_snodes_thresh = src->clean_snodes_thresh;
  }

  if(src->clean_kmers_thresh > 0 &&
     (dst->clean_kmers_thresh == 0 ||
      src->clean_kmers_thresh < dst->clean_kmers_thresh))
  {
    dst->clean_kmers_thresh = src->clean_kmers_thresh;
  }

  if(src->is_graph_intersection)
    graph_info_append_intersect(dst, src->intersection_name.buff);

  dst->is_graph_intersection |= src->is_graph_intersection;
}

void graph_info_init(GraphInfo *ginfo)
{
  strbuf_set(&ginfo->sample_name, "undefined");
  ginfo->total_sequence = 0;
  ginfo->mean_read_length = 0;
  ginfo->seq_err = 0.01;
  error_cleaning_init(&ginfo->cleaning);
}

void graph_info_alloc(GraphInfo *ginfo)
{
  strbuf_alloc(&ginfo->sample_name, 256);
  error_cleaning_alloc(&ginfo->cleaning);
  graph_info_init(ginfo);
}

void graph_info_dealloc(GraphInfo *ginfo)
{
  strbuf_dealloc(&ginfo->sample_name);
  error_cleaning_dealloc(&ginfo->cleaning);
}

void graph_info_make_intersect(const GraphInfo *ginfo, StrBuf *intersect_name)
{
  if(intersect_name->len > 0) strbuf_append_char(intersect_name, ',');
  strbuf_append_str(intersect_name, ginfo->sample_name.buff);

  if(ginfo->cleaning.is_graph_intersection) {
    strbuf_append_char(intersect_name, ',');
    strbuf_append_str(intersect_name, ginfo->cleaning.intersection_name.buff);
  }
}

void graph_info_append_intersect(ErrorCleaning *cleaning, const char *intersect_name)
{
  if(!cleaning->is_graph_intersection)
  {
    strbuf_set(&cleaning->intersection_name, intersect_name);
  }
  else
  {
    strbuf_append_char(&cleaning->intersection_name, ',');
    strbuf_append_str(&cleaning->intersection_name, intersect_name);
  }
  cleaning->is_graph_intersection = true;
}

void graph_info_cpy(GraphInfo *dst, const GraphInfo *src)
{
  dst->mean_read_length = src->mean_read_length;
  dst->total_sequence = src->total_sequence;
  dst->seq_err = src->seq_err;
  strbuf_set(&dst->sample_name, src->sample_name.buff);
  error_cleaning_cpy(&dst->cleaning, &src->cleaning);
}

static void graph_info_update_contigs(GraphInfo *ginfo,
                                      uint64_t added_seq, uint64_t num_contigs)
{
  if(!added_seq && !num_contigs) return;

  size_t ginfo_num_contigs = 0;

  if(ginfo->total_sequence && ginfo->mean_read_length)
     ginfo_num_contigs = ((double)ginfo->total_sequence/ginfo->mean_read_length)+0.5;

 if(ginfo_num_contigs + num_contigs > 0) {
    ginfo->mean_read_length
      = (uint32_t)(((ginfo->total_sequence + added_seq) /
                    (ginfo_num_contigs + num_contigs)) + 0.5);
  }

  ginfo->total_sequence += added_seq;
}

void graph_info_merge(GraphInfo *dst, const GraphInfo *src)
{
  // Update sample name
  if(strcmp(src->sample_name.buff,"undefined") != 0) {
    if(strcmp(dst->sample_name.buff,"undefined") == 0) {
      strbuf_set(&dst->sample_name, src->sample_name.buff);
    } else {
      strbuf_append_char(&dst->sample_name, ',');
      strbuf_append_str(&dst->sample_name, src->sample_name.buff);
    }
  }

  uint64_t total_sequence = dst->total_sequence + src->total_sequence;

  if(total_sequence > 0)
  {
    // Average error rates
    dst->seq_err
      = (dst->seq_err * dst->total_sequence +
         src->seq_err * src->total_sequence) /
        total_sequence;

    // Update mean read length
    size_t src_num_contigs = 0;

    if(src->total_sequence && src->mean_read_length)
       src_num_contigs = ((double)src->total_sequence/src->mean_read_length)+0.5;

    graph_info_update_contigs(dst, src->total_sequence, src_num_contigs);
 }

  // Update error cleaning
  error_cleaning_merge(&dst->cleaning, &src->cleaning,
                       dst->total_sequence, src->total_sequence);

  dst->total_sequence = total_sequence;
}

void graph_info_update_stats(GraphInfo *ginfo, const LoadingStats *stats)
{
  graph_info_update_contigs(ginfo, stats->total_bases_loaded, stats->contigs_loaded);
}
