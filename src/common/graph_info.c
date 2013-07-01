#include "global.h"
#include "graph_info.h"

void error_cleaning_init(ErrorCleaning *ec)
{
  ec->tip_clipping = ec->remv_low_cov_sups = ec->remv_low_cov_nodes = false;
  ec->remv_low_cov_sups_thresh = ec->remv_low_cov_nodes_thresh = 0;
  ec->cleaned_against_another_graph = false;
  strbuf_set(&ec->cleaned_against_graph_name, "undefined");
}

void error_cleaning_alloc(ErrorCleaning *ec)
{
  strbuf_alloc(&ec->cleaned_against_graph_name, 256);
  error_cleaning_init(ec);
}

void error_cleaning_dealloc(ErrorCleaning *ec)
{
  strbuf_dealloc(&ec->cleaned_against_graph_name);
}

void error_cleaning_overwrite(ErrorCleaning *tgt, const ErrorCleaning *src)
{
  tgt->tip_clipping |= src->tip_clipping;
  tgt->remv_low_cov_sups |= src->remv_low_cov_sups;
  tgt->remv_low_cov_nodes |= src->remv_low_cov_nodes;
  tgt->remv_low_cov_sups_thresh = MAX2(tgt->remv_low_cov_sups_thresh,
                                       src->remv_low_cov_sups_thresh);
  tgt->remv_low_cov_nodes_thresh = MAX2(tgt->remv_low_cov_nodes_thresh,
                                        src->remv_low_cov_nodes_thresh);

  if(src->cleaned_against_another_graph &&
     src->cleaned_against_graph_name.len > 0 &&
     strcmp(src->cleaned_against_graph_name.buff,"undefined") != 0)
  {
    strbuf_set(&tgt->cleaned_against_graph_name,
               src->cleaned_against_graph_name.buff);
  }
  tgt->cleaned_against_another_graph |= src->cleaned_against_another_graph;
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

void graph_info_update_seq_stats(GraphInfo *ginfo,
                                 uint32_t added_mean, uint64_t added_seq)
{
  uint64_t total_sequence = ginfo->total_sequence + added_seq;

  if(total_sequence > 0)
  {
    ginfo->mean_read_length = (ginfo->mean_read_length * ginfo->total_sequence +
                               added_mean * added_seq) / total_sequence;
  }

  ginfo->total_sequence = total_sequence;
}

void graph_info_merge(GraphInfo *dst, const GraphInfo *src)
{
  // Update sample name
  const StrBuf *new_sample_name = &src->sample_name;
  if(new_sample_name->len > 0 && strcmp(new_sample_name->buff, "undefined") != 0)
    strbuf_set(&dst->sample_name, new_sample_name->buff);

  uint64_t total_sequence = dst->total_sequence + src->total_sequence;

  if(total_sequence > 0)
  {
    // Average error rates
    dst->seq_err
      = (dst->seq_err * dst->total_sequence +
         src->seq_err * src->total_sequence) /
        total_sequence;

    // Update mean read length
    dst->mean_read_length
      = (dst->mean_read_length * dst->total_sequence +
         src->mean_read_length * src->total_sequence) /
        total_sequence;
  }

  dst->total_sequence = total_sequence;

  // Update error cleaning
  error_cleaning_overwrite(&dst->cleaning, &src->cleaning);
}
