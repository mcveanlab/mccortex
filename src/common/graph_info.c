#include "global.h"
#include "db_graph.h"

ErrorCleaning* error_cleaning_alloc(ErrorCleaning *ec)
{
  ec->cleaned_against_graph_name = strbuf_new();
  error_cleaning_init(ec);
  return ec;
}

void error_cleaning_dealloc(ErrorCleaning *ec)
{
  strbuf_free(ec->cleaned_against_graph_name);
}

void error_cleaning_init(ErrorCleaning *ec)
{
  ec->tip_clipping = ec->remv_low_cov_sups = ec->remv_low_cov_nodes = false;
  ec->remv_low_cov_sups_thresh = ec->remv_low_cov_nodes_thresh = 0;
  ec->cleaned_against_another_graph = false;
  strbuf_set(ec->cleaned_against_graph_name, "undefined");
}

GraphInfo *graph_info_alloc(GraphInfo *ginfo)
{
  int i;
  for(i = 0; i < NUM_OF_COLOURS; i++)
  {
    ginfo->sample_names[i] = strbuf_new();
    error_cleaning_alloc(ginfo->cleaning + i);
  }

  graph_info_init(ginfo);
  return ginfo;
}

void graph_info_dealloc(GraphInfo *ginfo)
{
  int i;
  for(i = 0; i < NUM_OF_COLOURS; i++)
  {
    strbuf_free(ginfo->sample_names[i]);
    error_cleaning_dealloc(ginfo->cleaning + i);
  }
}

void graph_info_init(GraphInfo *ginfo)
{
  uint32_t i;
  for(i = 0; i < NUM_OF_COLOURS; i++)
    graph_info_reset_one_colour(ginfo, i);

  ginfo->num_of_colours_loaded = 0;
  ginfo->num_of_shades_loaded = 0;
}

void graph_info_reset_one_colour(GraphInfo *ginfo, uint32_t colour)
{
  strbuf_set(ginfo->sample_names[colour], "undefined");
  graph_info_set_seq(ginfo, colour, 0);
  graph_info_set_mean_readlen(ginfo, colour, 0);
  ginfo->seq_err[colour] = 0.01;
  error_cleaning_init(ginfo->cleaning + colour);
}


void graph_info_set_tip_clipping(GraphInfo *ginfo, int colour)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->cleaning[colour].tip_clipping = true;
}

void graph_info_set_remv_low_cov_sups(GraphInfo *ginfo, int colour, Covg thresh)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->cleaning[colour].remv_low_cov_sups = true;
  ginfo->cleaning[colour].remv_low_cov_sups_thresh = thresh;
}

void graph_info_set_remv_low_cov_nodes(GraphInfo *ginfo, int colour, Covg thresh)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->cleaning[colour].remv_low_cov_nodes = true;
  ginfo->cleaning[colour].remv_low_cov_nodes_thresh = thresh;
}

void graph_info_set_seq_err(GraphInfo *ginfo, int colour, long double err)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->seq_err[colour] = err;
}

void graph_info_set_seq(GraphInfo *ginfo, int colour, uint64_t num_bp)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->total_sequence[colour] = num_bp;
}

void graph_info_set_mean_readlen(GraphInfo *ginfo, int colour, int len)
{
  assert(colour < NUM_OF_COLOURS);
  ginfo->mean_read_length[colour] = len;
}

void graph_info_update_mean_readlen_and_total_seq(GraphInfo *ginfo,
                                                  uint32_t colour,
                                                  uint32_t added_mean,
                                                  uint64_t added_seq)
{
  uint64_t total_sequence = ginfo->total_sequence[colour] + added_seq;

  if(total_sequence > 0)
  {
    // Update mean read length
    ginfo->mean_read_length[colour]
      = (ginfo->mean_read_length[colour] * ginfo->total_sequence[colour] +
         added_mean * added_seq) / total_sequence;
  }

  ginfo->total_sequence[colour] = total_sequence;
}


int get_mean_readlen_across_colours(GraphInfo *ginfo)
{
  int colour;

  long long alpha = 0;
  long long beta = 0;
  for(colour = 0; colour < NUM_OF_COLOURS; colour++)
  {
    alpha += ginfo->total_sequence[colour] * ginfo->mean_read_length[colour];
    beta += ginfo->total_sequence[colour];
  }
  return (int)(alpha / beta);
}

void read_estimated_seq_errors_from_file(GraphInfo *ginfo, FILE *fp)
{
  StrBuf *line = strbuf_new();
  int col;

  for(col = 0; strbuf_readline(line, fp) > 0; col++)
  {
    strbuf_chomp(line);
    if(line->len > 0)
    {
      if(col == NUM_OF_COLOURS) die("Too many colours in seq error file");

      // DEV: catch error if number is not valid
      ginfo->seq_err[col] = (long double)strtod(line->buff, NULL);
    }
  }

  strbuf_free(line);
}

void print_seq_err_rates_to_screen(GraphInfo *ginfo)
{
  printf("Setting the following per-colour sequencing error rates "
         "(used only for genotyping):\nColour\tRate\n");
  int i;
  for(i = 0; i < NUM_OF_COLOURS; i++)
  {
    printf("%d\t%.3Lf\n", i, ginfo->seq_err[i]);
  }
}
