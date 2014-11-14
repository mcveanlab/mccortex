#include "global.h"
#include "contig_confidence.h"
#include "util.h"
#include "file_util.h"

#include <fcntl.h> // O_CREAT et al
#include <math.h>

static void conf_table_capacity(ContigConfidenceTable *table, size_t max_read)
{
  size_t i, init_len = table->table.len;
  double_buf_extend(&table->table, (max_read+1)*table->ncols);
  for(i = init_len; i < table->table.len; i++) table->table.data[i] = 0.0;
}

// expl(x) is the same as doing powl(M_E, x) a.k.a. 2**x for type long double
static double calc_confid(double bp_covg_depth, size_t read_length_bp,
                          size_t kmer_size)
{
  // printf("  covg: %.2f read: %zu bp contig: %zu bp\n",
  //        bp_covg_depth, read_length_bp, kmer_size);
  double lambda = bp_covg_depth / read_length_bp;
  double read_kmers = read_length_bp - kmer_size + 1;
  double power = (1.0 - expl(-lambda * read_kmers)) *
                 expl(-lambda * expl(-lambda * read_kmers));
  return power;
}

static void _update_table_with_length_count(ContigConfidenceTable *table,
                                            size_t col, size_t genome_size,
                                            size_t contig_len, size_t num)
{
  ctx_assert(col < table->ncols);

  size_t i;
  conf_table_capacity(table, contig_len);

  for(i = 1; i <= contig_len; i++) {
    double covg = (double)(contig_len * num) / genome_size;
    double old_conf = table->table.data[i];
    double new_conf = 1.0 - ((1.0 - old_conf) *
                             (1.0 - calc_confid(covg, contig_len, i)));
    // printf("  contig: %zu num: %zu genome: %zu covg: %.2f\n", contig_len, num,
    //        genome_size, covg);
    // printf("updating %zu from %.2f -> %.2f [%.2f]\n", i, old_conf, new_conf,
    //        calc_confid(covg, contig_len, i));
    table->table.data[i*table->ncols+col] = new_conf;
  }
}

void conf_table_update_hist(ContigConfidenceTable *table,
                            size_t col, size_t genome_size,
                            size_t *contig_hist, size_t hist_len)
{
  ctx_assert(col < table->ncols);
  size_t i;

  for(i = 1; i < hist_len; i++) {
    if(contig_hist[i]) {
      _update_table_with_length_count(table, col, genome_size,
                                      i, contig_hist[i]);
    }
  }
}

// Call conf_table_dealloc to release memory after calling this function
void conf_table_calc(ContigConfidenceTable *table, size_t col,
                     size_t max_read_len, double avg_bp_covg)
{
  ctx_assert(col < table->ncols);
  size_t i;

  status("Confidences for max. read length %zu and expected coverage %.2fX",
         max_read_len, avg_bp_covg);

  conf_table_capacity(table, max_read_len);

  for(i = 1; i <= max_read_len; i++) {
    table->table.data[i*table->ncols+col] = calc_confid(avg_bp_covg, max_read_len, i);
  }
}

void conf_table_alloc(ContigConfidenceTable *table, size_t ncols)
{
  double_buf_alloc(&table->table, 512*ncols);
  table->ncols = ncols;
}

void conf_table_dealloc(ContigConfidenceTable *table)
{
  double_buf_dealloc(&table->table);
  memset(table, 0, sizeof(ContigConfidenceTable));
}

double conf_table_lookup(const ContigConfidenceTable *table,
                         size_t col, size_t dist)
{
  size_t idx = table->ncols * dist + col;
  ctx_assert(col < table->ncols);
  ctx_assert(idx < table->table.len);
  return table->table.data[idx];
}

void conf_table_print(const ContigConfidenceTable *table, FILE *fh)
{
  size_t i, j, endj, num_rows = table->table.len / table->ncols;

  fprintf(fh, "gap_dist");
  for(i = 0; i < table->ncols; i++)
    fprintf(fh, "\tconfidence_%zu", i);
  fprintf(fh, "\n");

  if(table->table.len == 0) return;

  for(i = 1, j = table->ncols; i < num_rows; i++) {
    fprintf(fh, "%zu", i);
    for(endj = j + table->ncols; j < endj; j++)
      fprintf(fh, "\t%.5f", table->table.data[j]);
    fprintf(fh, "\n");
  }
}

void conf_table_save(const ContigConfidenceTable *table, const char *path)
{
  int fd;
  FILE *fh;
  if(strcmp(path,"-") == 0) conf_table_print(table, stdout);
  else if((fd = futil_create_file(path, O_CREAT | O_EXCL | O_WRONLY)) < 0 ||
          (fh = fdopen(fd, "w")) == NULL)
  {
    warn("Cannot open file to save CSV table: %s", path);
    if(fd >= 0) close(fd);
  }
  else {
    conf_table_print(table, fh);
    fclose(fh);
  }
}
