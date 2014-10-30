#include "global.h"
#include "contig_confidence.h"
#include "util.h"
#include "file_util.h"

#include <fcntl.h> // O_CREAT et al
#include <math.h>

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

static void _update_table_with_length_count(ContigConfidenceTable *conf_table,
                                            size_t genome_size,
                                            size_t contig_len, size_t num)
{
  size_t i, init_len = conf_table->table.len;
  if(contig_len+1 > init_len) {
    double_buf_extend(&conf_table->table, contig_len+1);
    for(i = init_len; i <= contig_len; i++) conf_table->table.data[i] = 0;
  }

  for(i = 1; i <= contig_len; i++) {
    double covg = (double)(contig_len * num) / genome_size;
    double old_conf = conf_table->table.data[i];
    double new_conf = 1.0 - ((1.0 - old_conf) *
                             (1.0 - calc_confid(covg, contig_len, i)));
    // printf("  contig: %zu num: %zu genome: %zu covg: %.2f\n", contig_len, num,
    //        genome_size, covg);
    // printf("updating %zu from %.2f -> %.2f [%.2f]\n", i, old_conf, new_conf,
    //        calc_confid(covg, contig_len, i));
    conf_table->table.data[i] = new_conf;
  }
}

void conf_table_update_hist(ContigConfidenceTable *conf_table, size_t genome_size,
                            size_t *contig_hist, size_t hist_len)
{
  size_t i;

  for(i = 1; i < hist_len; i++) {
    if(contig_hist[i]) {
      _update_table_with_length_count(conf_table, genome_size,
                                      i, contig_hist[i]);
    }
  }
}

// Call conf_table_destroy to release memory after calling this function
void conf_table_calc(ContigConfidenceTable *conf_table,
                     size_t max_read_len, double avg_bp_covg)
{
  size_t i;

  status("Confidences for max. read length %zu and expected coverage %.2fX",
         max_read_len, avg_bp_covg);

  double_buf_extend(&conf_table->table, max_read_len+1);

  for(i = 1; i <= max_read_len; i++) {
    conf_table->table.data[i] = calc_confid(avg_bp_covg, max_read_len, i);
  }
}

void conf_table_destroy(ContigConfidenceTable *conf_table)
{
  double_buf_dealloc(&conf_table->table);
  memset(conf_table, 0, sizeof(ContigConfidenceTable));
}

double conf_table_lookup(const ContigConfidenceTable *table, size_t rlen)
{
  if(rlen >= table->table.len) die("%zu > %zu", rlen, table->table.len);
  return table->table.data[rlen];
}

void conf_table_print(const ContigConfidenceTable *table, FILE *fh)
{
  fprintf(fh, "gap_dist\tconfidence\n");
  if(table->table.len == 0) return;

  size_t i, max_read = table->table.len-1;
  for(i = 1; i <= max_read; i++)
    fprintf(fh, "%zu\t%f\n", i, table->table.data[i]);
  fprintf(fh, "\n");
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
