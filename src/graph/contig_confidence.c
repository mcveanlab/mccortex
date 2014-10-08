#include "global.h"
#include "contig_confidence.h"
#include "util.h"
#include "file_util.h"

#include <fcntl.h> // O_CREAT et al
#include <math.h>

static inline int contig_hist_cmp(const void *aa, const void *bb)
{
  const ContigHistItem *a = *(const ContigHistItem *const*)aa;
  const ContigHistItem *b = *(const ContigHistItem *const*)bb;
  if(a->contig_len != b->contig_len)
    return (a->contig_len < b->contig_len ? -1 : 1);
  if(a->num != b->num)
    return (a->num < b->num ? -1 : 1);
  return 0;
}

// Merge duplicate entries together
static void contig_hist_sort_merge(ContigConfidenceTable *table)
{
  size_t i, j;
  ContigHistItem *arr = table->hist.data;
  size_t len = table->hist.len;

  if(len == 0) return;

  qsort(arr, len, sizeof(ContigHistItem), contig_hist_cmp);

  // Merge entries where contig_len matches
  for(i = j = 1; j < len; j++) {
    if(arr[i-1].contig_len == arr[j].contig_len) {
      arr[i-1].num += arr[j].num;
    } else {
      arr[i] = arr[j];
      i++;
    }
  }

  table->hist.len = i;
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

static void _update_table_with_length_count(ContigConfidenceTable *conf_table,
                                            size_t genome_size,
                                            size_t contig_len, size_t num)
{
  size_t i, init_len = conf_table->table.len;
  if(contig_len+1 > init_len) {
    double_buf_capacity(&conf_table->table, contig_len+1);
    conf_table->table.len = contig_len+1;
    for(i = init_len; i <= contig_len; i++) conf_table->table.data[i] = 0;
  }

  for(i = 1; i <= contig_len; i++) {
    double covg = (double)(contig_len * num) / genome_size;
    double old_conf = conf_table->table.data[i];
    double new_conf = 1.0 - (1.0 - old_conf) *
                            (1.0 - calc_confid(covg, contig_len, i));
    // printf("  contig: %zu num: %zu genome: %zu covg: %.2f\n", contig_len, num,
    //        genome_size, covg);
    // printf("updating %zu from %.2f -> %.2f [%.2f]\n", i, old_conf, new_conf, calc_confid(covg, contig_len, i));
    conf_table->table.data[i] = new_conf;
  }
}

// Call conf_table_destroy to release memory after calling this function
// fh is not closed
// format should be:
//
//   # ignored line       <- comment line
//   contig_len_bp\tlen   <- header line
//   <num>\t<count>       <- entry
//   # another comment
//   <num>\t<count>
//   ...
//
// Call conf_table_destroy to release memory after calling this function
void conf_table_load_csv(ContigConfidenceTable *conf_table,
                         FILE *fh, const char *path)
{
  StrBuf line;
  strbuf_alloc(&line, 1024);
  bool seen_header = false;
  size_t lineno, prev_contig_len = 0, contig_len = 0, num = 0;

  for(lineno = 1; strbuf_reset_readline(&line, fh); lineno++) {
    strbuf_chomp(&line);
    if(line.end > 0 && line.b[0] != '#') {
      char *a = line.b, *b = strchr(a, '\t');
      if(!b) die("Invalid CSV line [%s:%zu]", path, lineno);
      *b = '\0'; b++;

      a = string_trim(a);
      b = string_trim(b);

      if(!*a || !*b) die("Invalid CSV line [%s:%zu]", path, lineno);

      // Try to parse
      if(!seen_header) seen_header = true;
      else if(!parse_entire_size(a, &contig_len) || !parse_entire_size(b, &num))
        die("Invalid CSV line [%s:%zu]", path, lineno);
      else if(contig_len <= prev_contig_len && prev_contig_len > 0)
        die("Invalid CSV line [%s:%zu]", path, lineno);
      else {
        ContigHistItem item = {.contig_len = contig_len, .num = num};
        contig_hist_add(&conf_table->hist, item);
        prev_contig_len = contig_len;
      }
    }
  }
  strbuf_dealloc(&line);

  contig_hist_sort_merge(conf_table);
}

void conf_table_calc_csv(ContigConfidenceTable *conf_table, size_t genome_size)
{
  size_t i;

  for(i = 0; i < conf_table->hist.len; i++)
  {
    ContigHistItem item = conf_table->hist.data[i];
    _update_table_with_length_count(conf_table, genome_size,
                                    item.contig_len, item.num);
  }

  contig_hist_dealloc(&conf_table->hist);
}

// Call conf_table_destroy to release memory after calling this function
void conf_table_calc(ContigConfidenceTable *conf_table,
                     size_t max_read_len, double avg_bp_covg)
{
  size_t i;

  status("Confidences for max. read length %zu and expected coverage %.2fX",
         max_read_len, avg_bp_covg);

  double_buf_capacity(&conf_table->table, max_read_len+1);
  conf_table->table.len = max_read_len+1;

  for(i = 0; i <= max_read_len; i++) {
    conf_table->table.data[i] = calc_confid(avg_bp_covg, max_read_len, i);
  }
}

void conf_table_destroy(ContigConfidenceTable *conf_table)
{
  contig_hist_dealloc(&conf_table->hist);
  double_buf_dealloc(&conf_table->table);
  memset(conf_table, 0, sizeof(ContigConfidenceTable));
}

double conf_table_lookup(const ContigConfidenceTable *table, size_t rlen)
{
  return rlen < table->table.len ? table->table.data[rlen] : -1;
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
