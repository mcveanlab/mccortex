
#include "global.h"

#include "seq_file.h"

#include "seq_reader.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "db_node.h"
#include "file_reader.h"

static boolean invalid_bases[256];
boolean fastq_qual_too_small = false, fastq_qual_too_big = false;

#define init_bases_check() do {                                                \
    memset(invalid_bases, 0, sizeof(invalid_bases));                           \
    fastq_qual_too_small = fastq_qual_too_big = false;                         \
  } while(0)

int prev_guess_fastq_format = -1;

// cut-offs:
//  > quality_cutoff valid
//  < homopolymer_cutoff valid

// Returns index of first kmer or r->seq.end if no kmers
size_t seq_contig_start(const read_t *r, long offset, uint32_t kmer_size,
                        int qual_cutoff, int hp_cutoff)
{
  size_t next_kmer, index = offset;
  while((next_kmer = index+kmer_size) <= r->seq.end)
  {
    // Check for invalid bases
    size_t i = next_kmer;
    while(i > index && is_base_char(r->seq.b[i-1])) i--;

    if(i > index) {
      index = i;
      continue;
    }

    // Check for low qual values
    if(qual_cutoff > 0 && r->qual.end > 0)
    {
      i = MIN2(next_kmer, r->qual.end);
      while(i > index && r->qual.b[i-1] > qual_cutoff) i--;

      if(i > index) {
        index = i;
        continue;
      }
    }

    // Check for homopolymer runs
    if(hp_cutoff > 0)
    {
      int run_length = 1;
      for(i = next_kmer-1; i > index; i--)
      {
        if(r->seq.b[i-1] == r->seq.b[i])
        {
          run_length++;
          if(run_length == hp_cutoff) break;
        }
        else run_length = 1;
      }

      if(i > index)
      {
        index = i;
        continue;
      }
    }

    return index;
  }

  return r->seq.end;
}

// *search_start is the next position to pass to seq_contig_start
size_t seq_contig_end(const read_t *r, size_t contig_start, uint32_t kmer_size,
                      int qual_cutoff, int hp_cutoff, size_t *search_start)
{
  size_t contig_end = contig_start+kmer_size;

  int hp_run = 1;
  if(hp_cutoff > 0)
  {
    // Get the length of the hp run at the end of the current kmer
    // kmer won't contain a run longer than hp_run-1
    while(r->seq.b[contig_end-1-hp_run] == r->seq.b[contig_end-1]) hp_run++;
  }

  for(; contig_end < r->seq.end; contig_end++)
  {
    if(!is_base_char(r->seq.b[contig_end]) ||
       (contig_end < r->qual.end && r->qual.b[contig_end] < qual_cutoff))
    {
      break;
    }

    // Check hp
    if(hp_cutoff > 0)
    {
      if(r->seq.b[contig_end] == r->seq.b[contig_end-1])
      {
        hp_run++;
        if(hp_run >= hp_cutoff) break;
      }
      else hp_run = 1;
    }
  }

  if(hp_cutoff > 0 && hp_run == hp_cutoff)
    *search_start = contig_end - hp_cutoff + 1;
  else
    *search_start = contig_end;

  return contig_end;
}

// returns offset of the first node found or -1 if no nodes were found
// Gaps collapsed down to a single HASH_NOT_FOUND
int get_nodes_from_read(const read_t *r, int qcutoff, int hp_cutoff,
                        const dBGraph *db_graph, dBNodeBuffer *list)
{
  size_t contig_start, contig_end, search_start = 0;
  uint32_t kmer_size = db_graph->kmer_size;

  BinaryKmer bkmer, tmp_key;
  Nucleotide nuc;
  hkey_t node, prev_node = HASH_NOT_FOUND;
  size_t next_base;
  int first_node_offset = -1;

  db_node_buf_ensure_capacity(list, list->len + r->seq.end);

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qcutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qcutoff, hp_cutoff, &search_start);

    if(prev_node != HASH_NOT_FOUND)
    {
      // Gap
      list->nodes[list->len] = HASH_NOT_FOUND;
      list->orients[list->len] = forward;
      list->len++;
      prev_node = HASH_NOT_FOUND;
    }

    const char *contig = r->seq.b + contig_start;
    size_t contig_len = contig_end - contig_start;

    binary_kmer_from_str(contig, kmer_size, bkmer);
    binary_kmer_right_shift_one_base(bkmer);

    for(next_base = kmer_size-1; next_base < contig_len; next_base++)
    {
      nuc = binary_nuc_from_char(contig[next_base]);
      binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
      db_node_get_key(bkmer, kmer_size, tmp_key);
      node = hash_table_find(&db_graph->ht, tmp_key);

      if(node != HASH_NOT_FOUND && first_node_offset == -1)
        first_node_offset = contig_start + next_base + 1 - kmer_size;

      if(prev_node != HASH_NOT_FOUND || node != HASH_NOT_FOUND)
      {
        list->nodes[list->len] = node;
        list->orients[list->len] = db_node_get_orientation(bkmer, tmp_key);
        list->len++;
        prev_node = node;
      }
    }
  }

  if(list->len > 0 && list->nodes[list->len-1] == HASH_NOT_FOUND)
    list->len--;

  return first_node_offset;
}

void process_new_read(read_t *r, int qmin, int qmax, const char *path)
{
  seq_read_to_uppercase(r);

  char *tmp;
  for(tmp = r->seq.b; *tmp != '\0'; tmp++)
  {
    if(!is_base_char(*tmp) && *tmp != 'N' && !invalid_bases[(size_t)*tmp])
    {
      warn("Invalid sequence [%c] [read: %s; path: %s]\n",
           *tmp, r->name.b, path);
      invalid_bases[(size_t)*tmp] = 1;
    }
  }

  // Check quality string
  if(r->qual.end > 0 && r->seq.end != r->qual.end)
  {
    warn("Quality string is not the same length as sequence [%i vs %i; read: %s; "
         "path: %s]", (int)r->qual.end, (int)r->seq.end, r->name.b, path);
  }

  // Check out-of-range qual string
  for(tmp = r->qual.b; *tmp != '\0'; tmp++)
  {
    int q = *tmp;
    if(q < qmin && !fastq_qual_too_small)
    {
      warn("FASTQ qual too small [%i < %i..%i; read: %s; path: %s]",
           q, qmin, qmax, r->name.b, path);
      fastq_qual_too_small = 1;
    }
    else if(q > qmax && !fastq_qual_too_big)
    {
      warn("FASTQ qual too big [%i > %i..%i; read: %s; path: %s]",
           q, qmin, qmax, r->name.b, path);
      fastq_qual_too_big = 1;
    }
  }
}

// Returns -1 if no quality scores
// Defaults to 0 if not recognisable (offset:33, min:33, max:126)
int guess_fastq_format(const char *path)
{
  // Detect fastq offset
  int min_qual = INT_MAX, max_qual = INT_MIN;
  int fmt = seq_guess_fastq_format(path, 500, &min_qual, &max_qual);

  if(prev_guess_fastq_format != fmt && fmt != -1)
  {
    int min = 0, max = 0;
    seq_get_qual_limits(path, 500, &min, &max);

    prev_guess_fastq_format = fmt;
    message("%s: Qual scores: %s [offset: %i, range: [%i,%i], sample: [%i,%i]]\n",
            path, FASTQ_FORMATS[fmt], FASTQ_OFFSET[fmt],
            FASTQ_MIN[fmt], FASTQ_MAX[fmt], min, max);
  }

  return fmt;
}

// Test min and max fastq scores for first 500bp
void test_file_qual_offset(const char *path, int qoffset, int qmax)
{
  int tmp_min = -1, tmp_max = -1;
  if(seq_get_qual_limits(path, 500, &tmp_min, &tmp_max) > 0)
  {
    if(tmp_min > qoffset + 20)
    {
      warn("Input file has min quality score %i but qoffset is set to %i: %s\n"
           "  Have you predefined an incorrect fastq offset? "
           "Or is cortex guessing it wrong?", tmp_min, qoffset, path);
    }
    else if(tmp_max > qmax + 20)
    {
      warn("Input file has max quality score %i but expected qmax is to %i: %s\n"
           "  Have you predefined an incorrect fastq offset? "
           "Or is cortex guessing it wrong?", tmp_min, qoffset, path);
    }
  }
}

void seq_parse_pe(const char *path1, const char *path2,
                  read_t *r1, read_t *r2, 
                  SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                  void (*read_func)(read_t *_r1, read_t *_r2,
                                    int _qoffset1, int _qoffset2,
                                    SeqLoadingPrefs *_prefs,
                                    SeqLoadingStats *_stats,
                                    void *_ptr),
                  void *reader_ptr)
{
  // Guess offset if needed
  int qoffset1 = prefs->ascii_fq_offset, qoffset2 = prefs->ascii_fq_offset;
  int qmin1 = prefs->ascii_fq_offset, qmin2 = prefs->ascii_fq_offset;
  int qmax1 = 126, qmax2 = 126;

  if(prefs->ascii_fq_offset == 0)
  {
    int fmt1, fmt2;
    if((fmt1 = guess_fastq_format(path1)) != -1) {
      qmin1 = FASTQ_MIN[fmt1];
      qmax1 = FASTQ_MAX[fmt1];
      qoffset1 = FASTQ_OFFSET[fmt1];
    }
    if((fmt2 = guess_fastq_format(path2)) != -1) {
      qmin2 = FASTQ_MIN[fmt2];
      qmax2 = FASTQ_MAX[fmt2];
      qoffset2 = FASTQ_OFFSET[fmt2];
    }
  }
  test_file_qual_offset(path1, qoffset1, qmax1);
  test_file_qual_offset(path2, qoffset2, qmax2);

  init_bases_check();

  seq_file_t *sf1, *sf2;
  if((sf1 = seq_open(path1)) == NULL)
  {
    fprintf(stderr, "Error cannot read: %s\n", path1);
    exit(EXIT_FAILURE);
  }
  if((sf2 = seq_open(path2)) == NULL)
  {
    fprintf(stderr, "Error cannot read: %s\n", path2);
    exit(EXIT_FAILURE);
  }

  int success1, success2;

  while(1)
  {
    success1 = seq_read(sf1, r1);
    success2 = seq_read(sf2, r2);

    if(success1 < 0) warn("input error: %s", sf1->path);
    if(success2 < 0) warn("input error: %s", sf2->path);
    if(!success1 != !success2) {
      warn("Different number of reads in pe files [%s; %s]\n",
           sf1->path, sf2->path);
    }
    if(success1 <= 0 || success2 <= 0) break;

    // pe
    process_new_read(r1, qmin1, qmax1, path1);
    process_new_read(r2, qmin2, qmax2, path2);
    // seq_read_reverse_complement(r2); // now done later
    read_func(r1, r2, qoffset1, qoffset2, prefs, stats, reader_ptr);
    stats->num_pe_reads += 2;
  }

  seq_close(sf1);
  seq_close(sf2);

  stats->num_files_loaded += 2;
}

void seq_parse_se(const char *path, read_t *r1, read_t *r2,
                  SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                  void (*read_func)(read_t *r1, read_t *r2,
                                    int qoffset1, int qoffset2,
                                    SeqLoadingPrefs *prefs,
                                    SeqLoadingStats *stats,
                                    void *ptr),
                  void *reader_ptr)
{
  // Guess offset if needed
  int qoffset = prefs->ascii_fq_offset;
  int qmin = prefs->ascii_fq_offset, qmax = 126;

  int format;
  if(prefs->ascii_fq_offset == 0 && (format = guess_fastq_format(path)) != -1)
  {
    qmin = FASTQ_MIN[format];
    qmax = FASTQ_MAX[format];
    qoffset = FASTQ_OFFSET[format];
  }
  test_file_qual_offset(path, qoffset, qmax);

  seq_file_t *sf;
  if((sf = seq_open(path)) == NULL)
  {
    fprintf(stderr, "Error cannot read: %s\n", path);
    exit(EXIT_FAILURE);
  }

  init_bases_check();

  if(!seq_is_sam(sf) && !seq_is_bam(sf))
  {
    // Single file with single-ended reads
    int s;
    while((s = seq_read(sf, r1)) > 0)
    {
      process_new_read(r1, qmin, qmax, path);
      read_func(r1, NULL, qoffset, 0, prefs, stats, reader_ptr);
      stats->num_se_reads++;
    }
    if(s < 0) warn("Input error: %s\n", sf->path);
  }
  else
  {
    // SAM / BAM single file may contain paired-reads
    char swap_reads = 0, have_mate = 0, mates_seen = 0;

    // For now assume sam format is 
    qmin=33;
    qmax=73;
    qoffset=33;

    while(1)
    {
      if(swap_reads)
      {
        // swap r2 into r1
        buffer_t tmp;
        SWAP(r1->name, r2->name, tmp);
        SWAP(r1->seq,  r2->seq,  tmp);
        SWAP(r1->qual, r2->qual, tmp);
      }
      else {
        int s = seq_read(sf, r1);
        if(s < 0) { warn("Input error: %s\n", sf->path); break; }
        if(s == 0) break;
      }

      if(seq_read(sf, r2) > 0)
      {
        have_mate = (strcmp(r1->name.b, r2->name.b) == 0);
        swap_reads = !have_mate;
      }
      else have_mate = swap_reads = 0;

      mates_seen = mates_seen || have_mate;

      if(have_mate)
      {
        // pe
        process_new_read(r1, qmin, qmax, path);
        process_new_read(r2, qmin, qmax, path);
        // seq_read_reverse_complement(r2); // now done later
        read_func(r1, r2, qoffset, qoffset, prefs, stats, reader_ptr);
        stats->num_pe_reads += 2;
      }
      else
      {
        // se
        process_new_read(r1, qmin, qmax, path);
        read_func(r1, NULL, qoffset, 0, prefs, stats, reader_ptr);
        stats->num_se_reads++;
      }
    }
  }

  seq_close(sf);
  stats->num_files_loaded++;
}

//
// Function to load into graph
//

char seq_reads_are_novel(const read_t *r1, const read_t *r2, dBGraph *db_graph,
                         int qual_cutoff1, int qual_cutoff2, int hp_cutoff)
{
  hkey_t node1 = HASH_NOT_FOUND, node2 = HASH_NOT_FOUND;
  Orientation or1 = forward, or2 = forward;
  BinaryKmer curr_kmer, tmp_key;
  boolean found1 = false, found2 = false;
  size_t start1, start2, kmer_size = db_graph->kmer_size;

  start1 = seq_contig_start(r1, 0, kmer_size, qual_cutoff1, hp_cutoff);
  start2 = seq_contig_start(r2, 0, kmer_size, qual_cutoff2, hp_cutoff);

  boolean got_kmer1 = (start1 < r1->seq.end);
  boolean got_kmer2 = (start2 < r2->seq.end);

  if(!got_kmer1 && !got_kmer2) return 1;

  // Look up first kmer
  if(got_kmer1)
  {
    binary_kmer_from_str(r1->seq.b+start1, kmer_size, curr_kmer);
    db_node_get_key(curr_kmer, kmer_size, tmp_key);
    node1 = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found1);
    or1 = db_node_get_orientation(curr_kmer, tmp_key);
  }

  // Look up second kmer
  if(got_kmer2)
  {
    binary_kmer_from_str(r2->seq.b+start2, kmer_size, curr_kmer);
    db_node_get_key(curr_kmer, kmer_size, tmp_key);
    node2 = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found2);
    or2 = db_node_get_orientation(curr_kmer, tmp_key);
  }

  // Each read gives no kmer or a duplicate kmer
  // used find_or_insert so if we have a kmer we have a graph node
  if((!got_kmer1 || db_node_has_read_start(db_graph, node1, or1)) &&
     (!got_kmer2 || db_node_has_read_start(db_graph, node2, or2)))
  {
    return 0;
  }

  // Read is novel
  if(got_kmer1) db_node_set_read_start(db_graph, node1, or1);
  if(got_kmer2) db_node_set_read_start(db_graph, node2, or2);
  return 1;
}

char seq_read_is_novel(const read_t *r, dBGraph *db_graph,
                       int qual_cutoff, int hp_cutoff)
{
  hkey_t node = HASH_NOT_FOUND;
  Orientation or;
  BinaryKmer bkmer, tmp_key;
  boolean found = false;

  size_t start = seq_contig_start(r, 0, db_graph->kmer_size,
                                  qual_cutoff, hp_cutoff);

  if(start == r->seq.end) return 1;

  binary_kmer_from_str(r->seq.b+start, db_graph->kmer_size, bkmer);
  db_node_get_key(bkmer, db_graph->kmer_size, tmp_key);
  node = hash_table_find_or_insert(&db_graph->ht, tmp_key, &found);
  or = db_node_get_orientation(bkmer, tmp_key);

  if(db_node_has_read_start(db_graph, node, or)) return 0;

  // Read is novel
  db_node_set_read_start(db_graph, node, or);
  return 1;
}

void load_read(const read_t *r, dBGraph *db_graph,
               int qual_cutoff, int hp_cutoff,
               Colour colour, SeqLoadingStats *stats)
{
  if(r->seq.end < db_graph->kmer_size) return;

  size_t search_start = 0;
  size_t contig_start, contig_end = 0;

  while((contig_start = seq_contig_start(r, search_start, db_graph->kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, db_graph->kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    size_t contig_len = contig_end - contig_start;

    // printf("contig: %.*s\n", (int)contig_len, r->seq.b+contig_start);

    // Load into graph
    BinaryKmer bkmer;
    binary_kmer_from_str(r->seq.b+contig_start, db_graph->kmer_size, bkmer);

    BinaryKmer tmp_key;
    hkey_t prev_node, curr_node;
    Orientation prev_or, curr_or;

    db_node_get_key(bkmer, db_graph->kmer_size, tmp_key);
    curr_node = db_graph_find_or_add_node(db_graph, tmp_key, colour);
    curr_or = db_node_get_orientation(bkmer, tmp_key);

    size_t i;
    for(i = contig_start+db_graph->kmer_size; i < contig_end; i++)
    {
      Nucleotide nuc = binary_nuc_from_char(r->seq.b[i]);
      binary_kmer_left_shift_add(bkmer, db_graph->kmer_size, nuc);

      prev_node = curr_node;
      prev_or = curr_or;

      db_node_get_key(bkmer, db_graph->kmer_size, tmp_key);
      curr_node = db_graph_find_or_add_node(db_graph, tmp_key, colour);
      curr_or = db_node_get_orientation(bkmer, tmp_key);

      db_graph_add_edge(db_graph, prev_node, curr_node, prev_or, curr_or);
    }

    // Update contig stats
    if(stats->readlen_count_array != NULL) {
      contig_len = MIN2(contig_len, stats->readlen_count_array_size-1);
      stats->readlen_count_array[contig_len]++;
    }
    stats->total_bases_loaded += contig_len;
    stats->kmers_loaded += contig_len + 1 - db_graph->kmer_size;
  }

  // contig_end == 0 if no contigs from this read
  if(contig_end == 0) stats->total_bad_reads++;
  else stats->total_good_reads++;
}

void seq_load_into_db_graph(read_t *r1, read_t *r2,
                            int qoffset1, int qoffset2,
                            SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                            void *ptr)
{
  (void)ptr;

  stats->total_bases_read += r1->seq.end;
  if(r2 != NULL) stats->total_bases_read += r2->seq.end;

  int qcutoff1 = prefs->quality_cutoff;
  int qcutoff2 = prefs->quality_cutoff;

  if(prefs->quality_cutoff > 0)
  {
    qcutoff1 += qoffset1;
    qcutoff2 += qoffset2;
  }

  int hp_cutoff = prefs->homopolymer_cutoff;

  if((r1->from_sam && r1->bam->core.flag & BAM_FDUP &&
      (r2 == NULL || (r2->from_sam && r2->bam->core.flag & BAM_FDUP))) ||
     (r2 == NULL && prefs->remove_dups_se &&
      !seq_read_is_novel(r1, prefs->db_graph, qcutoff1, hp_cutoff)) ||
     (r2 != NULL && prefs->remove_dups_pe &&
      !seq_reads_are_novel(r1, r2, prefs->db_graph, qcutoff1, qcutoff2, hp_cutoff)))
  {
    stats->total_dup_reads += (r2 == NULL ? 1 : 2);
    return;
  }

  load_read(r1, prefs->db_graph, qcutoff1, hp_cutoff, prefs->into_colour, stats);

  if(r2 != NULL)
  {
    load_read(r2, prefs->db_graph, qcutoff2, hp_cutoff, prefs->into_colour, stats);
  }
}
