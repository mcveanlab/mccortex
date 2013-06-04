
#include "global.h"

#include <time.h>
#include <ctype.h> // isspace

#include "khash.h"
#include "string_buffer.h"

#include "global.h"
#include "binary_kmer.h"
#include "util.h"
#include "file_util.h"
#include "file_reader.h"
#include "delta_array.h"

#define PROC_MAXK 100

static const char usage[] =
"usage: proc_calls [options] <path_calls.bubbles> <out.base>\n"
"  Produces files <out.base>.vcf and <out.base>.5pflank.fa\n\n"
"  Options:\n"
"    --nobubbles <col>  Filter out e.g. ref bubbles\n";

#define load_assert(cond,sb) \
  if(!(cond)) { die("Loading err ["QUOTE(cond)"]: '%s'", (sb)->buff); }

/* alleles */

typedef struct VarBranch VarBranch;

struct VarBranch {
  StrBuf *seq;
  uint32_t num_nodes;
  delta_array_t *covgs;
  VarBranch *next;
};

static VarBranch* branch_new(int num_samples)
{
  int i;
  assert(num_samples > 0);
  VarBranch *branch = malloc(sizeof(VarBranch));
  branch->seq = strbuf_new();
  branch->covgs = malloc(num_samples * sizeof(delta_array_t));
  for(i = 0; i < num_samples; i++) delta_arr_alloc(&(branch->covgs[i]));
  branch->num_nodes = 0;
  branch->next = NULL;
  return branch;
}

static void branch_free(VarBranch *branch, int num_samples)
{
  strbuf_free(branch->seq);
  int i;
  for(i = 0; i < num_samples; i++) delta_arr_dealloc(&(branch->covgs[i]));
  free(branch->covgs);
  free(branch);
}

/* Variants */

typedef struct
{
  char key[2*PROC_MAXK+1]; // concat of key0, key1 (keys of kmer0, kmer1)
  StrBuf *flank5p, *flank3p;
  VarBranch *first_allele; // Linkedlist
  uint32_t var_num;
  uint32_t shift_left, shift_right;
} Var;

static Var* var_new(uint32_t kmer_size)
{
  Var *var = malloc(sizeof(Var));
  var->key[0] = var->key[2*kmer_size] = 0;
  var->flank5p = strbuf_new();
  var->flank3p = strbuf_new();
  var->first_allele = NULL;
  var->var_num = var->shift_left = var->shift_right = 0;
  return var;
}

static void var_free(Var *var, int num_samples)
{
  strbuf_free(var->flank5p);
  strbuf_free(var->flank3p);

  VarBranch *allele, *next;
  for(allele = var->first_allele; allele != NULL; allele = next)
  {
    next = allele->next;
    branch_free(allele, num_samples);
  }
  free(var);
}

static void var_revcmp(Var *var, int num_samples)
{
  StrBuf *tmpbuf;
  SWAP(var->flank5p, var->flank3p, tmpbuf);
  reverse_complement_str(var->flank5p->buff, var->flank5p->len);
  reverse_complement_str(var->flank3p->buff, var->flank3p->len);

  uint32_t tmpshift;
  SWAP(var->shift_left, var->shift_right, tmpshift);

  VarBranch *allele;
  int i;
  for(allele = var->first_allele; allele != NULL; allele = allele->next) {
    reverse_complement_str(allele->seq->buff, allele->seq->len);
    for(i = 0; i < num_samples; i++) delta_arr_reverse(&(allele->covgs[i]));
  }
}

// Merge var1 into var0, free var1
static void var_merge(Var *var0, Var *var1, int num_samples)
{
  // Compare sorted linked lists of alleles to count number of novel alleles
  VarBranch *allele0 = var0->first_allele;
  VarBranch *allele1 = var1->first_allele;

  VarBranch null_allele;
  null_allele.next = NULL;

  VarBranch *prev = &null_allele;

  while(allele0 != NULL && allele1 != NULL)
  {
    int cmp = strcmp(allele0->seq->buff, allele1->seq->buff);

    if(cmp <= 0) {
      prev->next = allele0;
      allele0 = allele0->next;

      if(cmp == 0)
      {
        VarBranch *tmp = allele1->next;
        branch_free(allele1, num_samples);
        allele1 = tmp;
      }
    }
    else {
      prev->next = allele1;
      allele1 = allele1->next;
    }
    prev = prev->next;
  }

  if(allele0 != NULL) prev->next = allele0;
  else if(allele1 != NULL) prev->next = allele1;
  else prev->next = NULL;

  var0->first_allele = null_allele.next;
  var1->first_allele = NULL;
}

// Trim from right side, add to flank3p
static void var_trim_alleles(Var *var, StrBuf *flank3p)
{
  VarBranch *branch0 = var->first_allele, *branch1 = branch0->next, *branch;
  size_t trim;
  size_t len0 = branch0->seq->len;
  const char *allele0 = branch0->seq->buff;
  boolean match = 1;

  for(trim = 0; trim < len0; trim++)
  {
    char lastc = allele0[len0-1-trim];

    // loop over alleles the other alleles
    for(branch = branch1; branch != NULL; branch = branch->next)
    {
      if(trim >= branch->seq->len ||
         branch->seq->buff[branch->seq->len-1-trim] != lastc)
      {
        match = 0;
        break;
      }
    }
    if(!match) break;
  }

  if(trim > 0)
  {
    strbuf_insert(flank3p, 0, allele0+len0-trim, trim);

    for(branch = branch0; branch != NULL; branch = branch->next)
      strbuf_shrink(branch->seq, branch->seq->len - trim);
  }
}

static void var_set_flank_shifts(Var *var)
{
  // May need to left shift 3p start
  uint32_t min_left = UINT_MAX, min_right = UINT_MAX;
  uint32_t j, k, max_dist;
  VarBranch *branch;
  const StrBuf *fl5p = var->flank5p;
  const StrBuf *fl3p = var->flank3p;

  for(branch = var->first_allele;
      branch != NULL && min_left + min_right > 0;
      branch = branch->next)
  {
    StrBuf *allele = branch->seq;

    if(allele->len > 0)
    {
      // no point checking futher than min_left
      max_dist = MIN3(allele->len, fl5p->len, min_left);
      j = k = 0;

      while(j < max_dist &&
            allele->buff[allele->len-j-1] == fl5p->buff[fl5p->len-j-1])
      { j++; }

      if(j == allele->len)
      {
        // Can continue along 5p flank
        max_dist = MIN2(min_left, fl5p->len) - j;
        while(k < max_dist &&
              fl5p->buff[fl5p->len-j-k-1] == fl5p->buff[fl5p->len-k-1])
        { k++; }
      }

      min_left = MIN2(min_left, j+k);
    
      // Now go right
      max_dist = MIN3(allele->len, fl3p->len, min_right);
      j = k = 0;

      while(j < max_dist && allele->buff[j] == fl3p->buff[j]) j++;

      if(j == allele->len)
      {
        // Can continue along 3p flank
        max_dist = MIN2(min_right, fl3p->len) - j;
        while(k < max_dist && fl3p->buff[j+k] == fl3p->buff[k]) k++;
      }

      min_right = MIN2(min_right, j+k);
    }
  }

  var->shift_left = min_left;
  var->shift_right = min_right;
}

// Warning: does not null terminate
static void str_get_key(char *key, const char *kmer, uint32_t kmer_size)
{
  memcpy(key, kmer, kmer_size*sizeof(char));
  reverse_complement_str(key, kmer_size);

  if(strncmp(kmer, key, kmer_size) < 0)
    memcpy(key, kmer, kmer_size*sizeof(char));
}

static void var_set_keys(Var *var, uint32_t kmer_size)
{
  StrBuf *fl5p = var->flank5p, *fl3p = var->flank3p;

  printf("5p:%s 3p:%s\n", fl5p->buff, fl3p->buff);

  char *key0 = var->key;
  char *key1 = var->key + kmer_size;

  str_get_key(key0, fl5p->buff + fl5p->len - kmer_size, kmer_size);

  if(var->shift_left > 0)
  {
    // Take from 5p flank
    char kmer1[PROC_MAXK+1];
    size_t cpy_len = MIN2(var->shift_left, kmer_size);
    memcpy(kmer1, fl5p->buff + fl5p->len - var->shift_left, cpy_len);
    if(cpy_len < kmer_size) memcpy(kmer1+cpy_len, fl3p->buff, kmer_size-cpy_len);
    str_get_key(key1, kmer1, kmer_size);
  }
  else {
    str_get_key(key1, fl3p->buff, kmer_size);
  }

  var->key[kmer_size+kmer_size] = '\0';
}

int var_cmp_alleles(const void *a, const void *b)
{
  const VarBranch *const a0 = *(const VarBranch *const *)a;
  const VarBranch *const a1 = *(const VarBranch *const *)b;

  return strcmp(a0->seq->buff, a1->seq->buff);
}

// Only applies if alleles have been trimmed of padding bases
static char var_is_snp(Var *var)
{
  VarBranch *allele;
  for(allele = var->first_allele; allele != NULL; allele = allele->next)
    if(allele->seq->len != 1) return 0;
  return 1;
}

/* CallHeader */

typedef struct {
  char **tags;
  char **values;
  size_t hcap, hlines;
  size_t num_samples, kmer_size;
  char **sample_names;
  char is_old_bc;
} CallHeader;

static CallHeader* header_new()
{
  CallHeader *ch = malloc(sizeof(CallHeader));
  ch->num_samples = ch->kmer_size = 0;
  ch->sample_names = NULL;
  ch->hcap = 128;
  ch->hlines = 0;
  ch->tags = malloc(ch->hcap * sizeof(*ch->tags));
  ch->values = malloc(ch->hcap * sizeof(*ch->tags));
  ch->is_old_bc = 0;
  return ch;
}

static void header_free(CallHeader *ch)
{
  size_t i;
  for(i = 0; i < ch->hlines; i++) { free(ch->tags[i]); free(ch->values[i]); }
  for(i = 0; i < ch->num_samples; i++) { free(ch->sample_names[i]); }
  free(ch->tags);
  free(ch->values);
  free(ch);
}

static void header_add(CallHeader *ch, char *tag, char *value)
{
  if(ch->hlines == ch->hcap) {
    ch->hcap <<= 1;
    ch->tags = realloc(ch->tags, ch->hcap * sizeof(char*));
    ch->values = realloc(ch->values, ch->hcap * sizeof(char*));
  }

  ch->tags[ch->hlines] = tag;
  ch->values[ch->hlines] = value;
  ch->hlines++;
}

static void print_vcf_header(gzFile vcf, CallHeader *ch, size_t argc, char **argv)
{
  size_t i, j;
  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  char cwd[PATH_MAX + 1];

  gzprintf(vcf, "##fileformat=VCFv4.1\n");
  gzprintf(vcf, "##fileDate=%s\n", datestr);
  gzprintf(vcf, "##reference=unplaced\n");
  gzprintf(vcf, "##phasing=none\n");

  gzprintf(vcf, "##procCmd=%s", argv[0]);
  if(file_reader_get_current_dir(cwd) != NULL) {
    gzprintf(vcf, "##procCwd=%s", cwd);
  }
  for(j = 1; j < argc; j++) gzprintf(vcf, " %s", argv[j]);
  gzputc(vcf, '\n');
  gzprintf(vcf, "##procDate=%s\n", datestr);

  for(i = 0; i < ch->hlines; i++) {
    if(strcasecmp(ch->tags[i],"fileformat") != 0)
      gzprintf(vcf, "##%s=%s\n", ch->tags[i], ch->values[i]);
  }

  gzprintf(vcf, "##INFO=<ID=LF,Number=1,Type=String,Description=\"Left flank\">\n");
  gzprintf(vcf, "##INFO=<ID=RF,Number=1,Type=String,Description=\"Right flank\">\n");
  gzprintf(vcf, "##INFO=<ID=BN,Number=.,Type=Integer,Description=\"Branch nodes; length in kmers of the alleles\">\n");
  gzprintf(vcf, "##INFO=<ID=BUB,Number=1,Type=String,Description=\"Cortex bubble index\">\n");
  gzprintf(vcf, "##FORMAT=<ID=COVG,Number=.,Type=Integer,Description=\"Number of read arrivals\">\n");
  gzprintf(vcf, "##contig=<ID=un,length=1000000,assembly=None>\n");

  gzprintf(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  for(i = 0; i < ch->num_samples; i++) {
    if(strcmp(ch->sample_names[i], "undefined") == 0)
      gzprintf(vcf, "\tsample%zu", i);
    else
      gzprintf(vcf, "\t%s", ch->sample_names[i]);
  }

  gzprintf(vcf, "\n");
}

static void parse_shaded_header(gzFile in, CallHeader* ch)
{
  StrBuf *line = strbuf_new(), *tmp = strbuf_new();
  int c;

  size_t sampleid = 0;

  // Parse first line ##fileformat=CTXv1.0
  char fftag[] = "##fileformat=";
  load_assert(strbuf_gzreadline_nonempty(line, in) != 0, line);
  strbuf_chomp(line);
  load_assert(strncasecmp(line->buff,fftag,strlen(fftag)) == 0, line);
  char *fileformat = line->buff+strlen(fftag);
  if(strcmp(fileformat, "CTXv1.0"))
    warn("unexpected fileformat '%s'", fileformat);

  while((c = gzgetc(in)) != -1 && c == '#')
  {
    strbuf_reset(line);
    strbuf_append_char(line, c);
    strbuf_gzreadline(line, in);
    strbuf_chomp(line);

    char *eq = strchr(line->buff, '=');
    if(line->buff[1] != '#' || eq == NULL) die("parse error: %s", line->buff);

    *eq = 0;
    char *tag = strdup(line->buff+2);
    char *value = strdup(eq+1);
    *eq = '=';

    if(strcasecmp(tag, "ctxCmd") == 0)
    {
      // Get kmer size from cmd line
      char *str = strstr(value, "--kmer_size");
      if(str == NULL || sscanf(str, "--kmer_size %zu", &ch->kmer_size) != 1) {
        die("Couldn't read kmer size line: %s=%s", tag, value);
      }
    }
    if(strcasecmp(tag, "ctx_calling_in_colours") == 0)
    {
      if(ch->num_samples != 0) die("duplicate header line: %s", line->buff);

      ch->num_samples = 1;
      char *str = value;
      while((str = strchr(str, ',')) != NULL) { str++; ch->num_samples++; }

      ch->sample_names = malloc(ch->num_samples*sizeof(char*));
    }
    else if(strncasecmp(tag, "sample", 6) == 0)
    {
      if(ch->sample_names == NULL) die("missing 'ctx_calling_in_colours=' line");

      strbuf_ensure_capacity(tmp, line->len);
      if(sscanf(line->buff, "##SAMPLE=<ID=%[^,],", tmp->buff) != 1)
        die("parse error on line: %s", line->buff);

      if(sampleid >= ch->num_samples)
        die("more samples than expected in header");

      tmp->len = strlen(tmp->buff);
      ch->sample_names[sampleid++] = strbuf_as_str(tmp);
    }

    header_add(ch, tag, value);
  }

  if(c != -1) gzungetc(c, in);

  if(ch->kmer_size == 0) die("Cannot determine kmer_size in call file");
  if(ch->num_samples == 0) die("Cannot determine num_samples in call file");

  strbuf_free(tmp);
  strbuf_free(line);
}

#define header_assert(cond,sb)                                \
  if(!(cond)) { die("Loading err ["QUOTE(cond)"]: %s", (sb)->buff); }

// Fabricate header by looking at first entry
static void synthesize_bubble_caller_header(gzFile fh, CallHeader *ch)
{
  // Need to guess kmer_size, num_of_samples and set samplenames
  ch->kmer_size = 31;
  ch->num_samples = 0;

  StrBuf *line = strbuf_new();

  if(!strbuf_gzreadline_nonempty(line, fh)) {
    strbuf_free(line);
    return;
  }

  int varnum, flank5plen, kmer_size;
  int h = sscanf(line->buff, ">var_%i_5p_flank length:%i INFO:KMER=%i",
                 &varnum, &flank5plen, &kmer_size);

  header_assert(h == 3, line);
  ch->kmer_size = kmer_size;

  // Skip next 7 lines
  size_t i;
  for(i = 0; i < 7; i++) {
    header_assert(strbuf_gzskipline(fh) > 0, line);
  }

  strbuf_reset(line);
  if(!strbuf_gzreadline_nonempty(line, fh) || line->buff[0] == '>') {
    strbuf_free(line);
    return;
  }

  strbuf_chomp(line);
  header_assert(strcmp(line->buff, "branch1 coverages") == 0, line);
  strbuf_reset_gzreadline(line, fh);
  header_assert(strncasecmp(line->buff, "Covg in Colour", 14) == 0, line);

  do
  {
    ch->num_samples++;
    strbuf_gzskipline(fh);
    strbuf_reset_gzreadline(line, fh);
  } while(strncasecmp(line->buff, "Covg in Colour", 14) == 0);

  ch->is_old_bc = 1;

  ch->sample_names = malloc(ch->num_samples * sizeof(char*));
  for(i = 0; i < ch->num_samples; i++)
  {
    strbuf_reset(line);
    strbuf_sprintf(line, "sample%zu", i);
    ch->sample_names[i] = strbuf_dup(line);
  }

  strbuf_free(line);
}


/* CallReader */

typedef struct
{
  CallHeader *ch;
  // Temp variables used for loading
  StrBuf *line;
  uint32_t *covgs, covgs_cap;
  VarBranch **alleles;
  uint32_t num_alleles, alleles_cap;
} CallReader;

static CallReader* reader_new(CallHeader *ch)
{
  CallReader *cr = malloc(sizeof(CallReader));
  cr->ch = ch;
  cr->line = strbuf_new();
  cr->covgs_cap = 1024;
  cr->covgs = malloc(cr->covgs_cap * sizeof(*cr->covgs));
  cr->num_alleles = 0;
  cr->alleles_cap = 32;
  cr->alleles = malloc(cr->alleles_cap * sizeof(*cr->alleles));
  return cr;
}

static void reader_free(CallReader *cr)
{
  strbuf_free(cr->line);
  free(cr->covgs);
  free(cr->alleles);
  free(cr);
}

static void reader_alleles_array_to_linkedlist(const CallReader *cr, Var *var)
{
  size_t i;
  for(i = 0; i+1 < cr->num_alleles; i++) cr->alleles[i]->next = cr->alleles[i+1];
  cr->alleles[cr->num_alleles-1]->next = NULL;
  var->first_allele = cr->alleles[0];
}

static void reader_sort_and_link_alleles(CallReader *cr, Var *var)
{
  if(cr->num_alleles == 0) die("No alleles passed");
  qsort(cr->alleles, cr->num_alleles, sizeof(VarBranch*), var_cmp_alleles);
  reader_alleles_array_to_linkedlist(cr, var);
}

//
// Clean up variant
//
static void reader_clean_up_var(CallReader *cr, Var *var)
{
  // Set alleles
  reader_alleles_array_to_linkedlist(cr, var);

  #ifdef DEBUG
  printf("\n\n");
  printf("var_%u\n 5p:%s\n 3p:%s\n", var->var_num,
         var->flank5p->buff, var->flank3p->buff);
  size_t i;
  for(i = 0; i < cr->num_alleles; i++)
    printf(" %zu:%s\n", i, cr->alleles[i]->seq->buff);
  #endif

  uint32_t kmer_size = cr->ch->kmer_size;

  // remove matching bp at the ends of alleles
  var_trim_alleles(var, var->flank3p);

  var_set_flank_shifts(var);

  // Get kmer keys
  char *key0 = var->key, *key1 = key0 + kmer_size;
  var_set_keys(var, kmer_size);

  // Variant key is defined by <key1><key2>
  // where <key1> is less than <key2>
  // and the keys are the kmer keys of the start of the flanks
  // if <key1> == <key2>, then the variant is `forward` when <kmer1>==<key1>

  int cmp = strncmp(key0, key1, kmer_size);
  char *kmer0 = var->flank5p->buff + var->flank5p->len - kmer_size;
  boolean altered = 0;

  if(cmp > 0 || (cmp == 0 && strcmp(key0, kmer0) != 0))
  {
    var_revcmp(var, cr->ch->num_samples);
    altered = 1;
  }

  // We want our variants to be flush to the right
  if(var->shift_right > 0)
  {
    StrBuf *fl5p = var->flank5p, *fl3p = var->flank3p;

    strbuf_append_strn(fl5p, fl3p->buff, var->shift_right);

    VarBranch *branch;
    for(branch = var->first_allele; branch != NULL; branch = branch->next)
    {
      StrBuf *allele = branch->seq;
      size_t dist = var->shift_right;

      if(allele->len > 0) {
        // DEV: is this correct?
        if(allele->len < var->shift_right) { dist %= allele->len; }
        strbuf_delete(allele, 0, dist);
        strbuf_append_strn(allele, fl3p->buff, dist);
      }
    }

    strbuf_delete(fl3p, 0, var->shift_right);

    var->shift_left += var->shift_right;
    altered = 1;
  }

  if(altered) var_set_keys(var, kmer_size);

  reader_sort_and_link_alleles(cr, var);
}

// Returns 1 on success, 0 on EOF
static char reader_next(CallReader *cr, gzFile fh, Var *var)
{
  StrBuf *line = cr->line;

  // Find the first line
  strbuf_reset(line);
  if(!strbuf_gzreadline_nonempty(line, fh))
    return 0;

  uint32_t i, n, var_num1, var_num2;
  size_t num_samples = cr->ch->num_samples;

  load_assert(line->buff[0] == '>', line);
  n = sscanf(line->buff, ">var_%u_5p_flank", &var_num1);
  load_assert(n == 1, line);
  var->var_num = var_num1;

  load_assert(strbuf_reset_gzreadline(var->flank5p, fh) != 0, line);
  strbuf_chomp(var->flank5p);

  // Skip flank5p covg lines
  for(i = 0; i < num_samples; i++) strbuf_gzskipline(fh);

  load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
  strbuf_chomp(line);
  n = sscanf(line->buff, ">var_%u_3p_flank", &var_num2);
  load_assert(n == 1, line);
  load_assert(var_num1 == var_num2, line);

  load_assert(strbuf_reset_gzreadline(var->flank3p, fh) != 0, line);
  strbuf_chomp(var->flank3p);

  // Skip flank3p covg lines
  for(i = 0; i < num_samples; i++) load_assert(strbuf_gzskipline(fh) > 0, line);

  // Start reading alleles
  load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
  strbuf_chomp(line);

  cr->num_alleles = 0;

  while(line->len > 0 && line->buff[0] == '>')
  {
    if(cr->num_alleles == cr->alleles_cap) {
      cr->alleles_cap <<= 1;
      cr->alleles = realloc(cr->alleles, cr->alleles_cap * sizeof(*cr->alleles));
    }

    VarBranch *branch = branch_new(num_samples);

    uint32_t var_num3, branch_num, branch_nodes;
    n = sscanf(line->buff, ">var_%u_branch_%u length=%u",
               &var_num3, &branch_num, &branch_nodes);

    load_assert(n == 3, line);
    load_assert(var_num3 == var_num1, line);
    load_assert(branch_num == cr->num_alleles, line);

    branch->num_nodes = branch_nodes;

    if(branch_nodes > cr->covgs_cap)
    {
      cr->covgs_cap = ROUNDUP2POW(cr->covgs_cap);
      cr->covgs = realloc(cr->covgs, cr->covgs_cap * sizeof(*cr->alleles));
    }

    // Load sequence
    load_assert(strbuf_reset_gzreadline(branch->seq, fh) != 0, line);
    strbuf_chomp(branch->seq);

    // Read covgs in each sample
    for(i = 0; i < num_samples; i++)
    {
      load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
      strbuf_chomp(line);
      load_assert(parse_uint_liststr(line->buff, cr->covgs, branch_nodes) != 0, line);
      delta_arr_from_uint_arr(cr->covgs, branch_nodes, &(branch->covgs[i]));
    }

    cr->alleles[cr->num_alleles++] = branch;

    // Start reading next allele
    if(strbuf_reset_gzreadline(line, fh) == 0) break;
    strbuf_chomp(line);
  }

  // Expect an empty line between entries
  load_assert(line->len == 0, line);
  load_assert(cr->num_alleles >= 2, line);

  reader_clean_up_var(cr, var);

  return 1;
}

static char reader_next_old_bc(CallReader *cr, gzFile fh, Var *var)
{
  StrBuf *line = cr->line;

  // Find the first line
  strbuf_reset(line);
  if(!strbuf_gzreadline_nonempty(line, fh))
    return 0;

  uint32_t n, var_num1, var_num2;
  size_t num_samples = cr->ch->num_samples;
  cr->num_alleles = 0;

  // Read flank5p
  load_assert(line->buff[0] == '>', line);
  n = sscanf(line->buff, ">var_%u_5p_flank", &var_num1);
  var->var_num = var_num1;
  load_assert(n == 1, line);
  load_assert(strbuf_reset_gzreadline(var->flank5p, fh) != 0, line);
  strbuf_chomp(var->flank5p);

  // Read branches
  strbuf_reset_gzreadline(line, fh);
  strbuf_chomp(line);
  load_assert(line->len > 0, line);
  cr->num_alleles = 0;

  while(1)
  {
    if(cr->num_alleles == cr->alleles_cap) {
      cr->alleles_cap <<= 1;
      cr->alleles = realloc(cr->alleles, cr->alleles_cap * sizeof(*cr->alleles));
    }

    uint32_t var_num3, branch_num, branch_nodes;
    n = sscanf(line->buff, ">var_%u_branch_%u length:%u",
               &var_num3, &branch_num, &branch_nodes);
    // Remove 5p flank node from end of branch
    branch_nodes--;

    if(n != 3) break;

    load_assert(n == 3, line);
    load_assert(var_num3 == var_num1, line);
    load_assert(branch_num == cr->num_alleles+1, line);

    // Load sequence
    VarBranch *branch = branch_new(num_samples);
    branch->num_nodes = branch_nodes;
    load_assert(strbuf_reset_gzreadline(branch->seq, fh) != 0, line);
    strbuf_chomp(branch->seq);
    cr->alleles[cr->num_alleles++] = branch;

    // Read next line
    load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
    strbuf_chomp(line);
  }

  load_assert(cr->num_alleles >= 2, line);

  // Read flank3p
  n = sscanf(line->buff, ">var_%u_3p_flank", &var_num2);
  load_assert(n == 1, line);
  load_assert(var_num1 == var_num2, line);
  load_assert(strbuf_reset_gzreadline(var->flank3p, fh) != 0, line);
  strbuf_chomp(var->flank3p);

  // Read Covg
  size_t brnum, col;
  uint32_t tmp_brnum, tmp_colnum;
  for(brnum = 0; brnum < cr->num_alleles; brnum++)
  {
    VarBranch *branch = cr->alleles[brnum];
    strbuf_reset(line);
    load_assert(strbuf_gzreadline_nonempty(line, fh), line);
    strbuf_chomp(line);
    sscanf(line->buff, "branch%u coverages", &tmp_brnum);
    load_assert(n == 1, line);
    load_assert(tmp_brnum == brnum+1, line);

    for(col = 0; col < num_samples; col++)
    {
      load_assert(strbuf_reset_gzreadline(line, fh), line);
      strbuf_chomp(line);
      n = sscanf(line->buff, "Covg in Colour %u:", &tmp_colnum);
      load_assert(n == 1, line);

      // Read covgs
      load_assert(strbuf_reset_gzreadline(line, fh), line);
      strbuf_chomp(line);

      size_t branch_nodes = branch->num_nodes;
      if(branch_nodes > cr->covgs_cap) {
        cr->covgs_cap = ROUNDUP2POW(cr->covgs_cap);
        cr->covgs = realloc(cr->covgs, cr->covgs_cap * sizeof(*cr->alleles));
      }

      // remove extra whitespace on the end
      while(line->len > 0 && isspace(line->buff[line->len-1])) line->len--;
      line->buff[line->len] = '\0';

      if(branch_nodes == 0) {
        branch->covgs[col].len = 0;
      }
      else {
        // Skip first value
        char *tmp = strchr(line->buff, ' ');
        load_assert(tmp != NULL, line);
        printf("branch_nodes: %zu\n", branch_nodes);

        load_assert(parse_uint_liststr(tmp+1, cr->covgs, branch_nodes) != 0, line);
        delta_arr_from_uint_arr(cr->covgs, branch_nodes, &(branch->covgs[col]));
      }
    }
  }

  // Expect an empty line between entries
  strbuf_reset_gzreadline(line, fh);
  strbuf_chomp(line);
  load_assert(line->len == 0, line);

  reader_clean_up_var(cr, var);

  return 1;
}

// StrBuf tmp is used as temporary memory
static void print_vcf_entry(gzFile out_vcf, gzFile out_flank,
                            Var *var, CallHeader *ch,
                            size_t entry_num, StrBuf *tmp)
{
  size_t i;
  uint32_t kmer_size = ch->kmer_size;
  size_t num_samples = ch->num_samples;
  char is_snp = var_is_snp(var);
  strbuf_reset(tmp);

  char *flank5p = var->flank5p->buff + var->flank5p->len - kmer_size;
  char *flank3p = var->flank3p->buff;

  VarBranch *allele = var->first_allele;

  if(is_snp)
  {
    // SNPs -- only take first base of each allele
    strbuf_append_char(tmp, allele->seq->buff[0]);

    for(allele = allele->next; allele != NULL; allele = allele->next)
    {
      strbuf_append_char(tmp, ',');
      strbuf_append_char(tmp, allele->seq->buff[0]);
    }
  }
  else
  {
    // Use N as ref padding base
    strbuf_append_char(tmp, 'N');
    strbuf_append_str(tmp, allele->seq->buff);

    for(allele = allele->next; allele != NULL; allele = allele->next)
    {
      strbuf_append_str(tmp, ",N");
      strbuf_append_str(tmp, allele->seq->buff);
    }
  }

  gzprintf(out_vcf, "un\t1\tvar%zu\tN\t%s\t.\tPASS\tLF=%.*s;RF=%.*s",
           entry_num, tmp->buff, kmer_size, flank5p, kmer_size, flank3p);

  // Print BN (branch nodes)
  allele = var->first_allele;
  gzprintf(out_vcf, ";BN=%u", allele->num_nodes);
  for(allele = allele->next; allele != NULL; allele = allele->next)
    gzprintf(out_vcf, ",%u", allele->num_nodes);

  gzprintf(out_vcf, ";BUB=%u", var->var_num);

  gzprintf(out_vcf, "\tCOVG");
  for(i = 0; i < num_samples; i++)
  {
    allele = var->first_allele;
    gzputc(out_vcf, '\t');
    delta_arr_gzprint(&(allele->covgs[i]), out_vcf);

    for(allele = allele->next; allele != NULL; allele = allele->next) {
      gzputc(out_vcf, ';');
      delta_arr_gzprint(&(allele->covgs[i]), out_vcf);
    }
  }

  gzputc(out_vcf, '\n');

  // Print 5p flank FASTA
  gzprintf(out_flank, ">var%zu\n%s\n", entry_num, var->flank5p->buff);
}

KHASH_MAP_INIT_STR(vhsh, Var*);

int main(int argc, char **argv)
{
  if(argc != 3) print_usage(usage, NULL);

  // open file to read
  gzFile fh = gzopen(argv[1], "r");
  if(fh == NULL) die("Cannot open input file: %s", argv[1]);

  StrBuf *out_file = strbuf_new();

  strbuf_set(out_file, argv[2]);
  strbuf_append_str(out_file, ".vcf.gz");
  char *vcf_path = strbuf_as_str(out_file);

  gzFile vcf = gzopen(vcf_path, "w");
  if(vcf == NULL) die("Cannot open output file: %s", out_file->buff);

  strbuf_set(out_file, argv[2]);
  strbuf_append_str(out_file, ".5pflanks.fa.gz");
  char *flanks_path = strbuf_as_str(out_file);

  gzFile flankfh = gzopen(flanks_path, "w");
  if(flankfh == NULL) die("Cannot open output file: %s", out_file->buff);

  message("proc_calls\n  reading: %s\n  vcf: %s\n  fasta: %s\n\n",
          argv[1], vcf_path, flanks_path);

  // Deduce input filetype
  StrBuf *tmp = strbuf_new();
  strbuf_gzreadline(tmp, fh);

  if(gzseek(fh, 0, SEEK_SET) != 0)
    die("Read error on input file: %s", argv[1]);

  CallHeader *ch = header_new();

  char (*next_var)(CallReader *cr, gzFile fh, Var *var);

  if(strncasecmp(tmp->buff,"##fileformat=", strlen("##fileformat=")) == 0) {
    parse_shaded_header(fh, ch);
    next_var = reader_next;
  }
  else {
    synthesize_bubble_caller_header(fh, ch);
    next_var = reader_next_old_bc;
    if(gzseek(fh, 0, SEEK_SET) != 0)
      die("Read error on input file: %s", argv[1]);
  }

  print_vcf_header(vcf, ch, argc, argv);

  CallReader *cr = reader_new(ch);
  khash_t(vhsh) *varhash = kh_init(vhsh);
  khiter_t k;
  int hret;

  size_t bubbles_read = 0, entries_printed = 0;

  // read, add to hash
  Var *var = var_new(ch->kmer_size);

  while(next_var(cr, fh, var))
  {
    k = kh_put(vhsh, varhash, var->key, &hret);
    if(hret == 0) {
      var_merge(kh_value(varhash, k), var, ch->num_samples);
    } else {
      kh_value(varhash, k) = var;
      var = var_new(ch->kmer_size);
    }
    bubbles_read++;
  }

  var_free(var, ch->num_samples);

  // iterate hash, print variants

  for (k = kh_begin(varhash); k != kh_end(varhash); ++k) {
    if(kh_exist(varhash, k)) {
      var = kh_value(varhash, k);
      print_vcf_entry(vcf, flankfh, var, ch, entries_printed, tmp);
      entries_printed++;
      var_free(var, ch->num_samples);
    }
  }

  printf("  %zu bubbles loaded\n", bubbles_read);
  printf("  %zu vcf entries printed\n\n", entries_printed);

  strbuf_free(tmp);
  kh_destroy(vhsh, varhash);
  reader_free(cr);
  header_free(ch);

  gzclose(vcf);
  gzclose(flankfh);
  gzclose(fh);

  strbuf_set(out_file, argv[2]);
  strbuf_append_str(out_file, ".5pflanks.sam");

  printf("  Next steps - map with stampy and align to ref:\n");
  printf("    stampy.py -G hg19 hg19.fa\n");
  printf("    stampy.py -g hg19 -H hg19\n");
  printf("    stampy.py -g hg19 -h hg19 --inputformat=fasta -M %s > %s\n\n",
         flanks_path, out_file->buff);
  printf("    ./place_calls %s %s hg19.fa\n\n", vcf_path, out_file->buff);

  free(vcf_path);
  free(flanks_path);
  strbuf_free(out_file);

  return 0;
}
