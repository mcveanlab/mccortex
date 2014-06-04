#include "global.h"

#include <time.h>
#include <ctype.h> // isspace

#include "khash.h"

#include "commands.h"
#include "dna.h"
#include "util.h"
#include "file_util.h"

#define PROC_MAXK 100

const char unique_usage[] =
"usage: "CMD" unique [options] <path_calls.bubbles> <out.base>\n"
"  Produces files <out.base>.vcf and <out.base>.5pflank.fa\n";

#define load_assert(cond,sb) \
  if(!(cond)) { die("Loading err ["QUOTE_MACRO(cond)"]: '%s'", (sb)->buff); }

/* alleles */

typedef struct VarBranch VarBranch;

struct VarBranch {
  StrBuf seq;
  uint32_t num_nodes;
  VarBranch *next;
};

typedef struct
{
  char key[2*PROC_MAXK+1]; // concat of key0, key1 (keys of kmer0, kmer1)
  StrBuf flank5p, flank3p;
  VarBranch *first_allele; // Linkedlist
  StrBuf name;
  uint32_t shift_left, shift_right;
} Var;

static VarBranch* branch_new()
{
  VarBranch *branch = ctx_malloc(sizeof(VarBranch));
  strbuf_alloc(&branch->seq, 512);
  branch->num_nodes = 0;
  branch->next = NULL;
  return branch;
}

static void branch_free(VarBranch *branch)
{
  strbuf_dealloc(&branch->seq);
  ctx_free(branch);
}

/* Variants */

static Var* var_new(size_t kmer_size)
{
  Var *var = ctx_malloc(sizeof(Var));
  var->key[0] = var->key[2*kmer_size] = 0;
  strbuf_alloc(&var->flank5p, 512);
  strbuf_alloc(&var->flank3p, 512);
  strbuf_alloc(&var->name, 128);
  var->first_allele = NULL;
  var->shift_left = var->shift_right = 0;
  return var;
}

static void var_free(Var *var)
{
  strbuf_dealloc(&var->flank5p);
  strbuf_dealloc(&var->flank3p);
  strbuf_dealloc(&var->name);

  VarBranch *allele, *next;
  for(allele = var->first_allele; allele != NULL; allele = next)
  {
    next = allele->next;
    branch_free(allele);
  }
  ctx_free(var);
}

static void var_revcmp(Var *var)
{
  dna_reverse_complement_str(var->flank5p.buff, var->flank5p.len);
  dna_reverse_complement_str(var->flank3p.buff, var->flank3p.len);

  SWAP(var->flank5p, var->flank3p);
  SWAP(var->shift_left, var->shift_right);

  VarBranch *allele;
  for(allele = var->first_allele; allele != NULL; allele = allele->next) {
    dna_reverse_complement_str(allele->seq.buff, allele->seq.len);
  }
}

// Merge var1 into var0, free var1
static void var_merge(Var *var0, Var *var1)
{
  // Compare sorted linked lists of alleles to count number of novel alleles
  VarBranch *allele0 = var0->first_allele;
  VarBranch *allele1 = var1->first_allele;

  VarBranch null_allele;
  null_allele.next = NULL;

  VarBranch *prev = &null_allele;

  while(allele0 != NULL && allele1 != NULL)
  {
    int cmp = strcmp(allele0->seq.buff, allele1->seq.buff);

    if(cmp <= 0) {
      prev->next = allele0;
      allele0 = allele0->next;

      if(cmp == 0)
      {
        VarBranch *tmp = allele1->next;
        branch_free(allele1);
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

  strbuf_append_char(&var0->name, ',');
  strbuf_append_strn(&var0->name, var1->name.buff, var1->name.len);
}

// Trim from right side, add to flank3p
static void var_trim_alleles(Var *var, StrBuf *flank3p)
{
  VarBranch *branch0 = var->first_allele, *branch1 = branch0->next, *branch;
  size_t trim;
  size_t len0 = branch0->seq.len;
  const char *allele0 = branch0->seq.buff;
  bool match = 1;

  for(trim = 0; trim < len0; trim++)
  {
    char lastc = allele0[len0-1-trim];

    // loop over alleles the other alleles
    for(branch = branch1; branch != NULL; branch = branch->next)
    {
      if(trim >= branch->seq.len ||
         branch->seq.buff[branch->seq.len-1-trim] != lastc)
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
      strbuf_shrink(&branch->seq, branch->seq.len - trim);
  }
}

static void var_set_flank_shifts(Var *var)
{
  // May need to left shift 3p start
  size_t min_left = UINT_MAX, min_right = UINT_MAX;
  size_t j, k, max_dist;
  VarBranch *branch;
  const StrBuf *fl5p = &var->flank5p;
  const StrBuf *fl3p = &var->flank3p;

  for(branch = var->first_allele;
      branch != NULL && min_left + min_right > 0;
      branch = branch->next)
  {
    StrBuf *allele = &branch->seq;

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

  var->shift_left = (uint32_t)min_left;
  var->shift_right = (uint32_t)min_right;
}

// Warning: does not null terminate
static void str_get_key(char *key, const char *kmer, size_t kmer_size)
{
  memcpy(key, kmer, kmer_size*sizeof(char));
  dna_reverse_complement_str(key, kmer_size);

  if(strncmp(kmer, key, kmer_size) < 0)
    memcpy(key, kmer, kmer_size*sizeof(char));
}

static void var_set_keys(Var *var, size_t kmer_size)
{
  StrBuf *fl5p = &var->flank5p, *fl3p = &var->flank3p;

  #ifdef CTXVERBOSE
    message(" 5p:%s\n 3p:%s\n", fl5p->buff, fl3p->buff);
  #endif

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

  return strcmp(a0->seq.buff, a1->seq.buff);
}

// Only applies if alleles have been trimmed of padding bases
static char var_is_snp(Var *var)
{
  VarBranch *allele;
  for(allele = var->first_allele; allele != NULL; allele = allele->next)
    if(allele->seq.len != 1) return 0;
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
  CallHeader *ch = ctx_malloc(sizeof(CallHeader));
  ch->num_samples = ch->kmer_size = 0;
  ch->sample_names = NULL;
  ch->hcap = 128;
  ch->hlines = 0;
  ch->tags = ctx_malloc(ch->hcap * sizeof(*ch->tags));
  ch->values = ctx_malloc(ch->hcap * sizeof(*ch->tags));
  ch->is_old_bc = 0;
  return ch;
}

static void header_free(CallHeader *ch)
{
  size_t i;
  for(i = 0; i < ch->hlines; i++) { ctx_free(ch->tags[i]); ctx_free(ch->values[i]); }
  for(i = 0; i < ch->num_samples; i++) { ctx_free(ch->sample_names[i]); }
  ctx_free(ch->tags);
  ctx_free(ch->values);
  ctx_free(ch);
}

static void header_add(CallHeader *ch, char *tag, char *value)
{
  if(ch->hlines == ch->hcap) {
    ch->hcap <<= 1;
    ch->tags = ctx_realloc(ch->tags, ch->hcap * sizeof(char*));
    ch->values = ctx_realloc(ch->values, ch->hcap * sizeof(char*));
  }

  ch->tags[ch->hlines] = tag;
  ch->values[ch->hlines] = value;
  ch->hlines++;
}

static void print_vcf_header(gzFile vcf, CallHeader *ch)
{
  size_t i;
  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  gzprintf(vcf, "##fileformat=VCFv4.1\n");
  gzprintf(vcf, "##fileDate=%s\n", datestr);
  gzprintf(vcf, "##reference=unplaced\n");
  gzprintf(vcf, "##phasing=none\n");
  gzprintf(vcf, "##procCmd=%s\n", cmd_get_cmdline());
  gzprintf(vcf, "##procCwd=%s\n", cmd_get_cwd());

  gzprintf(vcf, "##procDate=%s\n", datestr);

  for(i = 0; i < ch->hlines; i++) {
    if(strcasecmp(ch->tags[i],"fileformat") != 0)
      gzprintf(vcf, "##%s=%s\n", ch->tags[i], ch->values[i]);
  }

  gzprintf(vcf, "##INFO=<ID=LF,Number=1,Type=String,Description=\"Left flank\">\n");
  gzprintf(vcf, "##INFO=<ID=RF,Number=1,Type=String,Description=\"Right flank\">\n");
  gzprintf(vcf, "##INFO=<ID=BN,Number=.,Type=Integer,Description=\"Branch nodes; length in kmers of the alleles\">\n");
  gzprintf(vcf, "##INFO=<ID=BUB,Number=1,Type=String,Description=\"Cortex bubble index\">\n");
  gzprintf(vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
  gzprintf(vcf, "##contig=<ID=un,length=1000000,assembly=None>\n");

  gzprintf(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  // for(i = 0; i < ch->num_samples; i++) {
  //   if(strcmp(ch->sample_names[i], "undefined") == 0)
  //     gzprintf(vcf, "\tsample%zu", i);
  //   else
  //     gzprintf(vcf, "\t%s", ch->sample_names[i]);
  // }

  gzprintf(vcf, "\n");
}

static void parse_bubble_header(gzFile in, CallHeader* ch)
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
    strbuf_append_char(line, (char)c);
    strbuf_gzreadline(line, in);
    strbuf_chomp(line);

    char *eq = strchr(line->buff, '=');
    if(line->buff[1] != '#' || eq == NULL) die("parse error: %s", line->buff);

    *eq = 0;
    char *tag = strdup(line->buff+2);
    char *value = strdup(eq+1);
    *eq = '=';

    if(strcasecmp(tag, "ctxKmerSize") == 0)
    {
      // Get kmer size from cmd line
      unsigned int kmer_size;
      if(!parse_entire_uint(value, &kmer_size))
        die("Couldn't read kmer size line: %s=%s", tag, value);
      ch->kmer_size = kmer_size;
    }
    else if(strcasecmp(tag, "ctxNumColoursUsedInCalling") == 0)
    {
      if(ch->num_samples != 0) die("duplicate header line: %s", line->buff);

      size_t num_samples;
      if(!parse_entire_size(value, &num_samples))
        die("parse error on line: %s", line->buff);
      ch->num_samples = num_samples;

      ch->sample_names = ctx_malloc(ch->num_samples*sizeof(char*));
    }
    else if(strncasecmp(tag, "colour", 6) == 0)
    {
      if(ch->sample_names == NULL) die("missing 'ctxNumColoursUsedInCalling=' line");

      strbuf_ensure_capacity(tmp, line->len);
      if(sscanf(line->buff, "##colour=<ID=%[^,],", tmp->buff) != 1)
        die("parse error on line: %s", line->buff);

      if(sampleid >= ch->num_samples)
        die("more samples than expected in header [%zu > %zu]", sampleid, ch->num_samples);

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

#define header_assert(cond,sb) \
  if(!(cond)) { die("Loading err ["QUOTE_MACRO(cond)"]: %s", (sb)->buff); }

// >[NAME][SUFFIX]
// Returns length of var name
// static size_t get_var_name(const char *line, size_t line_len,
//                            const char *suffix, char *name)
// {
//   const size_t suffix_len = strlen(suffix);
//   const char *match = line+line_len-suffix_len;
//   if(line[0] != '>') die("Expected '>': %s", line);
//   if(line_len < suffix_len+2 || strcmp(match,suffix) != 0)
//     die("Cannot find var name [line: %s]", line);
//   size_t name_len = line_len - suffix_len - 1;
//   memcpy(name, line+1, name_len);
//   name[name_len] = '\0';
//   return name_len;
// }

// This function is the same as above but does not expect the suffix to
// be at the end of the line.  This is needed for handling old output
// where meta data is on the name line
// >[NAME][SUFFIX]
// Returns length of var name
static size_t get_var_name(const char *line, size_t line_len,
                           const char *suffix, char *name)
{
  (void)line_len;
  if(line[0] != '>') die("Expected '>': %s", line);
  const char *match = strstr(line+1, suffix);
  if(match == NULL) die("Cannot find var name [line: %s]", line);
  size_t name_len = (size_t)(match-(line+1));
  memcpy(name, line+1, name_len);
  name[name_len] = '\0';
  return name_len;
}

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

  char name[line->len];
  size_t name_len = get_var_name(line->buff, line->len, "_5p_flank", name);

  int flank5plen, kmer_size;
  int n = sscanf(line->buff+1+name_len, "_5p_flank length:%i INFO:KMER=%i",
                 &flank5plen, &kmer_size);

  header_assert(n == 2, line);
  header_assert(kmer_size & 1, line);
  ch->kmer_size = (size_t)kmer_size;

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

  ch->sample_names = ctx_malloc(ch->num_samples * sizeof(char*));
  for(i = 0; i < ch->num_samples; i++)
  {
    strbuf_reset(line);
    strbuf_sprintf(line, "sample%zu", i);
    ch->sample_names[i] = strbuf_dup(line);
    strbuf_reset(line);
    strbuf_sprintf(line, "<ID=sample%zu,name=\"undefined\">", i);
    header_add(ch, strdup("colour"), strbuf_dup(line));
  }

  strbuf_free(line);
}


/* CallReader */

typedef struct
{
  CallHeader *ch;
  // Temp variables used for loading
  StrBuf *line;
  VarBranch **alleles;
  uint32_t num_alleles, alleles_cap;
} CallReader;

static CallReader* reader_new(CallHeader *ch)
{
  CallReader *cr = ctx_malloc(sizeof(CallReader));
  cr->ch = ch;
  cr->line = strbuf_new();
  cr->num_alleles = 0;
  cr->alleles_cap = 32;
  cr->alleles = ctx_malloc(cr->alleles_cap * sizeof(*cr->alleles));
  return cr;
}

static void reader_free(CallReader *cr)
{
  strbuf_free(cr->line);
  ctx_free(cr->alleles);
  ctx_free(cr);
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

  #ifdef CTXVERBOSE
    printf("\n\n");
    printf("%s\n 5p:%s\n 3p:%s\n", var->name.buff,
           var->flank5p.buff, var->flank3p.buff);
    size_t i;
    for(i = 0; i < cr->num_alleles; i++)
      printf(" %2zu:%s\n", i, cr->alleles[i]->seq.buff);
  #endif

  size_t kmer_size = cr->ch->kmer_size;

  // remove matching bp at the ends of alleles
  var_trim_alleles(var, &var->flank3p);

  var_set_flank_shifts(var);

  // Get kmer keys
  char *key0 = var->key, *key1 = key0 + kmer_size;
  var_set_keys(var, kmer_size);

  // Variant key is defined by <key1><key2>
  // where <key1> is less than <key2>
  // and the keys are the kmer keys of the start of the flanks
  // if <key1> == <key2>, then the variant is `forward` when <kmer1>==<key1>

  int cmp = strncmp(key0, key1, kmer_size);
  char *kmer0 = var->flank5p.buff + var->flank5p.len - kmer_size;
  bool altered = 0;

  if(cmp > 0 || (cmp == 0 && strcmp(key0, kmer0) != 0))
  {
    var_revcmp(var);
    altered = 1;
  }

  // We want our variants to be flush to the right
  if(var->shift_right > 0)
  {
    StrBuf *fl5p = &var->flank5p, *fl3p = &var->flank3p;

    strbuf_append_strn(fl5p, fl3p->buff, var->shift_right);

    VarBranch *branch;
    for(branch = var->first_allele; branch != NULL; branch = branch->next)
    {
      StrBuf *allele = &branch->seq;
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

  // Get var name
  StrBuf *name = &var->name;
  strbuf_ensure_capacity(name, line->len);
  name->len = get_var_name(line->buff, line->len, "_5p_flank", name->buff);

  load_assert(strbuf_reset_gzreadline(&var->flank5p, fh) != 0, line);
  strbuf_chomp(&var->flank5p);

  load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
  strbuf_chomp(line);

  // Check var name
  load_assert(line->buff[0] == '>', line);
  load_assert(strncmp(line->buff+1, name->buff, name->len) == 0, line);
  load_assert(strncmp(line->buff+1+name->len, "_3p_flank", 9) == 0, line);

  load_assert(strbuf_reset_gzreadline(&var->flank3p, fh) != 0, line);
  strbuf_chomp(&var->flank3p);

  // Start reading alleles
  load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
  strbuf_chomp(line);

  cr->num_alleles = 0;

  while(line->len > 0 && line->buff[0] == '>')
  {
    if(cr->num_alleles == cr->alleles_cap) {
      cr->alleles_cap <<= 1;
      cr->alleles = ctx_realloc(cr->alleles, cr->alleles_cap * sizeof(*cr->alleles));
    }

    VarBranch *branch = branch_new();

    uint32_t branch_num, branch_nodes;

    // Check name
    load_assert(line->buff[0] == '>', line);
    load_assert(strncmp(line->buff+1, name->buff, name->len) == 0, line);

    int n = sscanf(line->buff+1+name->len, "_branch_%u kmers=%u",
                   &branch_num, &branch_nodes);

    load_assert(n == 2, line);
    load_assert(branch_num == cr->num_alleles, line);

    branch->num_nodes = branch_nodes;

    // Load sequence
    load_assert(strbuf_reset_gzreadline(&branch->seq, fh) != 0, line);
    strbuf_chomp(&branch->seq);

    cr->alleles[cr->num_alleles++] = branch;

    // Start reading next allele
    if(strbuf_reset_gzreadline(line, fh) == 0) break;
    strbuf_chomp(line);
  }

  if(cr->num_alleles < 2) {
    printf("num_alleles: %zu\n", (size_t)cr->num_alleles);
    printf("var: %s\n", var->name.buff);
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

  int n;
  size_t num_samples = cr->ch->num_samples;
  cr->num_alleles = 0;
  strbuf_ensure_capacity(&var->name, line->len);

  // Get var name
  StrBuf *name = &var->name;
  strbuf_ensure_capacity(name, line->len);
  name->len = get_var_name(line->buff, line->len, "_5p_flank", name->buff);
  load_assert(strncmp(line->buff+1+name->len, "_5p_flank", 9) == 0, line);

  load_assert(strbuf_reset_gzreadline(&var->flank5p, fh) != 0, line);
  strbuf_chomp(&var->flank5p);

  // Read branches
  strbuf_reset_gzreadline(line, fh);
  strbuf_chomp(line);
  load_assert(line->len > 0, line);
  cr->num_alleles = 0;

  while(1)
  {
    if(cr->num_alleles == cr->alleles_cap) {
      cr->alleles_cap <<= 1;
      cr->alleles = ctx_realloc(cr->alleles, cr->alleles_cap * sizeof(*cr->alleles));
    }

    uint32_t branch_num, branch_nodes;

    load_assert(line->buff[0] == '>', line);
    load_assert(strncmp(line->buff+1, name->buff, name->len) == 0, line);
    n = sscanf(line->buff+1+name->len, "_branch_%u length:%u",
               &branch_num, &branch_nodes);

    if(n != 2) break;

    // Remove 5p flank node from end of branch
    branch_nodes--;
    load_assert(branch_num == cr->num_alleles+1, line);

    // Load sequence
    VarBranch *branch = branch_new();
    branch->num_nodes = branch_nodes;
    load_assert(strbuf_reset_gzreadline(&branch->seq, fh) != 0, line);
    strbuf_chomp(&branch->seq);
    cr->alleles[cr->num_alleles++] = branch;

    // Read next line
    load_assert(strbuf_reset_gzreadline(line, fh) != 0, line);
    strbuf_chomp(line);
  }

  load_assert(cr->num_alleles >= 2, line);

  // Read >name_3p_flank
  load_assert(line->buff[0] == '>', line);
  load_assert(strncmp(line->buff+1, name->buff, name->len) == 0, line);
  load_assert(strncmp(line->buff+1+name->len, "_3p_flank", 9) == 0, line);

  load_assert(strbuf_reset_gzreadline(&var->flank3p, fh) != 0, line);
  strbuf_chomp(&var->flank3p);

  // Read Covg
  size_t brnum, col;
  uint32_t tmp_brnum, tmp_colnum;
  for(brnum = 0; brnum < cr->num_alleles; brnum++)
  {
    // VarBranch *branch = cr->alleles[brnum];
    strbuf_reset(line);
    load_assert(strbuf_gzreadline_nonempty(line, fh), line);
    strbuf_chomp(line);
    n = sscanf(line->buff, "branch%u coverages", &tmp_brnum);
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

      // remove extra whitespace on the end
      while(line->len > 0 && isspace(line->buff[line->len-1])) line->len--;
      line->buff[line->len] = '\0';
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
  // size_t i;
  const size_t kmer_size = ch->kmer_size;
  // size_t num_samples = ch->num_samples;
  char is_snp = var_is_snp(var);
  strbuf_reset(tmp);

  char *flank5p = var->flank5p.buff + var->flank5p.len - kmer_size;
  char *flank3p = var->flank3p.buff;

  VarBranch *allele = var->first_allele;

  if(is_snp)
  {
    // SNPs -- only take first base of each allele
    strbuf_append_char(tmp, allele->seq.buff[0]);

    for(allele = allele->next; allele != NULL; allele = allele->next)
    {
      strbuf_append_char(tmp, ',');
      strbuf_append_char(tmp, allele->seq.buff[0]);
    }
  }
  else
  {
    // Use N as ref padding base
    strbuf_append_char(tmp, 'N');
    strbuf_append_strn(tmp, allele->seq.buff, allele->seq.len);

    for(allele = allele->next; allele != NULL; allele = allele->next)
    {
      strbuf_append_str(tmp, ",N");
      strbuf_append_strn(tmp, allele->seq.buff, allele->seq.len);
    }
  }

  gzprintf(out_vcf, "un\t1\tvar%zu\tN\t", entry_num);
  gzputs(out_vcf, tmp->buff);
  gzprintf(out_vcf, "\t.\tPASS\tLF=%.*s;RF=%.*s",
           (int)kmer_size, flank5p, (int)kmer_size, flank3p);

  // Print BN (branch nodes)
  allele = var->first_allele;
  gzprintf(out_vcf, ";BN=%u", allele->num_nodes);
  for(allele = allele->next; allele != NULL; allele = allele->next)
    gzprintf(out_vcf, ",%u", allele->num_nodes);

  gzprintf(out_vcf, ";BUB=%s", var->name.buff);

  gzprintf(out_vcf, "\t.");
  // for(i = 0; i < num_samples; i++) {
  //   gzprintf(out_vcf, "\t.");
  // }

  gzputc(out_vcf, '\n');

  // Print 5p flank FASTA
  gzprintf(out_flank, ">var%zu\n%s\n", entry_num, var->flank5p.buff);
}

KHASH_MAP_INIT_STR(vhsh, Var*);

int ctx_unique(CmdArgs *args)
{
  // hide unused function warnings
  (void)kh_clear_vhsh;
  (void)kh_get_vhsh;
  (void)kh_del_vhsh;

  char **argv = args->argv;
  // Have already checked we have exactly 2 arguments

  const char *input_path = argv[0];
  const char *output_path = argv[1];

  // open file to read
  gzFile fh = gzopen(input_path, "r");
  if(fh == NULL) die("Cannot open input file: %s", input_path);

  StrBuf *out_file = strbuf_new();

  strbuf_set(out_file, output_path);
  strbuf_append_str(out_file, ".vcf.gz");
  char *vcf_path = strbuf_as_str(out_file);

  gzFile vcf = gzopen(vcf_path, "w");
  if(vcf == NULL) die("Cannot open output file: %s", out_file->buff);

  strbuf_set(out_file, output_path);
  strbuf_append_str(out_file, ".5pflanks.fa.gz");
  char *flanks_path = strbuf_as_str(out_file);

  gzFile flankfh = gzopen(flanks_path, "w");
  if(flankfh == NULL) die("Cannot open output file: %s", out_file->buff);

  status("reading: %s\n", input_path);
  status("vcf: %s\n", vcf_path);
  status("fasta: %s\n", flanks_path);

  // Deduce input filetype
  StrBuf *tmp = strbuf_new();
  strbuf_gzreadline(tmp, fh);

  if(gzseek(fh, 0, SEEK_SET) != 0)
    die("Read error on input file: %s", input_path);

  CallHeader *ch = header_new();

  char (*next_var)(CallReader *cr, gzFile fh, Var *var);

  if(strncasecmp(tmp->buff,"##fileformat=", strlen("##fileformat=")) == 0) {
    parse_bubble_header(fh, ch);
    next_var = reader_next;
  }
  else {
    synthesize_bubble_caller_header(fh, ch);
    next_var = reader_next_old_bc;
    if(gzseek(fh, 0, SEEK_SET) != 0)
      die("Read error on input file: %s", input_path);
  }

  print_vcf_header(vcf, ch);

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
      var_merge(kh_value(varhash, k), var);
    } else {
      kh_value(varhash, k) = var;
      var = var_new(ch->kmer_size);
    }
    bubbles_read++;
  }

  var_free(var);

  // iterate hash, print variants

  for (k = kh_begin(varhash); k != kh_end(varhash); ++k) {
    if(kh_exist(varhash, k)) {
      var = kh_value(varhash, k);
      print_vcf_entry(vcf, flankfh, var, ch, entries_printed, tmp);
      entries_printed++;
      var_free(var);
    }
  }

  char nbubbles_str[100], nentries_str[100];
  ulong_to_str(bubbles_read, nbubbles_str);
  ulong_to_str(entries_printed, nentries_str);
  status("%s bubbles loaded\n", nbubbles_str);
  status("%s vcf entries printed\n", nentries_str);

  strbuf_free(tmp);
  kh_destroy(vhsh, varhash);
  reader_free(cr);
  header_free(ch);

  gzclose(vcf);
  gzclose(flankfh);
  gzclose(fh);

  strbuf_set(out_file, output_path);
  strbuf_append_str(out_file, ".5pflanks.sam");

  // message(...) is same as status(...) but without timestamp
  status("Next steps - map variants to a reference (e.g. with stampy):\n");
  message("  stampy.py -G hg19 hg19.fa\n");
  message("  stampy.py -g hg19 -H hg19\n");
  message("  stampy.py -g hg19 -h hg19 --inputformat=fasta -M %s > %s\n",
         flanks_path, out_file->buff);
  message("  "CMD" place %s %s hg19.fa\n", vcf_path, out_file->buff);

  ctx_free(vcf_path);
  ctx_free(flanks_path);
  strbuf_free(out_file);

  return EXIT_SUCCESS;
}
