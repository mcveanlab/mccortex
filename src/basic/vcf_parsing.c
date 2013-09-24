#include "global.h"
#include "vcf_parsing.h"
#include "util.h"
#include "binary_kmer.h"

static void strbuf_arr_resize(StrBuf **arr, size_t *cap, size_t newcap)
{
  size_t i;
  newcap = ROUNDUP2POW(newcap);
  *arr = realloc2(*arr, sizeof(StrBuf) * newcap);
  for(i = *cap; i < newcap; i++) strbuf_alloc(&(*arr)[i], 64);
  *cap = newcap;
}

// DP=asdf;TXT="a;b;c";F1;X=4
// ends:  ^           ^  ^   ^
char* info_tag_end(char *str)
{
  boolean speechmrks = false;
  for(; *str; str++) {
    if(*str == '"') speechmrks = !speechmrks;
    else if(*str == ';' && !speechmrks) return str;
  }
  return NULL;
}

StrBuf* info_tag_find(vcf_entry_t *entry, const char* str)
{
  size_t i, len = strlen(str);
  for(i = 0; i < entry->num_info; i++) {
    const char *buff = entry->info[i].buff;
    if(strcasecmp(buff, str) == 0 ||
       (buff[len] == '=' && strncasecmp(buff, str, len) == 0))
    {
      return &entry->info[i];
    }
  }
  return NULL;
}

void info_tag_add(vcf_entry_t *entry, const char* fmt, ...)
{
  if(entry->num_info+1 >= entry->info_capacity)
    strbuf_arr_resize(&entry->info, &entry->info_capacity, entry->num_info+1);

  StrBuf* sbuf = &entry->info[entry->num_info++];

  va_list argptr;
  va_start(argptr, fmt);
  strbuf_reset(sbuf);
  strbuf_vsprintf(sbuf, 0, fmt, argptr);
  va_end(argptr);
}

void vcf_entry_alloc(vcf_entry_t *entry, uint32_t num_samples)
{
  size_t i, j, num_cols = VCFSAMPLES+num_samples;
  entry->alts_capacity = entry->info_capacity = 2;
  entry->num_info = entry->num_alts = 0;
  entry->cols = malloc2(sizeof(StrBuf) * num_cols);
  entry->alts = malloc2(sizeof(StrBuf) * entry->alts_capacity);
  entry->info = malloc2(sizeof(StrBuf) * entry->info_capacity);

  for(i = 0; i < num_cols; i++)
    strbuf_alloc(&entry->cols[i], 1024);

  for(i = 0; i < entry->alts_capacity; i++) {
    strbuf_alloc(&entry->alts[i], 64);
    strbuf_alloc(&entry->info[i], 1024);
  }

  // samples
  entry->covgs = malloc2(num_samples * sizeof(DeltaArray*));

  for(i = 0; i < num_samples; i++)
  {
    entry->covgs[i] = malloc2(entry->alts_capacity * sizeof(DeltaArray));
    for(j = 0; j < entry->alts_capacity; j++)
      delta_arr_alloc(&entry->covgs[i][j]);
  }
}

void vcf_entry_dealloc(vcf_entry_t *entry, uint32_t num_samples)
{
  size_t i, j, num_cols = VCFSAMPLES+num_samples;
  for(i = 0; i < num_cols; i++) strbuf_dealloc(&entry->cols[i]);
  for(i = 0; i < entry->alts_capacity; i++) strbuf_dealloc(&entry->alts[i]);
  for(i = 0; i < entry->info_capacity; i++) strbuf_dealloc(&entry->info[i]);
  for(i = 0; i < num_samples; i++) {
    for(j = 0; j < entry->alts_capacity; j++)
      delta_arr_dealloc(&entry->covgs[i][j]);
    free(entry->covgs[i]);
  }
  free(entry->cols);
  free(entry->alts);
  free(entry->info);
  free(entry->covgs);
}

void vcf_entry_parse(StrBuf *line, vcf_entry_t *entry, uint32_t num_samples)
{
  strbuf_chomp(line);
  char *end, *tmp;
  StrBuf *buf;

  // Split vcf line into fields
  // size_t num_cols = count_char(line->buff, '\t') + 1;

  size_t num_cols = 0;
  tmp = line->buff;
  while(1)
  {
    if(num_cols == VCFSAMPLES+num_samples) die("Too many VCF columns");
    end = strchr(tmp, '\t');
    if(end != NULL) *end = '\0';
    strbuf_set(&entry->cols[num_cols++], tmp);
    if(end != NULL) *end = '\t';
    else break;
    tmp = end + 1;
  }

  if(num_cols != VCFSAMPLES+num_samples)
    die("Incorrect number of VCF columns: '%s'", line->buff);

  // Split ALT alleles
  size_t alts_count = count_char(entry->cols[VCFALT].buff, ',') + 1;
  if(alts_count > entry->alts_capacity)
  {
    // Increase in the number of alleles
    size_t old_alts_cap = entry->alts_capacity;
    strbuf_arr_resize(&entry->alts, &entry->alts_capacity, alts_count);

    size_t i, j, covgsize = entry->alts_capacity * sizeof(DeltaArray);
    for(i = 0; i < num_samples; i++)
    {
      entry->covgs[i] = realloc2(entry->covgs[i], covgsize);
      for(j = old_alts_cap; j < entry->alts_capacity; j++)
        delta_arr_alloc(&(entry->covgs[i][j]));
    }
  }

  entry->num_alts = 0;
  tmp = entry->cols[VCFALT].buff;
  while(1)
  {
    end = strchr(tmp, ',');
    if(end != NULL) *end = '\0';
    while(*tmp == 'N') tmp++; // Skip Ns at the beginning of alleles
    strbuf_set(&entry->alts[entry->num_alts++], tmp);
    if(end != NULL) *end = ',';
    else break;
    tmp = end + 1;
  }

  // Split INFO alleles
  entry->num_info = 0;
  tmp = entry->cols[VCFINFO].buff;
  entry->lf = entry->rf = NULL;

  size_t info_count = 0;
  while((tmp = info_tag_end(tmp)) != NULL) { info_count++; tmp++; }
  tmp = entry->cols[VCFINFO].buff;

  if(info_count > entry->info_capacity)
    strbuf_arr_resize(&entry->info, &entry->info_capacity, info_count);

  while(1)
  {
    end = info_tag_end(tmp);
    if(end != NULL) *end = '\0';
    buf = &entry->info[entry->num_info++];
    strbuf_set(buf, tmp);
    if(!strncmp(tmp, "LF=", 3)) entry->lf = buf;
    else if(!strncmp(tmp, "RF=", 3)) entry->rf = buf;
    if(end != NULL) *end = ';';
    else break;
    tmp = end + 1;
  }

  // Load sample covg info

  size_t i, j;
  for(i = 0; i < num_samples; i++)
  {
    tmp = entry->cols[VCFSAMPLES+i].buff;
    size_t count = count_char(tmp, ';') + 1;
    if(count != entry->num_alts)
      die("Invalid GT: %s [%zu vs %zu]", tmp, count, entry->num_alts);

    for(j = 0; j < entry->num_alts; j++)
    {
      if((end = strchr(tmp, ';')) != NULL) *end = '\0';
      delta_arr_from_str(tmp, &(entry->covgs[i][j]));
      delta_array_unpack(&(entry->covgs[i][j]));
      if(end != NULL) *end = ';';
      tmp = end+1;
    }
  }
}

void vcf_entry_revcmp(vcf_entry_t *entry)
{
  // For debugging
  info_tag_add(entry, "BUBREV");

  size_t i;
  for(i = 0; i < entry->num_alts; i++)
    reverse_complement_str(entry->alts[i].buff, entry->alts[i].len);

  // reverse complement and swap lflank and rflank
  // printf("lf:'%s'; rf:'%s'\n", entry->lf.buff, entry->rf.buff);
  reverse_complement_str(entry->lf->buff+3, entry->lf->len-3);
  reverse_complement_str(entry->rf->buff+3, entry->rf->len-3);
  entry->lf->buff[0] = 'R';
  entry->rf->buff[0] = 'L';
  StrBuf *tmpbuf; SWAP(entry->lf, entry->rf, tmpbuf);
}

size_t vcf_entry_longest_allele(vcf_entry_t *entry)
{
  size_t i, max = 0;
  for(i = 0; i < entry->num_alts; i++) {
    max = MAX2(max, entry->alts[i].len);
  }
  return max;
}

void vcf_entry_print(const vcf_entry_t *entry, FILE *out, uint32_t num_samples)
{
  size_t i, num_cols = VCFSAMPLES + num_samples;
  StrBuf *vcfalts = &entry->cols[VCFALT], *vcfinfo = &entry->cols[VCFINFO];

  // Convert alts back into alt
  strbuf_reset(vcfalts);
  for(i = 0; i < entry->num_alts; i++) {
    if(entry->alts[i].len > 0) {
      if(vcfalts->len > 0) strbuf_append_char(vcfalts, ',');
      strbuf_append_strn(vcfalts, entry->alts[i].buff, entry->alts[i].len);
    }
  }

  // Convert info tags back into info
  strbuf_reset(vcfinfo);
  for(i = 0; i < entry->num_info; i++) {
    if(entry->info[i].len > 0) {
      if(vcfinfo->len > 0) strbuf_append_char(vcfinfo, ';');
      strbuf_append_strn(vcfinfo, entry->info[i].buff, entry->info[i].len);
    }
  }

  // Print
  fputs(entry->cols[0].buff, out);
  for(i = 1; i < num_cols; i++) {
    fputc('\t', out);
    fputs(entry->cols[i].buff, out);
  }
  fputc('\n', out);
}

void vcf_entry_add_filter(vcf_entry_t *entry, const char *status)
{
  StrBuf *filter = &entry->cols[VCFFILTER];
  if(filter->len == 0 || strcmp(filter->buff, ".") == 0 ||
     strcasecmp(filter->buff, "PASS") == 0)
  {
    strbuf_set(filter, status);
  }
  else
  {
    strbuf_append_char(filter, ';');
    strbuf_append_str(filter, status);
  }
}
