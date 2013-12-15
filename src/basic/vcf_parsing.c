#include "global.h"
#include "vcf_parsing.h"
#include "util.h"

static void strbuf_arr_resize(StrBuf **arr, size_t *cap, size_t newcap)
{
  size_t i;
  newcap = ROUNDUP2POW(newcap);
  *arr = realloc2(*arr, sizeof(StrBuf) * newcap);
  for(i = *cap; i < newcap; i++) strbuf_alloc(&((*arr)[i]), 64);
  *cap = newcap;
}

// DP=asdf;TXT="a;b;c";F1;X=4
// ends:  ^           ^  ^   ^
char* vcf_info_tag_end(char *str)
{
  boolean speechmrks = false;
  for(; *str; str++) {
    if(*str == '"') speechmrks = !speechmrks;
    else if(*str == ';' && !speechmrks) return str;
  }
  return NULL;
}

StrBuf* vcf_info_tag_find(vcf_entry_t *entry, const char* str)
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

static inline void vcf_find_lf_rf_info(vcf_entry_t *entry)
{
  size_t i;
  entry->lf = entry->rf = NULL;
  for(i = 0; i < entry->num_info; i++) {
    if(!strncmp(entry->info[i].buff, "LF=", 3)) entry->lf = &entry->info[i];
    else if(!strncmp(entry->info[i].buff, "RF=", 3)) entry->rf = &entry->info[i];
  }
}

void vcf_info_tag_add(vcf_entry_t *entry, const char* fmt, ...)
{
  // Find emtpy tag
  size_t i; StrBuf *sbuf;
  for(i = 0; i < entry->num_info; i++) {
    sbuf = &entry->info[i];
    if(sbuf->buff[0] == '\0') break;
  }

  // Or add new tag
  if(i == entry->num_info) {
    vcf_entry_info_capacity(entry, entry->num_info+1);
    sbuf = &entry->info[entry->num_info++];
  }

  va_list argptr;
  va_start(argptr, fmt);
  strbuf_reset(sbuf);
  strbuf_vsprintf(sbuf, 0, fmt, argptr);
  va_end(argptr);
}

void vcf_info_tag_del(vcf_entry_t *entry, const char* tag)
{
  StrBuf *info = vcf_info_tag_find(entry, tag);
  if(info != NULL) strbuf_reset(info);
}

void vcf_entry_alt_capacity(vcf_entry_t *entry, size_t num_alts)
{
  if(num_alts > entry->alts_capacity)
    strbuf_arr_resize(&entry->alts, &entry->alts_capacity, num_alts);
}

void vcf_entry_info_capacity(vcf_entry_t *entry, size_t num_info)
{
  if(num_info > entry->info_capacity) {
    strbuf_arr_resize(&entry->info, &entry->info_capacity, num_info);
    // Re-find lf/rf
    vcf_find_lf_rf_info(entry);
  }
}

void vcf_entry_alloc(vcf_entry_t *entry, size_t num_samples)
{
  size_t i, num_cols = VCFSAMPLES+num_samples;
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
}

void vcf_entry_dealloc(vcf_entry_t *entry, size_t num_samples)
{
  size_t i, num_cols = VCFSAMPLES+num_samples;
  for(i = 0; i < num_cols; i++) strbuf_dealloc(&entry->cols[i]);
  for(i = 0; i < entry->alts_capacity; i++) strbuf_dealloc(&entry->alts[i]);
  for(i = 0; i < entry->info_capacity; i++) strbuf_dealloc(&entry->info[i]);
  
  free(entry->cols);
  free(entry->alts);
  free(entry->info);
}

void vcf_entry_cpy(vcf_entry_t *dst, const vcf_entry_t *src, size_t num_samples)
{
  vcf_entry_alt_capacity(dst, src->num_alts);
  vcf_entry_info_capacity(dst, src->num_info);

  size_t i;
  for(i = 0; i < VCFSAMPLES+num_samples; i++)
    strbuf_set(&dst->cols[i], src->cols[i].buff);

  // Copy alts
  dst->num_alts = src->num_alts;
  for(i = 0; i < src->num_alts; i++)
    strbuf_set(&dst->alts[i], src->alts[i].buff);

  // Copy info and set LF, RF
  dst->num_info = src->num_info;
  for(i = 0; i < src->num_info; i++)
    strbuf_set(&dst->info[i], src->info[i].buff);

  vcf_find_lf_rf_info(dst);
}

void vcf_entry_parse(StrBuf *line, vcf_entry_t *entry, size_t num_samples)
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

  vcf_entry_alt_capacity(entry, alts_count);

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

  size_t info_count = 0;
  while((tmp = vcf_info_tag_end(tmp)) != NULL) { info_count++; tmp++; }
  tmp = entry->cols[VCFINFO].buff;

  vcf_entry_info_capacity(entry, info_count);

  while(1)
  {
    end = vcf_info_tag_end(tmp);
    if(end != NULL) *end = '\0';
    buf = &entry->info[entry->num_info++];
    strbuf_set(buf, tmp);
    if(end != NULL) *end = ';';
    else break;
    tmp = end + 1;
  }

  vcf_find_lf_rf_info(entry);
}

void vcf_entry_revcmp(vcf_entry_t *entry)
{
  // For debugging
  vcf_info_tag_add(entry, "BUBREV");

  size_t i;
  for(i = 0; i < entry->num_alts; i++)
    reverse_complement_str(entry->alts[i].buff, entry->alts[i].len);

  // reverse complement and swap lflank and rflank
  // printf("lf:'%s'; rf:'%s'\n", entry->lf.buff, entry->rf.buff);
  reverse_complement_str(entry->lf->buff+3, entry->lf->len-3);
  reverse_complement_str(entry->rf->buff+3, entry->rf->len-3);
  StrBuf *tmpbuf;
  SWAP(entry->lf, entry->rf, tmpbuf);
  entry->lf->buff[0] = 'L';
  entry->rf->buff[0] = 'R';
}

size_t vcf_entry_longest_allele(const vcf_entry_t *entry)
{
  size_t i, max = 0;
  for(i = 0; i < entry->num_alts; i++) {
    max = MAX2(max, entry->alts[i].len);
  }
  return max;
}

void vcf_entry_print(const vcf_entry_t *entry, FILE *out, size_t num_samples)
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
