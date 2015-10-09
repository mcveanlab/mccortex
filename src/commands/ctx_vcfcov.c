#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
// #include "seq_reader.h"
#include "gpath_checks.h"
#include "seq_reader.h"
#include "genotyping.h"

#include "htslib/vcf.h"
#include "htslib/faidx.h"

// Don't attempt to genotype alleles bigger than this
size_t max_allele_len = 100;

// 2^8 = 256 possible haplotypes
uint32_t max_gt_vars = 8;

// Initial object buffer lengths
#define INIT_BUF_SIZE 128

// How many alt alleles to collect before printing
#define PRINT_BUF_LIMIT 100

const char vcfcov_usage[] =
"usage: "CMD" vcfcov [options] <in.vcf> <in.ctx> [in2.ctx ...]\n"
"\n"
"  Get coverage of a VCF in the cortex graphs. VCF must be sorted by position. \n"
"  It is recommended to use uncleaned graphs.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -o, --out <bub.txt.gz>  Output file [default: STDOUT]\n"
"  -r, --ref <ref.fa>      Reference file [required]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"ref",          required_argument, NULL, 'r'},
  {NULL, 0, NULL, 0}
};

char nkmers_ref_tag[10], nkmers_alt_tag[10];
char kcovgs_ref_tag[10], kcovgs_alt_tag[10];

typedef struct
{
  const char *path;
  const size_t *samplehdrids; // samplehdrids[x] is the sample id in VCF
  htsFile *vcffh;
  bcf_hdr_t *vcfhdr;
  // vpool are ready to be used
  GenoVCFPtrList vpool; // pool of vcf lines
  // index (vidx) of next GenoVCF to be printed
  // used to print VCF entries out in the correct (input) order
  size_t nxtprint, nextidx;
  // alist are current alleles; anchrom are on next chromosome
  // aprint are waiting to be printed; apool is a memory pool
  GenoVarPtrList alist, anchrom, aprint, apool; // alleles
  // Genotyping buffers
  int32_t *nkmers_r, *nkmers_a, *kcovgs_r, *kcovgs_a;
  size_t geno_buf_size;
  // Are annotations already in the graph?
  bool fetch_existing_tags;
  // Graph to get kmer coverage from
  const dBGraph *db_graph;
} VcfReader;


/**
 * Get coverage of ref and alt alleles from the de Bruijn graph in the given
 * colour.
 */
static inline void bkey_get_covg(BinaryKmer bkey, uint64_t altref_bits,
                                 GenoVar **gts, size_t ntgts,
                                 const dBGraph *db_graph)
{
  size_t i, col, ncols = db_graph->num_of_cols;
  dBNode node = db_graph_find(db_graph, bkey);
  uint64_t arbits;
  Covg covg;

  if(node.key != HASH_NOT_FOUND) {
    // printf("node: ");
    // db_nodes_print(&node, 1, db_graph, stdout);
    // printf("\n");

    for(col = 0; col < ncols; col++) {
      covg = db_node_get_covg(db_graph, node.key, col);
      if(!covg) continue;

      for(i = 0, arbits = altref_bits; i < ntgts; i++, arbits >>= 2) {
        if((arbits & 3) == 1) {
          gts[i]->c[col].nkmers[0]++;
          gts[i]->c[col].sumcovg[0] += covg;
        } else if((arbits & 3) == 2) {
          gts[i]->c[col].nkmers[1]++;
          gts[i]->c[col].sumcovg[1] += covg;
        } else { /* ignore kmer in both ref/alt */ }
      }
    }
  }
}

// return true if valid variant
static inline bool init_new_var(GenoVar *var, GenoVCF *parent, uint32_t aid)
{
  bcf1_t *v = &parent->v;
  uint32_t pos = v->pos, reflen, altlen;
  const char *ref = v->d.allele[0], *alt = v->d.allele[aid];
  reflen = strlen(ref);
  altlen = strlen(alt);

  // Left trim
  while(reflen && altlen && *ref == *alt) {
    ref++; alt++; pos++;
    reflen--; altlen--;
  }

  // Right trim
  while(reflen && altlen && ref[reflen-1] == alt[altlen-1]) {
    reflen--; altlen--;
  }

  if(reflen > max_allele_len || altlen > max_allele_len) return false;
  if(!reflen && !altlen) return false;

  // Initialise
  var->parent = parent;
  var->ref = ref;
  var->alt = alt;
  var->pos = pos;
  var->reflen = reflen;
  var->altlen = altlen;
  var->aid = aid;

  // Increment number of children
  parent->nchildren++;

  return true;
}

static void vcf_list_populate(GenoVCFPtrList *vlist, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) {
    GenoVCF *v = ctx_calloc(1, sizeof(GenoVCF));
    genovcf_ptr_list_append(vlist, v);
  }
}

static void var_list_populate(GenoVarPtrList *alist, size_t n, size_t ncols)
{
  size_t i;
  for(i = 0; i < n; i++) {
    GenoVar *a = ctx_calloc(1, sizeof(GenoVar));
    a->c = ctx_calloc(ncols, sizeof(VarCovg));
    genovar_ptr_list_append(alist, a);
  }
}

static void vcf_list_destroy(GenoVCFPtrList *vlist)
{
  size_t i;
  for(i = 0; i < genovcf_ptr_list_len(vlist); i++) {
    GenoVCF *v = genovcf_ptr_list_get(vlist, i);
    bcf_empty(&v->v);
    ctx_free(v);
  }
  genovcf_ptr_list_shift(vlist, NULL, genovcf_ptr_list_len(vlist));
}

static void var_list_destroy(GenoVarPtrList *alist)
{
  size_t i;
  for(i = 0; i < genovar_ptr_list_len(alist); i++) {
    GenoVar *a = genovar_ptr_list_get(alist, i);
    ctx_free(a->c);
    ctx_free(a);
  }
  genovar_ptr_list_shift(alist, NULL, genovar_ptr_list_len(alist));
}

static inline void vcfr_alloc(VcfReader *vcfr, const char *path,
                              htsFile *vcffh, bcf_hdr_t *vcfhdr,
                              const size_t *samplehdrids,
                              const dBGraph *db_graph)
{
  memset(vcfr, 0, sizeof(*vcfr));
  vcfr->path = path;
  vcfr->samplehdrids = samplehdrids;
  vcfr->vcffh = vcffh;
  vcfr->vcfhdr = vcfhdr;
  vcfr->db_graph = db_graph;
  vcfr->geno_buf_size = 0;

  genovcf_ptr_list_alloc(&vcfr->vpool, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->alist, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->anchrom, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->apool, INIT_BUF_SIZE);

  vcf_list_populate(&vcfr->vpool, INIT_BUF_SIZE);
  var_list_populate(&vcfr->apool, INIT_BUF_SIZE, db_graph->num_of_cols);
}

static inline void vcfr_dealloc(VcfReader *vcfr)
{
  ctx_assert(genovar_ptr_list_len(&vcfr->alist) == 0);
  ctx_assert(genovar_ptr_list_len(&vcfr->anchrom) == 0);
  vcf_list_destroy(&vcfr->vpool);
  var_list_destroy(&vcfr->apool);
  ctx_assert(genovcf_ptr_list_len(&vcfr->vpool) == 0);
  ctx_assert(genovar_ptr_list_len(&vcfr->apool) == 0);

  genovcf_ptr_list_dealloc(&vcfr->vpool);
  genovar_ptr_list_dealloc(&vcfr->alist);
  genovar_ptr_list_dealloc(&vcfr->anchrom);
  genovar_ptr_list_dealloc(&vcfr->apool);

  ctx_free(vcfr->nkmers_r);
  ctx_free(vcfr->nkmers_a);
  ctx_free(vcfr->kcovgs_r);
  ctx_free(vcfr->kcovgs_a);
}

static inline void fetch_chrom(bcf_hdr_t *hdr, bcf1_t *v,
                               faidx_t *fai, int *refid,
                               char **chrom, int *chromlen)
{
  if(*refid != v->rid) {
    free(*chrom);
    *chrom = fai_fetch(fai, bcf_seqname(hdr, v), chromlen);
    if(*chrom == NULL) die("Cannot find chrom '%s'", bcf_seqname(hdr, v));
    *refid = v->rid;
  }
}

// Returns:
//  -1 if EOF / error
//  0 if read would have added var on another chrom
//  Otherwise number of variants added
//
static int vcfr_fetch(VcfReader *vr)
{
  // Check we have GenoVCF entries to read into
  if(genovcf_ptr_list_len(&vr->vpool) == 0) vcf_list_populate(&vr->vpool, 16);

  // If we already have alleles from a diff chrom to use, return them
  if(genovar_ptr_list_len(&vr->alist) == 0 && genovar_ptr_list_len(&vr->anchrom))
  {
    genovar_ptr_list_push(&vr->alist, genovar_ptr_list_getptr(&vr->anchrom, 0),
                                      genovar_ptr_list_len(&vr->anchrom));
    genovar_ptr_list_reset(&vr->anchrom);
    return genovar_ptr_list_len(&vr->alist);
  }

  // Take vcf out of pool
  GenoVCF *ve;
  genovcf_ptr_list_pop(&vr->vpool, &ve, 1);
  ve->vidx = vr->nextidx++;
  ve->nchildren = 0;
  bcf1_t *v = &ve->v;

  while(1)
  {
    // Read VCF
    if(bcf_read(vr->vcffh, vr->vcfhdr, v) < 0) {
      // EOF
      genovcf_ptr_list_push(&vr->vpool, &ve, 1);
      return -1;
    }

    // Unpack all info
    bcf_unpack(v, BCF_UN_ALL);

    // Check we have enough vars to decompose
    size_t i, n = MAX2(v->n_allele, 16);
    if(genovar_ptr_list_len(&vr->apool) < n)
      var_list_populate(&vr->apool, n, vr->db_graph->num_of_cols);

    size_t nadded = 0, nexisting_alleles = genovar_ptr_list_len(&vr->alist);
    bool diff_chroms = false, overlap = false;
    GenoVar *lvar = NULL;

    if(nexisting_alleles) {
      lvar = genovar_ptr_list_get(&vr->alist, nexisting_alleles-1);
      diff_chroms = (lvar->parent->v.rid != v->rid);
      int32_t var_end = lvar->parent->v.pos + strlen(lvar->parent->v.d.allele[0]);
      overlap = (!diff_chroms && var_end > v->pos);

      // Check VCF is sorted
      if(!diff_chroms && lvar->parent->v.pos > v->pos) {
        die("VCF is not sorted: %s:%lli", vr->path, vr->vcffh->lineno);
      }
    }

    ctx_assert2(!diff_chroms || genovar_ptr_list_len(&vr->anchrom) == 0,
                "Already read diff chrom");

    // Load alleles into alist
    GenoVar *var;
    genovar_ptr_list_pop(&vr->apool, &var, 1);
    GenoVarPtrList *list = diff_chroms ? &vr->anchrom : &vr->alist;

    // i==0 is ref allele
    for(i = 1; i < v->n_allele; i++) {
      if(init_new_var(var, ve, i)) {
        genovar_ptr_list_push(list, &var, 1);
        genovar_ptr_list_pop(&vr->apool, &var, 1);
        nadded++;
      }
    }
    // Re-add unused var back into pool
    genovar_ptr_list_push(&vr->apool, &var, 1);

    if(nadded)
    {
      if(overlap || diff_chroms) {
        genovars_sort(genovar_ptr_list_getptr(list, 0),
                      genovar_ptr_list_len(list));
      } else {
        // Just sort the alleles we added to the end of alist
        genovars_sort(genovar_ptr_list_getptr(&vr->alist, nexisting_alleles),
                      genovar_ptr_list_len(&vr->alist) - nexisting_alleles);
      }

      return diff_chroms ? 0 : nadded;
    }
  }
}

static void vcfr_drop_var(VcfReader *vr, size_t idx)
{
  GenoVar *a = genovar_ptr_list_get(&vr->alist, idx);
  genovar_ptr_list_append(&vr->aprint, a);
  // Instead of removing `a` from sorted array vr->alist, set it to NULL
  genovar_ptr_list_set(&vr->alist, idx, NULL);
}

// Remove NULL entries from vr->alist
static void vcfr_shrink_vars(VcfReader *vr)
{
  size_t i, j, len = genovar_ptr_list_len(&vr->alist);
  GenoVar *var;
  for(i = j = 0; i < len; i++) {
    var = genovar_ptr_list_get(&vr->alist, i);
    if(var != NULL) {
      genovar_ptr_list_set(&vr->alist, j, var);
      j++;
    }
  }
  genovar_ptr_list_pop(&vr->alist, NULL, i-j);
}

static int _genovar_cmp_vidx(const void *aa, const void *bb)
{
  const GenoVar *a = *(const GenoVar*const*)aa, *b = *(const GenoVar*const*)bb;
  return (long)a->parent->vidx - (long)b->parent->vidx;
}

static void vcfr_print_entry(VcfReader *vfr, htsFile *outfh, bcf_hdr_t *outhdr,
                             GenoVar **alleles, size_t nalleles)
{
  GenoVCF *var = alleles[0]->parent;
  bcf1_t *v = &var->v;
  size_t nsamples = bcf_hdr_nsamples(outhdr);
  size_t i, col, aid, sid, nalts = v->n_allele-1, n = nsamples * nalts;
  size_t ncols = vfr->db_graph->num_of_cols;

  // Check if first time printing
  bool first_print = (vfr->geno_buf_size == 0);

  // _r ref, _a alt
  if(vfr->geno_buf_size < n)
  {
    vfr->kcovgs_r = ctx_reallocarray(vfr->kcovgs_r, n, sizeof(int32_t));
    vfr->kcovgs_a = ctx_reallocarray(vfr->kcovgs_a, n, sizeof(int32_t));
    vfr->nkmers_r = ctx_reallocarray(vfr->nkmers_r, n, sizeof(int32_t));
    vfr->nkmers_a = ctx_reallocarray(vfr->nkmers_a, n, sizeof(int32_t));
    // Initiate to missing
    for(i = vfr->geno_buf_size; i < n; i++) {
      vfr->kcovgs_r[i] = vfr->kcovgs_a[i] = bcf_int32_missing;
      vfr->nkmers_r[i] = vfr->nkmers_a[i] = bcf_int32_missing;
    }
    vfr->geno_buf_size = n;
  }

  // Fetch existing coverage from VCF
  if(vfr->fetch_existing_tags || first_print)
  {
    int nsize = vfr->geno_buf_size;
    bcf_get_format_int32(vfr->vcfhdr, v, nkmers_ref_tag, &vfr->nkmers_r, &nsize);
    bcf_get_format_int32(vfr->vcfhdr, v, nkmers_alt_tag, &vfr->nkmers_a, &nsize);
    bcf_get_format_int32(vfr->vcfhdr, v, kcovgs_ref_tag, &vfr->kcovgs_r, &nsize);
    bcf_get_format_int32(vfr->vcfhdr, v, kcovgs_alt_tag, &vfr->kcovgs_a, &nsize);
    if(first_print && nsize > 0) vfr->fetch_existing_tags = true;
    ctx_assert2(nsize == (int)vfr->geno_buf_size, "htslib resized our buffer!");
  }

  // Add coverage
  VarCovg *cov;

  for(i = 0; i < nalleles; i++) {
    aid = alleles[i]->aid;
    ctx_assert(aid < v->n_allele);

    if(!alleles[i]->has_covg) {
      for(col = 0; col < ncols; col++) {
        sid = vfr->samplehdrids[col];
        vfr->nkmers_r[sid*nalts+aid-1] = bcf_int32_missing;
        vfr->nkmers_a[sid*nalts+aid-1] = bcf_int32_missing;
        vfr->kcovgs_r[sid*nalts+aid-1] = bcf_int32_missing;
        vfr->kcovgs_a[sid*nalts+aid-1] = bcf_int32_missing;
      }
    }
    else {
      for(col = 0; col < ncols; col++) {
        sid = vfr->samplehdrids[col];
        cov = &alleles[i]->c[col];
        vfr->nkmers_r[sid*nalts+aid-1] = cov->nkmers[0];
        vfr->nkmers_a[sid*nalts+aid-1] = cov->nkmers[1];
        vfr->kcovgs_r[sid*nalts+aid-1] = cov->nkmers[0] ? cov->sumcovg[0] / cov->nkmers[0] : 0;
        vfr->kcovgs_a[sid*nalts+aid-1] = cov->nkmers[1] ? cov->sumcovg[1] / cov->nkmers[1] : 0;
      }
    }
  }

  // Update VCF entry
  bcf_update_format_int32(outhdr, v, nkmers_ref_tag, vfr->nkmers_r, n);
  bcf_update_format_int32(outhdr, v, nkmers_alt_tag, vfr->nkmers_a, n);
  bcf_update_format_int32(outhdr, v, kcovgs_ref_tag, vfr->kcovgs_r, n);
  bcf_update_format_int32(outhdr, v, kcovgs_alt_tag, vfr->kcovgs_a, n);

  if(bcf_write(outfh, outhdr, v) != 0) die("Cannot write record");
}

static void vcfr_print_waiting(VcfReader *vr, htsFile *outfh, bcf_hdr_t *outhdr,
                               bool force)
{
  // Sort waiting by vidx
  size_t num_a, alen = genovar_ptr_list_len(&vr->aprint);
  if(alen == 0 || (!force && alen < PRINT_BUF_LIMIT)) return;

  GenoVar **alleleptr = genovar_ptr_list_getptr(&vr->aprint, 0);
  qsort(alleleptr, alen, sizeof(GenoVar*), _genovar_cmp_vidx);

  while(alen > 0)
  {
    alleleptr = genovar_ptr_list_getptr(&vr->aprint, 0);
    alen = genovar_ptr_list_len(&vr->aprint);

    num_a = 0;
    while(num_a < alen && alleleptr[num_a]->parent->vidx == vr->nxtprint)
      num_a++;

    if(num_a < alleleptr[0]->parent->nchildren) break;
    ctx_assert(num_a == alleleptr[0]->parent->nchildren && num_a > 0);

    // print
    vcfr_print_entry(vr, outfh, outhdr, alleleptr, num_a);

    // Re-add to vpool
    genovcf_ptr_list_append(&vr->vpool, alleleptr[0]->parent);

    // Re-add to apool
    genovar_ptr_list_unshift(&vr->apool, alleleptr, num_a);
    genovar_ptr_list_shift(&vr->aprint, NULL, num_a);

    vr->nxtprint++;
    alen -= num_a;
  }
}

static void vcfcov_vars(GenoVar **vars, size_t nvars,
                          size_t tgtidx, size_t ntgts,
                          const char *chrom, size_t chromlen,
                          Genotyper *gtyper, const dBGraph *db_graph)
{
  ctx_assert(ntgts <= nvars);

  GenoKmer *kmers = NULL;
  size_t i, end, nkmers;

  // Zero coverage on alleles before querying graph
  for(i = tgtidx, end = tgtidx+ntgts; i < end; i++)
    genovar_wipe_covg(vars[i], db_graph->num_of_cols);

  if(nvars > max_gt_vars) { return; }

  nkmers = genotyping_get_kmers(gtyper, (const GenoVar *const*)vars,
                                nvars, tgtidx, ntgts,
                                chrom, chromlen,
                                db_graph->kmer_size, &kmers);

  for(i = 0; i < nkmers; i++) {
    bkey_get_covg(kmers[i].bkey, kmers[i].arbits,
                  vars+tgtidx, ntgts,
                  db_graph);
  }

  for(i = tgtidx, end = tgtidx+ntgts; i < end; i++)
    vars[i]->has_covg = true;
}

// Get index of first of vars (starting at index i), which starts after endpos
static inline int get_vcf_cov_end(GenoVar **vars, size_t nvars,
                                  size_t i, size_t endpos)
{
  while(i < nvars && endpos > vars[i]->pos) { i++; }
  return i;
}

// vars will be unordered after return
static void vcfcov_block(GenoVar **vars, size_t nvars,
                         const char *chrom, int chromlen,
                         Genotyper *gtyper, const dBGraph *db_graph)
{
  if(nvars <= max_gt_vars)
  {
    vcfcov_vars(vars, nvars, 0, nvars, chrom, chromlen, gtyper, db_graph);
  }
  else
  {
    // do one at a time
    const int ks = db_graph->kmer_size;
    int i, j;
    // genotype start/end, background start/end (end is not inclusive)
    size_t gs = 0, ge, bs, be, tmp_ge, tmp_be;

    while(gs < nvars)
    {
      // Get vars to the left of our genotyping-start (gs)
      for(i = j = (int)gs-1; i >= 0; i--) {
        if(genovar_end(vars[i]) + ks > vars[gs]->pos) {
          SWAP(vars[i], vars[j]);
          j--;
        }
      }

      bs = j+1;
      ge = gs+1;
      be = get_vcf_cov_end(vars, nvars, ge, genovar_end(vars[ge-1])+ks);

      // Try increasing ge, calculate be, accept if small enough
      for(tmp_ge = ge+1; tmp_ge < nvars; tmp_ge++) {
        tmp_be = get_vcf_cov_end(vars, nvars, tmp_ge, genovar_end(vars[tmp_ge-1])+ks);
        if(tmp_be - bs <= max_gt_vars) { ge = tmp_ge; be = tmp_be; }
        else break;
      }

      ctx_assert2(bs<=gs && gs<=ge && ge<=be, "%zu %zu %zu %zu",bs,gs,ge,be);

      // status("bs:%zu gs:%zu ge:%zu be:%zu", bs, gs, ge, be);
      vcfcov_vars(vars+bs, be-bs, gs-bs, ge-gs, chrom, chromlen, gtyper, db_graph);

      gs = ge;
    }
  }
}

// samplehdrids[col] is the index of a colour in the output vcf
static void vcfcov_file(htsFile *vcffh, bcf_hdr_t *vcfhdr,
                        htsFile *outfh, bcf_hdr_t *outhdr,
                        const char *path, faidx_t *fai,
                        const size_t *samplehdrids,
                        const dBGraph *db_graph)
{
  VcfReader vr;
  vcfr_alloc(&vr, path, vcffh, vcfhdr, samplehdrids, db_graph);

  Genotyper *gtyper = genotyper_init();

  // refid is id of chromosome currently loaded
  char *chrom = NULL;
  int refid = -1, chromlen = 0;

  int n;
  const size_t kmer_size = db_graph->kmer_size;
  size_t i, start, end, endpos;

  GenoVar **alist;
  size_t alen;

  while((n = vcfr_fetch(&vr)) >= 0)
  {
    // Get ref chromosome
    // Only loads if we don't currently have the right chrom
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chrom, &chromlen);

    alist = genovar_ptr_list_getptr(&vr.alist, 0);
    alen = genovar_ptr_list_len(&vr.alist);
    ctx_assert(alen > 0);

    // Break into blocks
    endpos = genovar_end(alist[0]) + kmer_size;
    for(start = 0, end = 1; end < alen; end++) {
      if(endpos <= alist[end]->pos) {
        // end of block
        vcfcov_block(alist+start, end-start, chrom, chromlen, gtyper, db_graph);
        for(i = start; i < end; i++) vcfr_drop_var(&vr, i);
        start = end;
      }
      endpos = MAX2(endpos, genovar_end(alist[end])+kmer_size);
    }

    if(n == 0) {
      // End of chromosome -- do all variants
      vcfcov_block(alist+start, alen-start, chrom, chromlen, gtyper, db_graph);
      for(i = start; i < alen; i++) vcfr_drop_var(&vr, i);
    }

    // Shrink array if we dropped any vars
    vcfr_shrink_vars(&vr);
    vcfr_print_waiting(&vr, outfh, outhdr, false);
  }

  // Deal with remainder
  alist = genovar_ptr_list_getptr(&vr.alist, 0);
  alen = genovar_ptr_list_len(&vr.alist);
  if(alen) {
    // Get ref chromosome
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chrom, &chromlen);
    vcfcov_block(alist, alen, chrom, chromlen, gtyper, db_graph);
    for(i = 0; i < alen; i++) vcfr_drop_var(&vr, i);
    vcfr_shrink_vars(&vr);
  }
  vcfr_print_waiting(&vr, outfh, outhdr, true);

  free(chrom);
  genotyper_destroy(gtyper);
  vcfr_dealloc(&vr);
}

int ctx_vcfcov(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;

  char *ref_path = NULL;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;
  size_t i;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'r': cmd_check(!ref_path, cmd); ref_path = optarg; break;
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" geno -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(ref_path == NULL) cmd_print_usage("Require a reference (-r,--ref <ref.fa>)");
  if(optind+2 > argc) cmd_print_usage("Require VCF and graph files");

  // open ref
  // index fasta with: samtools faidx ref.fa
  faidx_t *fai = fai_load(ref_path);
  if(fai == NULL) die("Cannot load ref index: %s / %s.fai", ref_path, ref_path);

  // Open input VCF file
  const char *vcf_path = argv[optind++];
  htsFile *vcffh = hts_open(vcf_path, "r");
  bcf_hdr_t *vcfhdr = bcf_hdr_read(vcffh);

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;
  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, NULL, 0, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Covg) * ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        -1, -1,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  futil_create_output(out_path);
  htsFile *outfh = hts_open(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_COVGS);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  size_t *samplehdrids = ctx_malloc(db_graph.num_of_cols * sizeof(size_t));

  // Add samples to vcf header
  bcf_hdr_t *outhdr = bcf_hdr_dup(vcfhdr);
  int sid;
  for(i = 0; i < db_graph.num_of_cols; i++) {
    char *sname = db_graph.ginfo[i].sample_name.b;
    if((sid = bcf_hdr_id2int(outhdr, BCF_DT_SAMPLE, sname)) < 0) {
      bcf_hdr_add_sample(outhdr, sname);
      sid = bcf_hdr_id2int(outhdr, BCF_DT_SAMPLE, sname);
    }
    samplehdrids[i] = sid;
    status("Colour %zu: %s [VCF column %zu]", i, sname, samplehdrids[i]);
  }

  sprintf(nkmers_ref_tag, "NK%zuR", db_graph.kmer_size);
  sprintf(nkmers_alt_tag, "NK%zuA", db_graph.kmer_size);
  sprintf(kcovgs_ref_tag, "CK%zuR", db_graph.kmer_size);
  sprintf(kcovgs_alt_tag, "CK%zuA", db_graph.kmer_size);

  // Add genotype format fields
  // One field per alternative allele
  char descr[200];
  sprintf(descr,
          "##FORMAT=<ID=%s,Number=A,Type=Integer,"
          "Description=\"Number of exclusive kmers on ref for each allele (k=%zu)\">\n",
          nkmers_ref_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);
  sprintf(descr,
          "##FORMAT=<ID=%s,Number=A,Type=Integer,"
          "Description=\"Number of exclusive kmers on alt for each allele (k=%zu)\">\n",
          nkmers_alt_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);
  sprintf(descr,
          "##FORMAT=<ID=%s,Number=A,Type=Integer,Description=\"Mean ref exclusive kmer coverage (k=%zu)\">\n",
          kcovgs_ref_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);
  sprintf(descr,
          "##FORMAT=<ID=%s,Number=A,Type=Integer,Description=\"Mean alt exclusive kmer coverage (k=%zu)\">\n",
          kcovgs_alt_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);

  bcf_hdr_set_version(outhdr, "VCFv4.2");

  if(bcf_hdr_write(outfh, outhdr) != 0)
    die("Cannot write header to: %s", futil_outpath_str(out_path));

  status("[vcfcov] Reading %s and adding coverage", vcf_path);
  vcfcov_file(vcffh, vcfhdr, outfh, outhdr, vcf_path, fai,
              samplehdrids, &db_graph);

  status("  saved to: %s\n", out_path);

  ctx_free(samplehdrids);

  bcf_hdr_destroy(vcfhdr);
  bcf_hdr_destroy(outhdr);
  hts_close(vcffh);
  hts_close(outfh);
  fai_destroy(fai);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
