#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_checks.h"
#include "genotyping.h"

#include "htslib/vcf.h"
#include "htslib/faidx.h"

#define DEFAULT_MAX_ALLELE_LEN 100
#define DEFAULT_MAX_GT_VARS 8

// Don't attempt to genotype alleles bigger than this
// defaults to DEFAULT_MAX_ALLELE_LEN
size_t max_allele_len = 0;

// 2^8 = 256 possible haplotypes
// defaults to DEFAULT_MAX_GT_VARS
uint32_t max_gt_vars = 0;

// Initial object buffer lengths
#define INIT_BUF_SIZE 128

// How many alt alleles to collect before printing
#define PRINT_BUF_LIMIT 100

// Stats
size_t num_vcf_lines = 0, num_alts_read = 0, num_alts_loaded = 0;
size_t num_alts_too_long = 0, num_alts_no_covg = 0, num_alts_with_covg = 0;

const char vcfcov_usage[] =
"usage: "CMD" vcfcov [options] <in.vcf> <in.ctx> [in2.ctx ...]\n"
"\n"
"  Get coverage of a VCF in the cortex graphs. VCF must be sorted by position. \n"
"  It is recommended to use uncleaned graphs.\n"
"\n"
"  -h, --help                This help message\n"
"  -q, --quiet               Silence status output normally printed to STDERR\n"
"  -f, --force               Overwrite output files\n"
"  -m, --memory <mem>        Memory to use\n"
"  -n, --nkmers <kmers>      Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -o, --out <bub.txt.gz>    Output file [default: STDOUT]\n"
"  -O, --out-fmt <f>         Format vcf|vcfgz|bcf|ubcf\n"
"  -r, --ref <ref.fa>        Reference file [required]\n"
"  -L, --max-var-len <A>     Only use alleles <= A bases long [default: "QUOTE_VALUE(DEFAULT_MAX_ALLELE_LEN)"]\n"
"  -N, --max-nvars <N>       Limit haplotypes to <= N variants [default: "QUOTE_VALUE(DEFAULT_MAX_GT_VARS)"]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"out-fmt",      required_argument, NULL, 'O'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"ref",          required_argument, NULL, 'r'},
  {"max-var-len",  required_argument, NULL, 'L'},
  {"max-nvars",    required_argument, NULL, 'N'},
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
  VcfCovLine *curr_line; // last line to be read into alist
  // vpool are ready to be used
  VcfCovLinePtrList vpool; // pool of vcf lines
  // index (vidx) of next VcfCovLine to be printed
  // used to print VCF entries out in the correct (input) order
  size_t nextidx, nxtprint;
  // alist are current alleles; anchrom are on next chromosome
  // aprint are waiting to be printed; apool is a memory pool
  VcfCovAltPtrList alist, anchrom, aprint, apool; // alleles
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
                                 VcfCovAlt **gts, size_t ntgts,
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
static inline bool init_new_alt(VcfCovAlt *var, VcfCovLine *line, uint32_t aid,
                                size_t ncols)
{
  bcf1_t *v = &line->v;
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

  if(reflen > max_allele_len || altlen > max_allele_len) {
    num_alts_too_long++;
    return false;
  }

  if(!reflen && !altlen) return false;

  // Initialise
  var->parent = line;
  var->ref = ref;
  var->alt = alt;
  var->pos = pos;
  var->reflen = reflen;
  var->altlen = altlen;
  var->aid = aid;
  var->has_covg = false;
  memset(var->c, 0, sizeof(var->c[0]) * ncols);

  // Increment number of children
  line->nchildren++;

  return true;
}

static void vcf_list_populate(VcfCovLinePtrList *vlist, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) {
    VcfCovLine *v = ctx_calloc(1, sizeof(VcfCovLine));
    vc_lines_append(vlist, v);
  }
}

static void var_list_populate(VcfCovAltPtrList *alist, size_t n, size_t ncols)
{
  size_t i;
  for(i = 0; i < n; i++) {
    VcfCovAlt *a = ctx_calloc(1, sizeof(VcfCovAlt));
    a->c = ctx_calloc(ncols, sizeof(VarCovg));
    vc_alts_append(alist, a);
  }
}

static void vcf_list_destroy(VcfCovLinePtrList *vlist)
{
  size_t i;
  for(i = 0; i < vc_lines_len(vlist); i++) {
    VcfCovLine *v = vc_lines_get(vlist, i);
    bcf_empty(&v->v);
    ctx_free(v);
  }
  vc_lines_shift(vlist, NULL, vc_lines_len(vlist));
}

static void var_list_destroy(VcfCovAltPtrList *alist)
{
  size_t i;
  for(i = 0; i < vc_alts_len(alist); i++) {
    VcfCovAlt *a = vc_alts_get(alist, i);
    ctx_free(a->c);
    ctx_free(a);
  }
  vc_alts_shift(alist, NULL, vc_alts_len(alist));
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

  vc_lines_alloc(&vcfr->vpool, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->alist, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->anchrom, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->apool, INIT_BUF_SIZE);

  vcf_list_populate(&vcfr->vpool, INIT_BUF_SIZE);
  var_list_populate(&vcfr->apool, INIT_BUF_SIZE, db_graph->num_of_cols);
}

static inline void vcfr_dealloc(VcfReader *vcfr)
{
  ctx_assert(vc_alts_len(&vcfr->alist) == 0);
  ctx_assert(vc_alts_len(&vcfr->anchrom) == 0);
  vcf_list_destroy(&vcfr->vpool);
  var_list_destroy(&vcfr->apool);
  ctx_assert(vc_lines_len(&vcfr->vpool) == 0);
  ctx_assert(vc_alts_len(&vcfr->apool) == 0);

  vc_lines_dealloc(&vcfr->vpool);
  vc_alts_dealloc(&vcfr->alist);
  vc_alts_dealloc(&vcfr->anchrom);
  vc_alts_dealloc(&vcfr->apool);

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
  // Check we have VcfCovLine entries to read into
  if(vc_lines_len(&vr->vpool) == 0) vcf_list_populate(&vr->vpool, 16);

  // If we already have alleles from a diff chrom to use, return them
  if(vc_alts_len(&vr->alist) == 0 && vc_alts_len(&vr->anchrom))
  {
    vc_alts_push(&vr->alist, vc_alts_getptr(&vr->anchrom, 0),
                             vc_alts_len(&vr->anchrom));
    vc_alts_reset(&vr->anchrom);
    VcfCovAlt **alts = vc_alts_getptr(&vr->alist, 0);
    vr->curr_line = alts[0]->parent;
    return vc_alts_len(&vr->alist);
  }

  // Take vcf out of pool
  VcfCovLine *line;
  vc_lines_pop(&vr->vpool, &line, 1);
  line->vidx = vr->nextidx++;
  line->nchildren = 0;
  bcf1_t *v = &line->v;

  while(1)
  {
    // Read VCF
    if(bcf_read(vr->vcffh, vr->vcfhdr, v) < 0) {
      // EOF
      vc_lines_push(&vr->vpool, &line, 1);
      return -1;
    }

    // Unpack all info
    bcf_unpack(v, BCF_UN_ALL);
    num_vcf_lines++;
    vr->curr_line = line;

    // Check we have enough vars to decompose
    size_t i, n = MAX2(v->n_allele, 16);
    if(vc_alts_len(&vr->apool) < n)
      var_list_populate(&vr->apool, n, vr->db_graph->num_of_cols);

    size_t nadded = 0, nprev_alts = vc_alts_len(&vr->alist);
    bool diff_chroms = false, overlap = false;
    VcfCovAlt *lvar = NULL;

    if(nprev_alts) {
      lvar = vc_alts_get(&vr->alist, nprev_alts-1);
      diff_chroms = (lvar->parent->v.rid != v->rid);
      int32_t var_end = lvar->parent->v.pos + strlen(lvar->parent->v.d.allele[0]);
      overlap = (!diff_chroms && var_end > v->pos);

      // Check VCF is sorted
      if(!diff_chroms && lvar->parent->v.pos > v->pos) {
        die("VCF is not sorted: %s:%lli", vr->path, (long long)vr->vcffh->lineno);
      }
    }

    ctx_assert2(!diff_chroms || vc_alts_len(&vr->anchrom) == 0,
                "Already read diff chrom");

    // Load alleles into alist
    VcfCovAlt *alt;
    vc_alts_pop(&vr->apool, &alt, 1);
    VcfCovAltPtrList *list = diff_chroms ? &vr->anchrom : &vr->alist;

    // i==0 is ref allele
    for(i = 1; i < v->n_allele; i++) {
      if(init_new_alt(alt, line, i, vr->db_graph->num_of_cols)) {
        vc_alts_push(list, &alt, 1);
        vc_alts_pop(&vr->apool, &alt, 1);
        nadded++;
        num_alts_loaded++;
      }
    }

    // Re-add unused var back into pool
    vc_alts_push(&vr->apool, &alt, 1);

    num_alts_read += v->n_allele - 1; // v_allele includes ref

    if(nadded)
    {
      if(overlap || diff_chroms) {
        vcfcov_alts_sort(vc_alts_getptr(list, 0),
                         vc_alts_len(list));
      } else {
        // Just sort the alleles we added to the end of alist
        vcfcov_alts_sort(vc_alts_getptr(&vr->alist, nprev_alts),
                         vc_alts_len(&vr->alist) - nprev_alts);
      }

      return diff_chroms ? 0 : nadded;
    }
  }
}

// Move a variant from alist -> print list
static void vcfr_move_to_print(VcfReader *vr, size_t idx)
{
  VcfCovAlt *a = vc_alts_get(&vr->alist, idx);
  vc_alts_append(&vr->aprint, a);
  // Instead of removing `a` from sorted array vr->alist, set it to NULL
  vc_alts_set(&vr->alist, idx, NULL);
}

// Remove NULL entries from vr->alist
static void vcfr_repack_alts(VcfReader *vr)
{
  size_t i, j, len = vc_alts_len(&vr->alist);
  VcfCovAlt *var;
  for(i = j = 0; i < len; i++) {
    var = vc_alts_get(&vr->alist, i);
    if(var != NULL) {
      vc_alts_set(&vr->alist, j, var);
      j++;
    }
  }
  vc_alts_pop(&vr->alist, NULL, i-j);
}

static int _vcfcov_alt_cmp_vidx(const void *aa, const void *bb)
{
  const VcfCovAlt *a = *(const VcfCovAlt*const*)aa;
  const VcfCovAlt *b = *(const VcfCovAlt*const*)bb;
  return (long)a->parent->vidx - (long)b->parent->vidx;
}

static void vcfr_print_entry(VcfReader *vfr, htsFile *outfh, bcf_hdr_t *outhdr,
                             VcfCovAlt **alleles, size_t nalleles)
{
  VcfCovLine *var = alleles[0]->parent;
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
      num_alts_no_covg++;
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
        num_alts_with_covg++;
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
  size_t num_a, alen = vc_alts_len(&vr->aprint);
  if(alen == 0 || (!force && alen < PRINT_BUF_LIMIT)) return;

  VcfCovAlt **alleleptr = vc_alts_getptr(&vr->aprint, 0);
  qsort(alleleptr, alen, sizeof(VcfCovAlt*), _vcfcov_alt_cmp_vidx);

  while(alen > 0)
  {
    alleleptr = vc_alts_getptr(&vr->aprint, 0);
    alen = vc_alts_len(&vr->aprint);

    num_a = 0;
    while(num_a < alen && alleleptr[num_a]->parent->vidx == vr->nxtprint)
      num_a++;

    if(num_a < alleleptr[0]->parent->nchildren) break;
    ctx_assert(num_a == alleleptr[0]->parent->nchildren && num_a > 0);

    // print
    vcfr_print_entry(vr, outfh, outhdr, alleleptr, num_a);

    // Re-add to vpool
    vc_lines_append(&vr->vpool, alleleptr[0]->parent);

    // Re-add to apool
    vc_alts_unshift(&vr->apool, alleleptr, num_a);
    vc_alts_shift(&vr->aprint, NULL, num_a);

    vr->nxtprint++;
    alen -= num_a;
  }
}

static void vcfcov_vars(VcfCovAlt **vars, size_t nvars,
                        size_t tgtidx, size_t ntgts,
                        const char *chrom, size_t chromlen,
                        Genotyper *gtyper, const dBGraph *db_graph)
{
  ctx_assert(ntgts <= nvars);

  HaploKmer *kmers = NULL;
  size_t i, end, nkmers;

  if(nvars > max_gt_vars) { return; }

  nkmers = genotyping_get_kmers(gtyper, (const VcfCovAlt *const*)vars,
                                nvars, tgtidx, ntgts,
                                chrom, chromlen,
                                db_graph->kmer_size, &kmers);

  for(i = 0; i < nkmers; i++) {
    bkey_get_covg(kmers[i].bkey, kmers[i].arbits,
                  vars+tgtidx, ntgts,
                  db_graph);
  }

  for(i = tgtidx, end = tgtidx+ntgts; i < end; i++) {
    ctx_assert(!vars[i]->has_covg);
    vars[i]->has_covg = true;
  }
}

// Get index of first of vars (starting at index i), which starts at/after endpos
// otherwise return nvars
static inline int vc_alts_starts_after(VcfCovAlt **vars, size_t nvars,
                                      size_t i, size_t endpos)
{
  while(i < nvars && vars[i]->pos < endpos) { i++; }
  return i;
}

// Get index of first of vars (starting at index i), which ends at/after endpos
// otherwise return nvars
static inline int vc_alts_ends_after(VcfCovAlt **vars, size_t nvars,
                                     size_t i, size_t endpos, size_t kmer_size)
{
  while(i < nvars && vcfcov_alt_end(vars[i])+kmer_size <= endpos) { i++; }
  return i;
}

// This will replace vcfcov_block
static void vcfcov_block2(VcfCovAlt **vars, size_t nvars,
                          size_t tgtidx, size_t ntgts,
                          const char *chrom, int chromlen,
                          Genotyper *gtyper, const dBGraph *db_graph)
{
  ctx_assert(tgtidx+ntgts<=nvars);
  if(!ntgts) { return; }
  else if(nvars <= max_gt_vars)
  {
    vcfcov_vars(vars, nvars, tgtidx, ntgts, chrom, chromlen, gtyper, db_graph);
  }
  else
  {
    // do a few at a time
    const size_t ks = db_graph->kmer_size;
    // genotype start/end, background start/end (end is not inclusive)
    size_t i, gs = tgtidx, ge, bs, be, tmp_ge, tmp_be;

    while(gs < tgtidx+ntgts)
    {
      // Get vars to the left of our genotyping-start (gs)
      for(i = bs = gs; i > 0; i--) {
        if(vcfcov_alt_end(vars[i-1]) + ks > vars[gs]->pos) {
          --bs; // don't deincrement bs in SWAP macro
          SWAP(vars[i-1], vars[bs]);
        }
      }

      ge = gs+1;
      be = vc_alts_starts_after(vars, nvars, ge, vcfcov_alt_end(vars[ge-1])+ks);

      // Try increasing ge, calculate be, accept if small enough
      for(tmp_ge = ge+1; tmp_ge < tgtidx+ntgts; tmp_ge++) {
        tmp_be = vc_alts_starts_after(vars, nvars, tmp_ge,
                                      vcfcov_alt_end(vars[tmp_ge-1])+ks);
        if(tmp_be - bs <= max_gt_vars) { ge = tmp_ge; be = tmp_be; }
        else break;
      }

      ctx_assert2(bs<=gs && gs<ge && ge<=be, "%zu %zu %zu %zu",bs,gs,ge,be);
      ctx_assert2(ge<=tgtidx+ntgts, "%zu %zu %zu %zu",ge,tgtidx,ntgts,nvars);
      ctx_assert2(be<=nvars, "%zu %zu %zu %zu",be,tgtidx,ntgts,nvars);

      // status("bs:%zu gs:%zu ge:%zu be:%zu", bs, gs, ge, be);
      vcfcov_vars(vars+bs, be-bs, gs-bs, ge-gs,
                  chrom, chromlen, gtyper, db_graph);

      gs = ge;
    }
  }
}

// return number of alts that have been genotyped but not removed
static size_t vcfcov_block3(VcfReader *vr, bool flush, size_t tgtidx,
                            const char *chr, int chrlen,
                            Genotyper *gtyper)
{
  const dBGraph *db_graph = vr->db_graph;
  const size_t ks = db_graph->kmer_size;

  VcfCovAlt **vars = vc_alts_getptr(&vr->alist, 0);
  size_t nvars = vc_alts_len(&vr->alist);
  ctx_assert(nvars > 0);

  // 0. get end pos and last index -> set to end if flush
  size_t lastpos = vr->curr_line->v.pos;
  size_t lastidx = flush ? nvars : vc_alts_ends_after(vars, nvars, 0, lastpos, ks);
  size_t bs, be, gs, ge, endpos = 0;

  // 1. Find breaks of >=kmer_size -> break into blocks
  for(bs=0, gs=tgtidx, ge=gs+1; ge <= lastidx; ge++)
  {
    endpos = MAX2(endpos, vcfcov_alt_end(vars[ge-1]) + ks);
    if(endpos <= vars[ge-1]->pos) {
      // end of block
      be = ge;
      vcfcov_block2(vars+bs, be-bs, gs-bs, ge-gs,
                    chr, chrlen, gtyper, db_graph);
      bs = gs = ge;
    }
  }

  // 2. Pass remainder as block
  be = nvars;
  ge = lastidx;
  vcfcov_block2(vars+bs, be-bs, gs-bs, ge-gs, chr, chrlen, gtyper, db_graph);

  // 3. Find start of background required for next time
  size_t i, nremove = nvars;
  if(lastidx < nvars)
    nremove = vc_alts_ends_after(vars, nvars, 0, vars[lastidx]->pos, ks);
  for(i = 0; i < nremove; i++) vcfr_move_to_print(vr, i);

  // return number of un-genotyped variants
  return lastidx - nremove;
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
  char *chr = NULL;
  int refid = -1, chrlen = 0;

  int n;
  size_t nrem = 0;

  while((n = vcfr_fetch(&vr)) >= 0)
  {
    // Get ref chromosome
    // Only loads if we don't currently have the right chrom
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chr, &chrlen);

    nrem = vcfcov_block3(&vr, n == 0, nrem, chr, chrlen, gtyper);

    // Shrink array if we dropped any vars
    vcfr_repack_alts(&vr);
    vcfr_print_waiting(&vr, outfh, outhdr, false);
  }

  // Deal with remainder
  if(vc_alts_len(&vr.alist) > 0) {
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chr, &chrlen);
    nrem = vcfcov_block3(&vr, true, nrem, chr, chrlen, gtyper);
  }
  vcfr_repack_alts(&vr);
  vcfr_print_waiting(&vr, outfh, outhdr, true);

  free(chr);
  genotyper_destroy(gtyper);
  vcfr_dealloc(&vr);
}

int ctx_vcfcov(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL, *out_type = NULL;

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
      case 'O': cmd_check(!out_type, cmd); out_type = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'r': cmd_check(!ref_path, cmd); ref_path = optarg; break;
      case 'L': cmd_check(!max_allele_len,cmd); max_allele_len = cmd_size(cmd,optarg); break;
      case 'N': cmd_check(!max_gt_vars,cmd); max_gt_vars = cmd_uint32(cmd,optarg); break;
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

  if(!max_allele_len) max_allele_len = DEFAULT_MAX_ALLELE_LEN;
  if(!max_gt_vars) max_gt_vars = DEFAULT_MAX_GT_VARS;

  status("[vcfcov] max allele length: %zu; max number of variants: %u",
         max_allele_len, max_gt_vars);

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
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  // v=>vcf, z=>compressed vcf, b=>bcf, bu=>uncompressed bcf
  const char *modes[] = {"wv","wz","wbu","wb"};
  const char *hsmodes[] = {"uncompressed VCF", "compressed VCF",
                           "uncompressed BCF", "compressed BCF"};
  int mode = 1;

  if(out_type) {
    const char *ot = out_type;
    if(!strcmp(ot,"vcf") || !strcmp(ot,"uvcf") || !strcmp(ot,"v")) mode = 0;
    else if(!strcmp(ot,"vcfgz") || !strcmp(ot,"z")) mode = 1;
    else if(!strcmp(ot,"ubcf")  || !strcmp(ot,"u")) mode = 2;
    else if(!strcmp(ot,"bcf")   || !strcmp(ot,"b")) mode = 3;
  } else {
    if(     futil_path_has_extension(out_path,".vcf"))    mode = 0;
    else if(futil_path_has_extension(out_path,".vcfgz"))  mode = 1;
    else if(futil_path_has_extension(out_path,".vcf.gz")) mode = 1;
    else if(futil_path_has_extension(out_path,".ubcf"))   mode = 2;
    else if(futil_path_has_extension(out_path,".bcf"))    mode = 3;
  }

  futil_create_output(out_path);
  htsFile *outfh = hts_open(out_path, modes[mode]);
  status("[vcfcov] Output format: %s", hsmodes[mode]);


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
    status("[vcfcov] Colour %zu: %s [VCF column %zu]", i, sname, samplehdrids[i]);
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

  status("[vcfcov]  saved to: %s\n", out_path);

  // Print statistics
  char ns0[50], ns1[50];
  status("[vcfcov] Read %s VCF lines", ulong_to_str(num_vcf_lines, ns0));
  status("[vcfcov] Read %s ALTs", ulong_to_str(num_alts_read, ns0));
  status("[vcfcov] ALTs used: %s / %s (%.2f%%)",
         ulong_to_str(num_alts_loaded, ns0), ulong_to_str(num_alts_read, ns1),
         num_alts_read ? (100.0*num_alts_loaded) / num_alts_read : 0.0);
  status("[vcfcov] ALTs too long (>%zubp): %s / %s (%.2f%%)", max_allele_len,
         ulong_to_str(num_alts_too_long, ns0), ulong_to_str(num_alts_read, ns1),
         num_alts_read ? (100.0*num_alts_too_long) / num_alts_read : 0.0);
  status("[vcfcov] ALTs too dense (>%u within %zubp): %s / %s (%.2f%%)",
         max_gt_vars, db_graph.kmer_size,
         ulong_to_str(num_alts_no_covg, ns0), ulong_to_str(num_alts_read, ns1),
         num_alts_read ? (100.0*num_alts_no_covg) / num_alts_read : 0.0);
  status("[vcfcov] ALTs printed with coverage: %s / %s (%.2f%%)",
         ulong_to_str(num_alts_with_covg, ns0), ulong_to_str(num_alts_read, ns1),
         num_alts_read ? (100.0*num_alts_with_covg) / num_alts_read : 0.0);

  ctx_free(samplehdrids);

  bcf_hdr_destroy(vcfhdr);
  bcf_hdr_destroy(outhdr);
  hts_close(vcffh);
  hts_close(outfh);
  fai_destroy(fai);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
