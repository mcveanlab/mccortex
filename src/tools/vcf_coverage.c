#include "global.h"
#include "vcf_coverage.h"
#include "util.h"
#include "db_node.h"
#include "genotyping.h"
#include "vcf_misc.h"

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(covg_buf, CovgBuffer, Covg);

// Initial object buffer lengths
#define INIT_BUF_SIZE 128

// How many alt alleles to collect before printing
#define PRINT_BUF_LIMIT 100

// for debugging
#ifdef DEBUG_VCFCOV
  bcf_hdr_t *globalhdr;
#endif

typedef struct
{
  const char *path;
  const size_t *samplehdrids; // samplehdrids[x] is the sample id in VCF
  size_t ncols, kmer_size; // number of colours loaded in the de Bruijn graph
  VcfCovStats stats;
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
  // Format for VCF output
  int32_t *kcovgs_r, *kcovgs_a;
  size_t geno_buf_size; // nsamples * nalts
} VcfReader;

// Genotyping buffers
typedef struct
{
  Genotyper *gtyper;
  // Fetch coverage from the graph
  CovgBuffer *covgs;
  size_t clen; // number of buffer is nalts*ncols*2 (2=>ref/alt for each alt)
  Uint32Buffer nrkmers;
  // Graph to get kmer coverage from or add kmers to
  dBGraph *db_graph;
} VcfCovBuffers;

static void covbuf_alloc(VcfCovBuffers *covbuf, dBGraph *db_graph)
{
  memset(covbuf, 0, sizeof(*covbuf));
  covbuf->db_graph = db_graph;
  covbuf->gtyper = genotyper_init();
  uint32_buf_alloc(&covbuf->nrkmers, 16);
}

static void covbuf_dealloc(VcfCovBuffers *covbuf)
{
  size_t i;
  for(i = 0; i < covbuf->clen; i++) covg_buf_dealloc(&covbuf->covgs[i]);
  ctx_free(covbuf->covgs);
  genotyper_destroy(covbuf->gtyper);
  uint32_buf_dealloc(&covbuf->nrkmers);
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
                              size_t ncols, size_t kmer_size)
{
  memset(vcfr, 0, sizeof(*vcfr));
  vcfr->path = path;
  vcfr->samplehdrids = samplehdrids;
  vcfr->ncols = ncols;
  vcfr->kmer_size = kmer_size;
  vcfr->vcffh = vcffh;
  vcfr->vcfhdr = vcfhdr;

  vc_lines_alloc(&vcfr->vpool, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->alist, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->anchrom, INIT_BUF_SIZE);
  vc_alts_alloc(&vcfr->apool, INIT_BUF_SIZE);

  vcf_list_populate(&vcfr->vpool, INIT_BUF_SIZE);
  var_list_populate(&vcfr->apool, INIT_BUF_SIZE, ncols);
}

static inline void vcfr_dealloc(VcfReader *vcfr)
{
  ctx_assert(vc_alts_len(&vcfr->alist) == 0);
  ctx_assert(vc_alts_len(&vcfr->anchrom) == 0);
  ctx_assert(vc_alts_len(&vcfr->aprint) == 0);

  vcf_list_destroy(&vcfr->vpool);
  var_list_destroy(&vcfr->apool);
  ctx_assert(vc_alts_len(&vcfr->apool) == 0);
  ctx_assert(vc_lines_len(&vcfr->vpool) == 0);

  vc_alts_dealloc(&vcfr->alist);
  vc_alts_dealloc(&vcfr->anchrom);
  vc_alts_dealloc(&vcfr->aprint);
  vc_alts_dealloc(&vcfr->apool);
  vc_lines_dealloc(&vcfr->vpool);

  ctx_free(vcfr->kcovgs_r);
  ctx_free(vcfr->kcovgs_a);
}

// return true if valid variant
static inline void init_new_alt(VcfCovAlt *var, VcfCovLine *line, uint32_t aid,
                                size_t ncols)
{
  bcf1_t *v = &line->v;
  size_t reflen = 0, altlen = 0, shift;
  const char *ref = v->d.allele[0], *alt = v->d.allele[aid];

  // Trim bases that match with the ref
  shift = trimmed_alt_lengths(v, aid, &reflen, &altlen);

  // Initialise
  var->parent = line;
  var->ref = ref + shift;
  var->alt = alt + shift;
  var->pos = v->pos + shift;
  var->reflen = reflen;
  var->altlen = altlen;
  var->aid = aid;
  vcfcov_alt_wipe_covg(var, ncols);
}

static inline void fetch_chrom(bcf_hdr_t *hdr, bcf1_t *v,
                               faidx_t *fai, int *refid,
                               char **chr, int *chrlen)
{
  if(*refid != v->rid) {
    free(*chr);
    *chr = fai_fetch(fai, bcf_seqname(hdr, v), chrlen);
    if(*chr == NULL) die("Cannot find chr '%s'", bcf_seqname(hdr, v));
    *refid = v->rid;
  }
  // Check ref allele is within the reference and matches
  if(v->pos + v->rlen > *chrlen) {
    die("Ref allele goes out of bounds: %s %i %s %s [chrlen: %i]",
        bcf_seqname(hdr, v), v->pos+1, v->d.id, v->d.allele[0], *chrlen);
  }
  if(strncasecmp((*chr)+v->pos,v->d.allele[0],v->rlen) != 0) {
    die("Alleles don't match: %s %i %s %s; ref: %.*s [%i]",
        bcf_seqname(hdr, v), v->pos+1, v->d.id, v->d.allele[0],
        v->rlen, *chr+v->pos, v->rlen);
  }
}

// Returns:
//  -1 if EOF / error
//  0 if read would have added var on another chrom
//  Otherwise number of variants added
//
static int vcfr_fetch(VcfReader *vr, const VcfCovPrefs *prefs)
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
  bcf1_t *v;

  while(1)
  {
    // Get VcfCovLine to read into
    vc_lines_pop(&vr->vpool, &line, 1);
    line->vidx = vr->nextidx++;
    v = &line->v;

    // Read VCF
    if(bcf_read(vr->vcffh, vr->vcfhdr, v) < 0) {
      // EOF
      vc_lines_push(&vr->vpool, &line, 1);
      return -1;
    }

    // Unpack all info
    bcf_unpack(v, BCF_UN_ALL);
    vr->stats.nvcf_lines++;
    vr->curr_line = line;

    ctx_assert(strlen(v->d.allele[0]) == (size_t)v->rlen);
    ctx_assert2(v->n_allele > 1, "n_allele: %i", v->n_allele);

    // Check we have enough vars to decompose
    size_t i, n = MAX2(v->n_allele, 16);
    if(vc_alts_len(&vr->apool) < n)
      var_list_populate(&vr->apool, n, vr->ncols);

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
    VcfCovAltPtrList *list = diff_chroms ? &vr->anchrom : &vr->alist;

    // i==0 is ref allele
    vc_alts_pop(&vr->apool, &alt, 1);
    for(i = 1; i < v->n_allele; i++) {
      init_new_alt(alt, line, i, vr->ncols);
      // straight to print queue if not valid
      if(MAX2(alt->reflen, alt->altlen) > prefs->max_allele_len) {
        vr->stats.nalts_too_long++;
        vc_alts_append(&vr->aprint, alt);
      } else if(alt->reflen == 0 && alt->altlen == 0) {
        vc_alts_append(&vr->aprint, alt);
      } else {
        vc_alts_append(list, alt);
        nadded++;
      }
      vc_alts_pop(&vr->apool, &alt, 1);
    }
    // Re-add unused var back into pool
    vc_alts_push(&vr->apool, &alt, 1);

    vr->stats.nalts_loaded += nadded;
    vr->stats.nalts_read += v->n_allele - 1; // v_allele includes ref

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

// Sort by variant index then by allele index
static int _vcfcov_alt_cmp_vidx(const void *aa, const void *bb)
{
  const VcfCovAlt *a = *(const VcfCovAlt*const*)aa;
  const VcfCovAlt *b = *(const VcfCovAlt*const*)bb;
  int c = cmp(a->parent->vidx, b->parent->vidx);
  return c ? c : cmp(a->aid, b->aid);
}

// alleles should be sorted by parent->vidx, then by aid
//  [See _vcfcov_alt_cmp_vidx(a,b).]
static void vcfr_print_entry(VcfReader *vr, htsFile *outfh, bcf_hdr_t *outhdr,
                             VcfCovAlt **alleles, size_t nalleles,
                             const VcfCovPrefs *prefs)
{
  VcfCovLine *var = alleles[0]->parent;
  bcf1_t *v = &var->v;
  size_t nsamples = bcf_hdr_nsamples(outhdr);
  size_t i, col, sid, nalts = v->n_allele-1, n = nsamples * nalts;
  size_t ncols = vr->ncols, nalts_covgs = 0;
  VarCovg *cov;

  ctx_assert2(nalts == nalleles, "%zu vs %zu", nalts, nalleles);

  // _r ref, _a alt
  if(vr->geno_buf_size < n)
  {
    vr->kcovgs_r = ctx_reallocarray(vr->kcovgs_r, n, sizeof(int32_t));
    vr->kcovgs_a = ctx_reallocarray(vr->kcovgs_a, n, sizeof(int32_t));
    for(i = vr->geno_buf_size; i < n; i++)
      vr->kcovgs_r[i] = vr->kcovgs_a[i] = bcf_int32_missing;
    vr->geno_buf_size = n;
  }

  // Fetch existing coverage from VCF
  if(nsamples > ncols)
  {
    int nsize = vr->geno_buf_size;
    bcf_get_format_int32(vr->vcfhdr, v, prefs->kcov_ref_tag, &vr->kcovgs_r, &nsize);
    bcf_get_format_int32(vr->vcfhdr, v, prefs->kcov_alt_tag, &vr->kcovgs_a, &nsize);
    ctx_assert2(nsize == (int)vr->geno_buf_size, "htslib resized our buffer!");
  }

  // Initiate new samples to missing
  for(i = bcf_hdr_nsamples(vr->vcfhdr)*nalts; i < n; i++)
    vr->kcovgs_r[i] = vr->kcovgs_a[i] = bcf_int32_missing;

  // Add coverage
  for(i = 0; i < nalts; i++) {
    if(alleles[i]->has_covg) {
      nalts_covgs++;
      for(col = 0; col < ncols; col++) {
        sid = vr->samplehdrids[col];
        cov = &alleles[i]->c[col];
        vr->kcovgs_r[sid*nalts+i] = cov->covg[0];
        vr->kcovgs_a[sid*nalts+i] = cov->covg[1];
      }
    }
  }

  // Update stats
  vr->stats.nalts_no_covg += nalts - nalts_covgs;
  vr->stats.nalts_with_covg += nalts_covgs;

  //
  // Update VCF entry
  //

  // TODO: reset sample fields on new samples
  for(i = bcf_hdr_nsamples(vr->vcfhdr); i < nsamples; i++) {
    /* reset fields */
  }

  // Update sample info
  int a,b;
  a = bcf_update_format_int32(outhdr, v, prefs->kcov_ref_tag, vr->kcovgs_r, n);
  b = bcf_update_format_int32(outhdr, v, prefs->kcov_alt_tag, vr->kcovgs_a, n);

  if(a || b) die("Cannot add format info");
  if(bcf_write(outfh, outhdr, v) != 0) die("Cannot write record");
}

static void vcfr_print_waiting(VcfReader *vr,
                               htsFile *outfh, bcf_hdr_t *outhdr,
                               const bcf_hdr_t *vcfhdr,
                               const VcfCovPrefs *prefs, bool force)
{
  // Sort waiting by vidx
  size_t num_a, alen = vc_alts_len(&vr->aprint);
  if(alen == 0 || (!force && alen < PRINT_BUF_LIMIT)) return;

  // Input header may have been modified by reading an entry
  // for instance adding a missing contig= entry
  bcf_hdr_merge(outhdr, vcfhdr);

  VcfCovLine *line;
  VcfCovAlt **alleleptr = vc_alts_getptr(&vr->aprint, 0);
  qsort(alleleptr, alen, sizeof(VcfCovAlt*), _vcfcov_alt_cmp_vidx);

  while(alen > 0)
  {
    alleleptr = vc_alts_getptr(&vr->aprint, 0);
    alen = vc_alts_len(&vr->aprint);

    num_a = 0;
    while(num_a < alen && alleleptr[num_a]->parent->vidx == vr->nxtprint)
      num_a++;

    line = alleleptr[0]->parent;
    if(num_a+1 < line->v.n_allele) break;
    ctx_assert(num_a+1 == line->v.n_allele && num_a > 0);

    // print
    if(!prefs->load_kmers_only)
      vcfr_print_entry(vr, outfh, outhdr, alleleptr, num_a, prefs);

    // Re-add to vpool
    vc_lines_append(&vr->vpool, line);

    // Re-add to apool
    vc_alts_unshift(&vr->apool, alleleptr, num_a);
    vc_alts_shift(&vr->aprint, NULL, num_a);

    vr->nxtprint++;
    alen -= num_a;
  }
}

#define covgbufidx(var,col,ncols,isalt) ((var)*(ncols)*2 + (col)*2 + (isalt))

// Returns new length
static inline void resize_covg_bufs(CovgBuffer **covgsp, size_t *lenp,
                                    size_t ncols, size_t nvars, size_t nkmers)
{
  // We have a buffer for each variant, in each colour, in ref and alt
  size_t i, nlen = nvars*ncols*2;

  // Resize existing buffers
  if(*lenp && (*covgsp)[0].size < nkmers) {
    for(i = 0; i < *lenp; i++)
      covg_buf_capacity(&(*covgsp)[i], nkmers);
  }

  if(nlen > *lenp) {
    *covgsp = ctx_realloc(*covgsp, sizeof((*covgsp)[0]) * nlen);
    for(i = *lenp; i < nlen; i++)
      covg_buf_alloc(&(*covgsp)[i], nkmers);
    *lenp = nlen;
  }
}

/**
 * Get coverage of ref and alt alleles from the de Bruijn graph in the given
 * colour.
 */
static inline void bkey_get_covg(BinaryKmer bkey,
                                 uint64_t altref_bits, size_t ntgts,
                                 CovgBuffer *covgs, // covgs[nvar*ncols*2]
                                 const dBGraph *db_graph)
{
  size_t i, col, ncols = db_graph->num_of_cols;
  dBNode node = db_graph_find(db_graph, bkey);
  uint64_t arbits;
  Covg covg;

  ctx_assert(altref_bits);

  if(node.key != HASH_NOT_FOUND) {
    // printf("node: ");
    // db_nodes_print(&node, 1, db_graph, stdout);
    // printf(" bits:%#06x\n", (unsigned int)altref_bits);

    for(col = 0; col < ncols; col++) {
      covg = db_node_get_covg(db_graph, node.key, col);
      // printf("  col%zu: %u\n", col, covg);
      if(!covg) continue;

      for(i = 0, arbits = altref_bits; i < ntgts; i++, arbits >>= 2) {
        switch(arbits & 3UL) {
          case 0: break; /* ignore kmer in neither ref nor alt */
          case 1: covg_buf_add(&covgs[covgbufidx(i,col,ncols,0)], covg); break;
          case 2: covg_buf_add(&covgs[covgbufidx(i,col,ncols,1)], covg); break;
          case 3: break; /* ignore kmer in both ref and alt */
        }
      }
    }
  }
}

// +0.5 to round correctly
#define vmeancovg(tot,nk) ((nk) ? ((double)(tot)) / (nk) + 0.5 : bcf_int32_missing)

/**
 * @param nkrmers  for each variant, lists number of kmers possible from ref
 */
static void vcfcov_update_covg(const HaploKmer *kmers, size_t nkmers,
                               VcfCovAlt **vars, size_t nvars,
                               const uint32_t *nrkmers,
                               CovgBuffer *covgs,
                               const dBGraph *db_graph)
{
  VcfCovAlt *var;
  CovgBuffer *rcovg, *acovg;
  size_t i, j, rk, ak, col, ncols = db_graph->num_of_cols;
  uint64_t rtot, atot;

  for(i = 0; i < nkmers; i++) {
    bkey_get_covg(kmers[i].bkey, kmers[i].arbits, nvars,
                  covgs, db_graph);
  }

  for(i = 0; i < nvars; i++)
  {
    var = vars[i];
    ctx_assert(!var->has_covg);

    // rk = hap_num_exp_kmers(var->pos, var->reflen, kmer_size);
    // ak = hap_num_exp_kmers(var->pos, var->altlen, kmer_size);
    rk = nrkmers[i];
    ak = vcfcovalt_akmers(var,nrkmers[i]);

    // status("nrkmers: %zu/%zu vs r+ks-1: %zu/%zu",
    //        (size_t)nrkmers[i], vcfcovalt_akmers(var,nrkmers[i]), rk, ak);

    for(col = 0; col < ncols; col++)
    {
      rcovg = &covgs[covgbufidx(i,col,ncols,0)];
      acovg = &covgs[covgbufidx(i,col,ncols,1)];
      for(j = rtot = 0; j < rcovg->len; j++) rtot += rcovg->b[j];
      for(j = atot = 0; j < acovg->len; j++) atot += acovg->b[j];

#ifdef DEBUG_VCFCOV
      printf("ref:");
      for(j = 0; j < rcovg->len; j++) printf(" %u", rcovg->b[j]);
      printf(" = %zu / %zu = %f\n", (size_t)rtot, rk, vmeancovg(rtot, rk));
      printf("alt:");
      for(j = 0; j < acovg->len; j++) printf(" %u", acovg->b[j]);
      printf(" = %zu / %zu = %f\n", (size_t)atot, ak, vmeancovg(atot, ak));
#endif

      // rkmers/akmers is est. of num. of ref/alt kmers
      var->c[col].covg[0] = vmeancovg(rtot, rk);
      var->c[col].covg[1] = vmeancovg(atot, ak);
    }
    var->has_covg = true;
  }
}

static void vcfcov_vars(VcfCovAlt **vars, size_t nvars,
                        size_t tgtidx, size_t ntgts,
                        const char *chrom, size_t chromlen,
                        VcfCovBuffers *covbuf,
                        const VcfCovPrefs *prefs, VcfCovStats *stats)
{
  ctx_assert(ntgts <= nvars);

  HaploKmer *kmers = NULL;
  size_t i, nkmers;
  dBGraph *db_graph = covbuf->db_graph;

  if(nvars > prefs->max_gt_vars) { return; }

#ifdef DEBUG_VCFCOV
  // debug printing
  kstring_t s = {0,0,NULL};
  printf("  %zu-%zu / %zu\n", tgtidx, tgtidx+ntgts-1, nvars);
  for(i = 0; i < nvars; i++) {
    vcf_format(globalhdr, &vars[i]->parent->v, &s);
    fprintf(stdout, "%zu: %s", i, s.s);
    s.s[s.l=0] = 0; // reset kstr
  }
  free(s.s);
#endif

  // Number of ref kmers for each variant
  uint32_buf_capacity(&covbuf->nrkmers, ntgts);
  uint32_t *nrkmers = covbuf->nrkmers.b;

  // returns kmer keys
  nkmers = genotyping_get_kmers(covbuf->gtyper, (const VcfCovAlt *const*)vars,
                                nvars, tgtidx, ntgts,
                                chrom, chromlen,
                                db_graph->kmer_size,
                                &kmers, nrkmers);

  stats->ngt_kmers += nkmers;

  if(prefs->load_kmers_only)
  {
    // Add kmers to graph
    // don't need binary_kmer_get_key(), genotyping returns keys only
    bool found = false;
    for(i = 0; i < nkmers; i++) {
      hash_table_find_or_insert(&db_graph->ht, kmers[i].bkey, &found);

#ifdef DEBUG_VCFCOV
      char tmpstr[MAX_KMER_SIZE+1], binstr[65];
      fprintf(stdout, "=%s %s\n",
              binary_kmer_to_str(kmers[i].bkey, db_graph->kmer_size, tmpstr),
              bin64_to_str(kmers[i].arbits, ntgts*2, binstr));
#endif
    }
  }
  else
  {
    // Reset coverage buffers
    for(i = 0; i < covbuf->clen; i++) covbuf->covgs[i].len = 0;
    resize_covg_bufs(&covbuf->covgs, &covbuf->clen,
                     db_graph->num_of_cols, nvars, nkmers);

    // Add coverage to vars
    vcfcov_update_covg(kmers, nkmers, vars + tgtidx, ntgts, nrkmers,
                       covbuf->covgs, db_graph);
  }
}

// Get index of first of vars (starting at index i), which starts at/after endpos
// otherwise return nvars
static inline size_t vc_alts_starts_after(VcfCovAlt **vars, size_t nvars,
                                          size_t i, size_t endpos)
{
  while(i < nvars && vars[i]->pos < endpos) { i++; }
  return i;
}

// Get index of first of vars (starting at index i), whose haplotype includes
// endpos. Otherwise return nvars
static inline size_t vc_alts_ends_after(VcfCovAlt **vars, size_t nvars,
                                        size_t i, size_t endpos, size_t kmer_size)
{
  while(i < nvars && vcfcovalt_hap_end(vars[i],kmer_size) <= endpos) { i++; }
  return i;
}

static void vcfcov_block(VcfCovAlt **vars, size_t nvars,
                         size_t tgtidx, size_t ntgts,
                         const char *chrom, int chromlen,
                         VcfCovBuffers *covbuf,
                         const VcfCovPrefs *prefs, VcfCovStats *stats)
{
  // printf("nvars: %zu tgtidx: %zu ntgts: %zu\n", nvars, tgtidx, ntgts);

  ctx_assert(tgtidx+ntgts<=nvars);
  if(!ntgts) { return; }
  else if(nvars <= prefs->max_gt_vars)
  {
    vcfcov_vars(vars, nvars, tgtidx, ntgts, chrom, chromlen, covbuf, prefs, stats);
  }
  else
  {
    // do a few at a time
    const size_t ks = covbuf->db_graph->kmer_size;
    // genotype start/end, background start/end (end is not inclusive)
    size_t i, gs = tgtidx, ge, bs, be, tmp_ge, tmp_be, endpos;

    while(gs < tgtidx+ntgts)
    {
      // Get vars to the left of our genotyping-start (gs)
      for(i = bs = gs; i > 0; i--) {
        if(vcfcovalt_hap_end(vars[i-1],ks) > vars[gs]->pos) {
          --bs; // don't deincrement bs in SWAP macro
          SWAP(vars[i-1], vars[bs]);
        }
      }

      ge = gs+1;
      endpos = vcfcovalt_hap_end(vars[ge-1],ks);
      be = vc_alts_starts_after(vars, nvars, ge, endpos);

      // Try increasing ge, calculate be, accept if small enough
      for(tmp_ge = ge+1; tmp_ge < tgtidx+ntgts; tmp_ge++) {
        endpos = MAX2(endpos, vcfcovalt_hap_end(vars[tmp_ge-1],ks));
        tmp_be = vc_alts_starts_after(vars, nvars, tmp_ge, endpos);
        if(tmp_be - bs <= prefs->max_gt_vars) { ge = tmp_ge; be = tmp_be; }
        else break;
      }

      ctx_assert2(bs<=gs && gs<ge && ge<=be, "%zu %zu %zu %zu",bs,gs,ge,be);
      ctx_assert2(ge<=tgtidx+ntgts, "%zu %zu %zu %zu",ge,tgtidx,ntgts,nvars);
      ctx_assert2(be<=nvars, "%zu %zu %zu %zu",be,tgtidx,ntgts,nvars);

      // status("bs:%zu gs:%zu ge:%zu be:%zu", bs, gs, ge, be);
      vcfcov_vars(vars+bs, be-bs, gs-bs, ge-gs,
                  chrom, chromlen, covbuf, prefs, stats);

      gs = ge;
    }
  }
}

// return number of alts that have been genotyped but not removed
static size_t vcfcov_block2(VcfReader *vr, bool flush, size_t tgtidx,
                            const char *chr, int chrlen,
                            VcfCovBuffers *covbuf,
                            const VcfCovPrefs *prefs)
{
  const size_t ks = vr->kmer_size;

  VcfCovAlt **vars = vc_alts_getptr(&vr->alist, 0);
  size_t nvars = vc_alts_len(&vr->alist);
  ctx_assert(nvars > 0);

  // 0. get end pos and last index -> set to end if flush
  size_t lastpos = vr->curr_line->v.pos;
  size_t lastidx = flush ? nvars-1 : vc_alts_ends_after(vars, nvars, 0, lastpos, ks);
  size_t bs, be, gs, ge, endpos = 0;

  // if !flush, variant at lastidx will be kept.
  // between lastidx-1,lastidx is the last position to check for a gap of k bp
  ctx_assert(lastidx < nvars);
  // status("lastidx: %zu / %zu; lastpos: %zu", lastidx, nvars, lastpos);

  // 1. Find breaks of >=kmer_size -> break into blocks
  for(bs=0, gs=tgtidx, ge=gs+1; ge <= lastidx; ge++)
  {
    endpos = MAX2(endpos, vcfcovalt_hap_end(vars[ge-1],ks));
    if(endpos <= vars[ge]->pos) {
      // end of block
      be = ge;
      vcfcov_block(vars+bs, be-bs, gs-bs, ge-gs,
                   chr, chrlen, covbuf, prefs, &vr->stats);
      bs = gs = ge;
    }
  }

  // 2. Pass remainder as block
  be = nvars;
  ge = flush ? nvars : lastidx;
  // printf("bs: %zu-%zu gs: %zu-%zu lastidx: %zu\n", bs, be, gs, ge, lastidx);
  vcfcov_block(vars+bs, be-bs, gs-bs, ge-gs,
               chr, chrlen, covbuf, prefs, &vr->stats);

  // 3. Find start of background required for next time
  size_t i, j, nxttgt = 0, nxtpos;
  VcfCovAlt *alt;

  if(ge == nvars) {
    // Move all variants from alist -> print list
    vc_alts_push(&vr->aprint, vc_alts_getptr(&vr->alist, 0), nvars);
    vc_alts_shift(&vr->alist, NULL, nvars);
    nxttgt = 0;
  }
  else {
    // Move only those that come before the next variant to be genotyped
    nxtpos = vars[ge]->pos;
    for(i = j = 0; i < nvars; i++) {
      alt = vc_alts_get(&vr->alist, i);
      if(vcfcovalt_hap_end(alt, ks) <= nxtpos) {
        vc_alts_append(&vr->aprint, alt);
      }
      else {
        vc_alts_set(&vr->alist, j, alt);
        j++;
        nxttgt += (i < ge);
      }
    }
    // remove from end
    vc_alts_pop(&vr->alist, NULL, nvars-j);
  }

  // return number of genotyped variants that have not been removed
  return nxttgt;
}

/**
 * @param outfh        Only required for printing
 * @param outhdr       Only required for printing
 * @param samplehdrids [col] is the index of a colour in the output vcf.
 *                     Only required for printing.
 */
void vcfcov_file(htsFile *vcffh, bcf_hdr_t *vcfhdr,
                 htsFile *outfh, bcf_hdr_t *outhdr,
                 const char *path, faidx_t *fai,
                 const size_t *samplehdrids,
                 const VcfCovPrefs *prefs,
                 VcfCovStats *stats,
                 dBGraph *db_graph)
{
  VcfReader vr;
  vcfr_alloc(&vr, path, vcffh, vcfhdr, samplehdrids,
             db_graph->num_of_cols, db_graph->kmer_size);

  // for debugging
#ifdef DEBUG_VCFCOV
  globalhdr = vcfhdr;
#endif

  VcfCovBuffers covbuf;
  covbuf_alloc(&covbuf, db_graph);

  // refid is id of chromosome currently loaded
  char *chr = NULL;
  int refid = -1, chrlen = 0;

  int n;
  size_t tgtidx = 0, max_len = 0;

  // TODO: try reading more than one variant at once

  while((n = vcfr_fetch(&vr, prefs)) >= 0)
  {
    max_len = MAX2(max_len, vc_alts_len(&vr.alist));

    // Get ref chromosome
    // Only loads if we don't currently have the right chrom
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chr, &chrlen);

    tgtidx = vcfcov_block2(&vr, n == 0, tgtidx, chr, chrlen, &covbuf, prefs);

    vcfr_print_waiting(&vr, outfh, outhdr, vcfhdr, prefs, false);
  }

  // Deal with remainder
  if(vc_alts_len(&vr.alist) > 0) {
    fetch_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                &refid, &chr, &chrlen);
    tgtidx = vcfcov_block2(&vr, true, tgtidx, chr, chrlen, &covbuf, prefs);
    ctx_assert(tgtidx == 0);
  }
  vcfr_print_waiting(&vr, outfh, outhdr, vcfhdr, prefs, true);

  status("[vcfcov] max alleles in buffer: %zu", max_len);

  memcpy(stats, &vr.stats, sizeof(*stats));

  free(chr);
  covbuf_dealloc(&covbuf);
  vcfr_dealloc(&vr);
}
