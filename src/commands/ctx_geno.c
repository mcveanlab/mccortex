#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "seq_reader.h"
#include "gpath_checks.h"
#include "seq_reader.h"
#include "genotyping.h"

#include "htslib/vcf.h"
#include "htslib/faidx.h"

// Don't attempt to genotype alleles bigger than this
size_t max_allele_len = 100;

// 2^8 = 256 possible haplotypes
uint32_t max_gt_vars = 8;

// TODO: remove duplicate variants?
// TODO: print VCF with genotype
// TODO: rename vcfcov

const char geno_usage[] =
"usage: "CMD" geno [options] <in.vcf> <in.ctx> [in2.ctx ...]\n"
"\n"
"  Genotype a VCF using cortex graphs. VCF must be sorted by position. \n"
"  VCF must be a file, not piped in. It is recommended to use uncleaned graphs.\n"
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

/**
 * Get coverage of ref and alt alleles from the de Bruijn graph in the given
 * colour.
 */
static inline void bkey_get_covg(BinaryKmer bkey, uint64_t altref_bits,
                                 GenoVar **gts, size_t ntgts,
                                 int colour, const dBGraph *db_graph)
{
  size_t i;
  dBNode node = db_graph_find(db_graph, bkey);

  if(node.key != HASH_NOT_FOUND) {
    Covg covg = colour >= 0 ? db_node_get_covg(db_graph, node.key, colour)
                            : db_node_sum_covg(db_graph, node.key);

    for(i = 0; i < ntgts; i++, altref_bits >>= 2) {
      if((altref_bits & 3) == 1) { gts[i]->nkmers[0]++; gts[i]->sumcovg[0] += covg; }
      else if((altref_bits & 3) == 2) { gts[i]->nkmers[1]++; gts[i]->sumcovg[1] += covg; }
      else { /* ignore kmer in both ref/alt */ }
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

static void var_list_populate(GenoVarPtrList *alist, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) {
    GenoVar *v = ctx_calloc(1, sizeof(GenoVar));
    genovar_ptr_list_append(alist, v);
  }
}

void bcf_empty1(bcf1_t *v);
static void vcf_list_destroy(GenoVCFPtrList *vlist)
{
  size_t i;
  for(i = 0; i < genovcf_ptr_list_len(vlist); i++) {
    GenoVCF *v = genovcf_ptr_list_get(vlist, i);
    bcf_empty1(&v->v);
    ctx_free(v);
  }
}

static void var_list_destroy(GenoVarPtrList *alist)
{
  size_t i;
  for(i = 0; i < genovar_ptr_list_len(alist); i++) {
    GenoVar *a = genovar_ptr_list_get(alist, i);
    ctx_free(a);
  }
}

typedef struct
{
  const char *path;
  htsFile *vcffh;
  bcf_hdr_t *vcfhdr;
  // vpool are ready to be used
  // vwait are waiting to be printed
  GenoVCFPtrList vpool, vwait; // pool of vcf lines
  // index (vidx) of next GenoVCF to be printed
  // used to print VCF entries out in the correct (input) order
  size_t nxtprint, nextidx;
  // alist are current alleles, anchrom are on next chromosome
  // apool is a memory pool
  GenoVarPtrList alist, anchrom, apool; // alleles
  const char *cov_tag, *nkmer_tag;
} VcfReader;

#define INIT_BUF_SIZE 128

static inline void vcfr_alloc(VcfReader *vcfr, const char *path,
                              const char *nkmer_tag, const char *cov_tag,
                              htsFile *vcffh, bcf_hdr_t *vcfhdr)
{
  memset(vcfr, 0, sizeof(*vcfr));
  vcfr->path = path;
  vcfr->nkmer_tag = nkmer_tag;
  vcfr->cov_tag = cov_tag;
  vcfr->vcffh = vcffh;
  vcfr->vcfhdr = vcfhdr;

  genovcf_ptr_list_alloc(&vcfr->vwait, INIT_BUF_SIZE);
  genovcf_ptr_list_alloc(&vcfr->vpool, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->alist, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->anchrom, INIT_BUF_SIZE);
  genovar_ptr_list_alloc(&vcfr->apool, INIT_BUF_SIZE);

  vcf_list_populate(&vcfr->vpool, INIT_BUF_SIZE);
  var_list_populate(&vcfr->apool, INIT_BUF_SIZE);
}

static inline void vcfr_dealloc(VcfReader *vcfr)
{
  ctx_assert(genovcf_ptr_list_len(&vcfr->vwait) == 0);
  ctx_assert(genovar_ptr_list_len(&vcfr->alist) == 0);
  ctx_assert(genovar_ptr_list_len(&vcfr->anchrom) == 0);
  vcf_list_destroy(&vcfr->vpool);
  var_list_destroy(&vcfr->apool);

  genovcf_ptr_list_dealloc(&vcfr->vwait);
  genovcf_ptr_list_dealloc(&vcfr->vpool);
  genovar_ptr_list_dealloc(&vcfr->alist);
  genovar_ptr_list_dealloc(&vcfr->anchrom);
  genovar_ptr_list_dealloc(&vcfr->apool);
}

static inline void vcfr_get_chrom(bcf_hdr_t *vcfhdr, bcf1_t *ve, faidx_t *fai,
                                  int *refid, char **chrom, int *chromlen)
{
  if(*refid != ve->rid) {
    free(*chrom);
    *chrom = fai_fetch(fai, bcf_seqname(vcfhdr, ve), chromlen);
    if(*chrom == NULL) die("Cannot find chrom '%s'", bcf_seqname(vcfhdr, ve));
    *refid = ve->rid;
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
    if(genovar_ptr_list_len(&vr->apool) < n) var_list_populate(&vr->apool, n);

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

    // Load alleles alist, using insert sort
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
  genovar_ptr_list_append(&vr->apool, a);
  a->parent->nchildren--;
  // Re-add parent vcf to pool if no longer used
  if(a->parent->nchildren == 0) {
    printf("dropping var\n");
    genovcf_ptr_list_append(&vr->vwait, a->parent);
  }
  else printf("keeping var %zu\n", a->parent->nchildren);
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
  printf("shrinking vars %zu -> %zu\n", len, genovar_ptr_list_len(&vr->alist));
}

static int _genovcf_cmp(const void *aa, const void *bb)
{
  const GenoVCF *a = *(const GenoVCF*const*)aa, *b = *(const GenoVCF*const*)bb;
  return (long)a->vidx - (long)b->vidx;
}

static void vcfr_print_waiting(VcfReader *vr, htsFile *outfh, bcf_hdr_t *outhdr)
{
  // Sort waiting by vidx
  size_t len = genovcf_ptr_list_len(&vr->vwait);
  if(len == 0) return;

  printf("vwait: %zu\n", len);

  GenoVCF **vcfptr = genovcf_ptr_list_getptr(&vr->vwait, 0);
  qsort(vcfptr, len, sizeof(GenoVCF*), _genovcf_cmp);

  // Set cov to something
  // float *cov = calloc(2 * nsamples, sizeof(float));
  // for(i = 0; i < 2*nsamples; i++) cov[i] = i;
  size_t i, nsamples = bcf_hdr_nsamples(outhdr);
  float kcov[nsamples*10];
  int32_t nkmers[nsamples*10];

  for(i = 0; i < nsamples*10; i++) { nkmers[i] = i; kcov[i] = i*1.1; }

  // print
  size_t start = vr->nxtprint, end = vr->nxtprint + len;
  while(vr->nxtprint < end && vr->nxtprint == (*vcfptr)->vidx)
  {
    bcf1_t *v = &(*vcfptr)->v;

    bcf_update_format_int32(outhdr, v, vr->nkmer_tag, nkmers, nsamples * v->n_allele);
    bcf_update_format_float(outhdr, v, vr->cov_tag, kcov, nsamples * v->n_allele);

    if(bcf_write(outfh, outhdr, v) != 0) die("Cannot write record");
    vcfptr++;
    vr->nxtprint++;
  }
  // Take off those that were printed
  size_t nprinted = vr->nxtprint - start;
  printf("shift off %zu\n", nprinted);
  genovcf_ptr_list_shift(&vr->vwait, NULL, nprinted);

  if(genovcf_ptr_list_len(&vr->vwait) > 0) {
    vcfptr = genovcf_ptr_list_getptr(&vr->vwait, 0);
    printf("%zu vs %zu\n", vr->nxtprint, (*vcfptr)->vidx);
  }
}

static void genotype_vars(GenoVar **vars, size_t nvars,
                          size_t tgtidx, size_t ntgts,
                          const char *chrom, size_t chromlen,
                          Genotyper *gtyper, const dBGraph *db_graph)
{
  ctx_assert(ntgts <= nvars);

  GenoKmer *kmers = NULL;
  size_t i, end, nkmers, col;

  for(col = 0; col < db_graph->num_of_cols; col++)
  {
    nkmers = 0;
    if(nvars < max_gt_vars) {
      nkmers = genotyping_get_kmers(gtyper, (const GenoVar *const*)vars,
                                    nvars, tgtidx, ntgts,
                                    chrom, chromlen,
                                    db_graph->kmer_size, &kmers);
    }

    // Zero coverage
    for(i = tgtidx, end = tgtidx+ntgts; i < end; i++) {
      vars[i]->nkmers[0] = vars[i]->nkmers[1] = 0;
      vars[i]->sumcovg[0] = vars[i]->sumcovg[1] = 0;
    }

    for(i = 0; i < nkmers; i++) {
      bkey_get_covg(kmers[i].bkey, kmers[i].arbits,
                    vars+tgtidx, ntgts,
                    col, db_graph);
    }

    // TODO: update vars[]->parent coverage
    for(i = 0; i < nvars; i++) {
      
    }
  }
}

static inline int get_vcf_cov_end(GenoVar **vars, size_t nvars,
                                  size_t i, size_t endpos)
{
  while(i < nvars && endpos > vars[i]->pos) { i++; }
  return i;
}

// vars will be unordered after return
static void genotype_block(GenoVar **vars, size_t nvars,
                           const char *chrom, int chromlen,
                           Genotyper *gtyper, const dBGraph *db_graph)
{
  if(nvars < max_gt_vars)
  {
    genotype_vars(vars, nvars, 0, nvars, chrom, chromlen, gtyper, db_graph);
  }
  else
  {
    // do one at a time
    const int kmer_size = db_graph->kmer_size;
    int i, j;
    // genotype start/end, background start/end
    size_t gs = 0, ge, bs, be, tmp_be;

    while(gs < nvars)
    {
      for(i = j = (int)gs-1; i >= 0; i--) {
        if(genovar_end(vars[i]) + kmer_size > vars[gs]->pos) {
          SWAP(vars[i], vars[j]);
          j--;
        }
      }

      bs = j+1;
      ge = gs+1;
      be = get_vcf_cov_end(vars, nvars, ge, genovar_end(vars[ge-1]));

      // Try increasing ge, calculate be, accept if small enough
      for(i = ge+1; i < (int)nvars; i++) {
        tmp_be = get_vcf_cov_end(vars, nvars, i, genovar_end(vars[i-1]));
        if(tmp_be - bs < max_gt_vars) { ge = i; be = tmp_be; }
        else break;
      }

      ctx_assert2(bs<=gs && gs<=ge && ge<=be, "%zu %zu %zu %zu",bs,gs,ge,be);

      genotype_vars(vars+bs, be-bs, gs-bs, ge-gs, chrom, chromlen,
                    gtyper, db_graph);

      gs = ge;
    }
  }
}

static void genotype_vcf(htsFile *vcffh, bcf_hdr_t *vcfhdr,
                         htsFile *outfh, bcf_hdr_t *outhdr,
                         const char *path, faidx_t *fai,
                         const char *nkmer_tag, const char *cov_tag,
                         const dBGraph *db_graph)
{
  VcfReader vr;
  vcfr_alloc(&vr, path, nkmer_tag, cov_tag, vcffh, vcfhdr);

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
    vcfr_get_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                   &refid, &chrom, &chromlen);

    alist = genovar_ptr_list_getptr(&vr.alist, 0);
    alen = genovar_ptr_list_len(&vr.alist);
    ctx_assert(alen > 0);

    // Break into blocks
    endpos = genovar_end(alist[0]) + kmer_size;
    for(start = 0, end = 1; end < alen; end++) {
      if(endpos <= alist[end]->pos) {
        // end of block
        genotype_block(alist+start, end-start, chrom, chromlen, gtyper, db_graph);
        for(i = start; i < end; i++) vcfr_drop_var(&vr, i);
        start = end;
      }
      endpos = MAX2(endpos, genovar_end(alist[end])+kmer_size);
    }

    if(n == 0) {
      // End of chromosome -- do all variants
      genotype_block(alist+start, alen-start, chrom, chromlen, gtyper, db_graph);
      for(i = start; i < alen; i++) vcfr_drop_var(&vr, i);
    }

    // Shrink array if we dropped any vars
    vcfr_shrink_vars(&vr);
    vcfr_print_waiting(&vr, outfh, outhdr);
  }

  // Deal with remainder
  alist = genovar_ptr_list_getptr(&vr.alist, 0);
  alen = genovar_ptr_list_len(&vr.alist);
  if(alen) {
    // Get ref chromosome
    vcfr_get_chrom(vcfhdr, &mdc_list_get(&vr.alist, 0)->parent->v, fai,
                   &refid, &chrom, &chromlen);
    genotype_block(alist, alen, chrom, chromlen, gtyper, db_graph);
    for(i = 0; i < alen; i++) vcfr_drop_var(&vr, i);
    vcfr_shrink_vars(&vr);
  }
  vcfr_print_waiting(&vr, outfh, outhdr);

  free(chrom);
  genotyper_destroy(gtyper);
  vcfr_dealloc(&vr);
}

int ctx_geno(int argc, char **argv)
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
                              .must_exist_in_graph = true,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Add samples to vcf header
  bcf_hdr_t *outhdr = bcf_hdr_dup(vcfhdr);
  for(i = 0; i < db_graph.num_of_cols; i++) {
    status("Name %zu: %s", i, db_graph.ginfo[i].sample_name.b);
    bcf_hdr_add_sample(outhdr, db_graph.ginfo[i].sample_name.b);
  }

  char nkmer_tag[50], cov_tag[50];
  sprintf(nkmer_tag, "NK%zu", db_graph.kmer_size);
  sprintf(cov_tag,   "CK%zu", db_graph.kmer_size);

  // Add genotype format fields
  char descr[200];
  sprintf(descr,
          "##FORMAT=<ID=%s,Number=G,Type=Integer,"
          "Description=\"Number of kmers counted (k=%zu)\">\n",
          nkmer_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);

  sprintf(descr,
          "##FORMAT=<ID=%s,Number=G,Type=Float,Description=\"Mean kmer coverage (k=%zu)\">\n",
          cov_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, descr);

  if(bcf_hdr_write(outfh, outhdr) != 0)
    die("Cannot write header to: %s", futil_outpath_str(out_path));

  status("[vcfcov] Reading %s and adding coverage", vcf_path);
  genotype_vcf(vcffh, vcfhdr, outfh, outhdr, vcf_path, fai,
               nkmer_tag, cov_tag, &db_graph);

  status("  saved to: %s\n", out_path);

  bcf_hdr_destroy(vcfhdr);
  bcf_hdr_destroy(outhdr);
  hts_close(vcffh);
  hts_close(outfh);
  fai_destroy(fai);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
