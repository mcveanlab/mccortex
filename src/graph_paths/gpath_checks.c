#include "global.h"
#include "gpath_checks.h"
#include "db_node.h"
#include "binary_seq.h"

typedef struct {
  size_t col, from_col;
  const char *name, *src;
} ColourSample;

static int _col_sample_cmp(const void *aa, const void *bb)
{
  const ColourSample *a = (const ColourSample *)aa;
  const ColourSample *b = (const ColourSample *)bb;
  if(a->col != b->col) return (a->col < b->col ? -1 : 1);
  return strcmp(a->name, b->name);
}

// Similar to path_file_reader.c:path_file_load_check()
// Check kmer size matches and sample names match
void graphs_gpaths_compatible(const GraphFileReader *graphs, size_t num_graphs,
                              const GPathReader *gpaths, size_t num_gpaths)
{
  size_t g, p, kmer_size, kmer_size2;
  size_t ctx_max_cols = 0, ctp_max_cols = 0;
  size_t ctx_max_kmers = 0, ctp_max_kmers = 0;
  size_t colours_loaded = 0;

  if(num_graphs) kmer_size = graphs[0].hdr.kmer_size;
  else if(num_gpaths) kmer_size = gpath_reader_get_kmer_size(&gpaths[0]);
  else return; // no graph or path files

  for(g = 0; g < num_graphs; g++) {
    kmer_size2 = graphs[g].hdr.kmer_size;
    if(kmer_size2 != kmer_size) {
      die("Kmer-size doesn't match between files [%zu vs %zu]: %s",
          kmer_size, kmer_size2, graphs[g].fltr.orig_path.buff);
    }
    ctx_max_cols = MAX2(ctx_max_cols, file_filter_usedcols(&graphs[g].fltr));
    ctx_max_kmers = MAX2(ctx_max_kmers, graphs[g].num_of_kmers);
    colours_loaded += graphs[g].fltr.ncols;
  }

  for(p = 0; p < num_gpaths; p++) {
    kmer_size2 = gpath_reader_get_kmer_size(&gpaths[p]);
    if(kmer_size2 != kmer_size) {
      die("Kmer-size doesn't match between files [%zu vs %zu]: %s",
          kmer_size, kmer_size2, gpaths[p].fltr.orig_path.buff);
    }
    ctp_max_cols = MAX2(ctp_max_cols, file_filter_usedcols(&gpaths[p].fltr));
    ctp_max_kmers = MAX2(ctp_max_kmers, gpath_reader_get_num_kmers(&gpaths[p]));
    colours_loaded += gpaths[g].fltr.ncols;
  }

  const FileFilter *fltr = (num_graphs ? &graphs[0].fltr : &gpaths[0].fltr);
  db_graph_check_kmer_size(kmer_size, fltr->orig_path.buff);

  // ctx_max_kmers may be zero if reading from a stream
  if(ctx_max_kmers > 0 && ctp_max_kmers > ctx_max_kmers) {
    die("More kmers in path files than in graph files! (%zu > %zu)",
        ctp_max_kmers, ctx_max_kmers);
  }

  if(num_graphs > 0 && ctp_max_cols > ctx_max_cols) {
    die("More colours in path files than in graph files! (%zu > %zu)",
        ctp_max_cols, ctx_max_cols);
  }

  // Check sample names
  ColourSample *samples = ctx_calloc(colours_loaded, sizeof(ColourSample));

  size_t pinto, pfrom, ginto, gfrom, i;
  const char *pname, *gname;
  colours_loaded = 0;

  for(p = 0; p < num_gpaths; p++) {
    for(i = 0; i < gpaths[p].fltr.ncols; i++) {
      pinto = file_filter_intocol(&gpaths[p].fltr, i);
      pfrom = file_filter_fromcol(&gpaths[p].fltr, i);
      pname = gpath_reader_get_sample_name(&gpaths[p], pfrom);
      samples[colours_loaded++]
        = (ColourSample){.col = pinto, .name = pname,
                         .src = gpaths[p].fltr.orig_path.buff,
                         .from_col = gpaths[p].fltr.cols[i]};
    }
  }

  for(g = 0; g < num_graphs; g++) {
    for(i = 0; i < graphs[g].fltr.ncols; i++) {
      ginto = file_filter_intocol(&graphs[g].fltr, i);
      gfrom = file_filter_fromcol(&graphs[g].fltr, i);
      gname = graphs[g].hdr.ginfo[gfrom].sample_name.buff;
      samples[colours_loaded++]
        = (ColourSample){.col = ginto, .name = gname,
                         .src = graphs[g].fltr.orig_path.buff,
                         .from_col = graphs[g].fltr.cols[i]};
    }
  }

  // Sort by colour number than sample name
  qsort(samples, colours_loaded, sizeof(ColourSample), _col_sample_cmp);

  for(i = 0; i+1 < colours_loaded; i++) {
    if(samples[i].col == samples[i+1].col &&
       strcmp(samples[i].name, samples[i+1].name) != 0)
    {
      die("Sample names don't match\n%s:%zu%s\n%s:%zu%s\n",
          samples[i].src, samples[i].from_col, samples[i].name,
          samples[i+1].src, samples[i+1].from_col, samples[i+1].name);
    }
  }

  ctx_free(samples);
}

//
// Integrity checks on graph+paths
//

// 1) check dBNode following `node` has indegree >1 in sample ctxcol
// 2) follow path, check each junction matches up with a node with outdegree >1
// col is graph colour
bool gpath_checks_path_col(dBNode node, const GPath *gpath, int exp_klen,
                           size_t ctxcol, const dBGraph *db_graph)
{
  ctx_assert_ret(db_graph->num_edge_cols == db_graph->num_of_cols ||
                 db_graph->node_in_cols != NULL);

  BinaryKmer bkmer;
  Edges edges;
  dBNode nodes[4];
  Nucleotide nucs[4];
  size_t i, j, n, edgecol = db_graph->num_edge_cols > 1 ? ctxcol : 0;
  // length is kmers and junctions
  size_t klen, plen;

  // fprintf(stderr, " == graph_paths_check_valid()\n");

  for(klen = 0, plen = 0; plen < gpath->num_juncs; klen++)
  {
    bkmer = db_node_get_bkmer(db_graph, node.key);
    edges = db_node_get_edges(db_graph, node.key, edgecol);

    // Check this node is in this colour
    if(db_graph->node_in_cols != NULL) {
      ctx_assert_ret(db_node_has_col(db_graph, node.key, ctxcol));
    } else if(db_graph->col_covgs != NULL) {
      ctx_assert_ret(db_node_get_covg(db_graph, node.key, ctxcol) > 0);
    }

    #ifdef CTXVERBOSE
      char bkmerstr[MAX_KMER_SIZE+1];
      binary_kmer_to_str(bkmer, db_graph->kmer_size, bkmerstr);
      fprintf(stderr, "klen: %zu plen: %zu %zu:%i %s\n",
             klen, plen, (size_t)node.key, node.orient, bkmerstr);
    #endif

    if(klen == 1) {
      // Check nodes[1] has indegree > 1
      dBNode rnode = db_node_reverse(node);
      Edges backedges = db_node_edges_in_col(rnode, ctxcol, db_graph);
      int outdegree = edges_get_outdegree(backedges, rnode.orient);
      if(outdegree <= 1) {
        char bkstr[MAX_KMER_SIZE+1];
        binary_kmer_to_str(db_node_get_bkmer(db_graph, node.key), db_graph->kmer_size, bkstr);
        status("outdegree: %i col: %zu kmer: %s", (int)outdegree, ctxcol, bkstr);
      }
      ctx_assert_ret(outdegree > 1);
    }

    n = db_graph_next_nodes(db_graph, bkmer, node.orient,
                            edges, nodes, nucs);

    ctx_assert_ret(n > 0);

    // Reduce to nodes in our colour if edges limited
    if(db_graph->num_edge_cols == 1 && db_graph->node_in_cols != NULL) {
      for(i = 0, j = 0; i < n; i++) {
        if(db_node_has_col(db_graph, nodes[i].key, ctxcol)) {
          nodes[j] = nodes[i];
          nucs[j] = nucs[i];
          j++;
        }
      }
      n = j; // update number of next nodes
      ctx_assert_ret(n > 0);
    }

    // If fork check nucleotide
    if(n > 1) {
      Nucleotide expbase = binary_seq_get(gpath->seq, plen);

      for(i = 0; i < n && nucs[i] != expbase; i++);
      if(i == n) {
        fprintf(stderr, "plen: %zu expected: %c\n", plen, dna_nuc_to_char(expbase));
        fprintf(stderr, "Got: ");
        for(i = 0; i < n; i++) fprintf(stderr, " %c", dna_nuc_to_char(nucs[i]));
        fprintf(stderr, "\n");
      }
      ctx_assert_ret(i < n && nucs[i] == expbase);
      node = nodes[i];
      plen++;
    }
    else {
      node = nodes[0];
    }
  }

  // We exit before the last kmer, need to add it
  klen++;

  // Check kmer length is correct
  if(exp_klen != -1) {
    ctx_assert2(klen == (size_t)exp_klen, "act: %zu vs %zu", klen, (size_t)exp_klen);
  }

  return true;
}

// Returns false on first error
bool gpath_checks_path(hkey_t hkey, const GPath *gpath, int exp_klen,
                       const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  size_t i, ncols = gpstore->gpset.ncols;

  dBNode node = {.key = hkey, .orient = gpath->orient};

  // Check at least one colour is set
  uint8_t *colset = gpath_get_colset(gpath, ncols), cumm = 0;
  for(i = 0; i < ncols; i++) cumm |= colset[i];
  ctx_assert(cumm != 0);

  // Check for each colour the path has
  for(i = 0; i < ncols; i++) {
    if(bitset_get(colset, i)) {
      ctx_assert(gpath_checks_path_col(node, gpath, exp_klen, i, db_graph));
    }
  }

  return true;
}

static bool _kmer_check_paths(hkey_t hkey, const dBGraph *db_graph,
                              size_t *npaths_ptr, size_t *nkmers_ptr)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  size_t num_gpaths = 0;
  GPath *gpath;
  int exp_klen = -1;

  for(gpath = gpstore->paths_all[hkey]; gpath != NULL; gpath = gpath->next)
  {
    if(gpath_set_has_nseen(gpset)) exp_klen = gpath_set_get_klen(gpset, gpath);
    ctx_assert_ret(gpath_checks_path(hkey, gpath, exp_klen, db_graph));
    num_gpaths++;
  }

  *npaths_ptr += num_gpaths;
  *nkmers_ptr += (num_gpaths > 0);
  return true;
}

// Returns false on first error
bool gpath_checks_all_paths(const dBGraph *db_graph)
{
  status("[GPathCheck] Running paths check...");

  size_t num_gpaths = 0, num_kmers = 0, act_num_gpaths, act_num_kmers;

  HASH_ITERATE(&db_graph->ht, _kmer_check_paths, db_graph,
               &num_gpaths, &num_kmers);

  act_num_gpaths = db_graph->gpstore.gpset.entries.len;
  act_num_kmers = db_graph->gpstore.num_kmers_with_paths;

  ctx_assert_ret2(num_gpaths == act_num_gpaths, "%zu vs %zu", num_gpaths, act_num_gpaths);
  ctx_assert_ret2(num_kmers == act_num_kmers, "%zu vs %zu", num_kmers, act_num_kmers);
  return true;
}

// For debugging
static void _gpstore_update_counts(hkey_t hkey, const dBGraph *db_graph,
                                   size_t *nvisited_ptr, size_t *nkmers_ptr,
                                   size_t *npaths_ptr)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  GPath *gpath = gpath_store_fetch(gpstore, hkey);
  size_t npaths = 0;

  (*nvisited_ptr)++;

  if(gpath == NULL) return;
  (*nkmers_ptr)++;

  // Count paths and coloured paths
  for(npaths = 0; gpath != NULL; gpath = gpath->next, npaths++) {}

  (*npaths_ptr) += npaths;
}

void gpath_checks_counts(const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  size_t nvisited = 0, nkmers = 0, npaths = 0;

  HASH_ITERATE(&db_graph->ht, _gpstore_update_counts,
               db_graph, &nvisited, &nkmers, &npaths);

  ctx_assert2(nvisited == db_graph->ht.num_kmers, "%zu vs %zu",
              nvisited, (size_t)db_graph->ht.num_kmers);
  ctx_assert2(nkmers == gpstore->num_kmers_with_paths, "%zu vs %zu",
              nkmers, (size_t)gpstore->num_kmers_with_paths);
  ctx_assert2(npaths == gpstore->gpset.entries.len, "%zu vs %zu",
              npaths, (size_t)gpstore->gpset.entries.len);
}
