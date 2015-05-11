#include "global.h"
#include "gpath_set.h"
#include "util.h"

// Save 16 bytes at the end of the sequence store
// This is relied on by GPathFollow
#define SEQ_STORE_PADDING 16

// If resize true, cannot do multithreaded but can resize array
// If resize false, die if out of mem, but can multithread
void gpath_set_alloc2(GPathSet *gpset, size_t ncols,
                      size_t initpaths, size_t initmem,
                      bool resize, bool keep_path_counts)
{
  GPathSet tmp = {.ncols = ncols, .can_resize = resize};

  // 1:1 split between seq and paths (6 bytes each)
  size_t entry_size = 0, klen_size = 0;

  entry_size = sizeof(GPath) + (ncols+7)/8;

  if(keep_path_counts)
    klen_size = sizeof(uint8_t)*ncols + sizeof(uint32_t);

  size_t mem_used = initpaths * (entry_size + klen_size);

  if(initmem < mem_used) {
    die("[GPathSet] Not enough memory for number of paths (%zu < %zu)",
        initmem, mem_used);
  }

  size_t seq_mem = initmem - mem_used;
  size_t col_mem = initpaths * ((ncols+7)/8);
  size_t seq_col_mem = col_mem + seq_mem;
  size_t total_mem = initpaths * (sizeof(GPath)+klen_size) + seq_col_mem;
  ctx_assert(total_mem <= initmem);

  char npathstr[50], colmemstr[50], seqmemstr[50], totalmemstr[50];
  ulong_to_str(initpaths, npathstr);
  bytes_to_str(col_mem, 1, colmemstr);
  bytes_to_str(seq_mem, 1, seqmemstr);
  bytes_to_str(total_mem, 1, totalmemstr);
  status("[GPathSet] Allocating for %s paths, %s colset, %s seq => %s total",
         npathstr, colmemstr, seqmemstr, totalmemstr);

  gpath_buf_alloc(&tmp.entries, initpaths);
  byte_buf_alloc(&tmp.seqs, seq_col_mem + SEQ_STORE_PADDING);

  if(keep_path_counts) {
    // madcrowlib allocates with calloc, so these are all zero'd
    byte_buf_alloc(&tmp.nseen_buf, initpaths * ncols);
    uint32_buf_alloc(&tmp.klen_buf, initpaths);
  } else {
    memset(&tmp.nseen_buf, 0, sizeof(tmp.nseen_buf));
    memset(&tmp.klen_buf, 0, sizeof(tmp.klen_buf));
  }

  memcpy(gpset, &tmp, sizeof(GPathSet));
}

// If resize true, cannot do multithreaded but can resize array
// If resize false, die if out of mem, but can multithread
void gpath_set_alloc(GPathSet *gpset, size_t ncols, size_t initmem,
                     bool resize, bool keep_path_counts)
{
  // Assume 8 bytes of sequence per path
  size_t entry_size
    = sizeof(GPath) + 8 + (ncols+7)/8 +
      (keep_path_counts ? sizeof(uint8_t)*ncols + sizeof(uint32_t) : 0);

  size_t nentries = initmem / entry_size;

  gpath_set_alloc2(gpset, ncols, nentries, initmem, resize, keep_path_counts);
}

void gpath_set_dealloc(GPathSet *gpset)
{
  gpath_buf_dealloc(&gpset->entries);
  byte_buf_dealloc(&gpset->seqs);
  byte_buf_dealloc(&gpset->nseen_buf);
  uint32_buf_dealloc(&gpset->klen_buf);
  memset(gpset, 0, sizeof(GPathSet));
}

void gpath_set_reset(GPathSet *gpset)
{
  gpath_buf_reset(&gpset->entries);
  byte_buf_reset(&gpset->seqs);
  byte_buf_reset(&gpset->nseen_buf);
}

void gpath_set_print_stats(const GPathSet *gpset)
{
  char paths_str[50], paths_cap_str[50], seq_str[50], seq_cap_str[50];
  ulong_to_str(gpset->entries.len, paths_str);
  ulong_to_str(gpset->entries.size, paths_cap_str);
  bytes_to_str(gpset->seqs.len, 1, seq_str);
  bytes_to_str(gpset->seqs.size, 1, seq_cap_str);
  status("[GPathSet] Paths: %s / %s [%.2f%%], seqs: %s / %s [%.2f%%] (%zu / %zu)",
         paths_str, paths_cap_str,
         (100.0 * gpset->entries.len) / gpset->entries.size,
         seq_str, seq_cap_str,
         (100.0 * gpset->seqs.len) / gpset->seqs.size,
         gpset->seqs.len, gpset->seqs.size);
}

// Get kmer length of a GPath
uint32_t gpath_set_get_klen(const GPathSet *gpset, const GPath *gpath)
{
  pkey_t pkey = gpset_get_pkey(gpset, gpath);
  return gpath_set_has_nseen(gpset) ? gpset->klen_buf.b[pkey] : 0;
}

uint8_t* gpath_set_get_nseen(const GPathSet *gpset, const GPath *gpath)
{
  pkey_t pkey = gpset_get_pkey(gpset, gpath);
  return gpath_set_has_nseen(gpset) ? &gpset->nseen_buf.b[pkey*gpset->ncols] : NULL;
}

// Copy nseen counts to dst from src
void gpath_set_nseen_sum_mt(const GPath *dst, GPathSet *dstset,
                            const GPath *src, const GPathSet *srcset)
{
  ctx_assert2(dstset->ncols == srcset->ncols, "ncols don't match");

  uint8_t *src_nseen = gpath_set_get_nseen(srcset, src);
  uint8_t *dst_nseen = gpath_set_get_nseen(dstset, dst);
  size_t i;

  if(src_nseen && dst_nseen) {
    for(i = 0; i < dstset->ncols; i++)
      safe_add_uint8(&dst_nseen[i], src_nseen[i]);
  }
}

// Resize buffers if needed and update internal pointers
void _check_resize(GPathSet *gpset, size_t req_num_bytes)
{
  const size_t ncols = gpset->ncols;
  size_t i;

  GPath *old_entries = gpset->entries.b;
  size_t old_num_entries = gpset->entries.size;
  gpath_buf_capacity(&gpset->entries, gpset->entries.len+1);

  if(old_entries != gpset->entries.b) {
    // Correct all path pointers
    for(i = 0; i < gpset->entries.len; i++) {
      if(gpset->entries.b[i].next != NULL) {
        gpset->entries.b[i].next = gpset->entries.b +
                                      (gpset->entries.b[i].next - old_entries);
      }
    }
  }

  if(old_num_entries != gpset->entries.size)
  {
    if(gpath_set_has_nseen(gpset)) {
      // Increase size of nseen buffer to match (zero'd by default)
      byte_buf_capacity(&gpset->nseen_buf, gpset->entries.size * ncols);
      uint32_buf_capacity(&gpset->klen_buf, gpset->entries.size);
    }
  }

  uint8_t *old_seq = gpset->seqs.b;
  byte_buf_capacity(&gpset->seqs, gpset->seqs.len+req_num_bytes+SEQ_STORE_PADDING);

  if(gpset->seqs.b != old_seq) {
    // Correct all seq pointers
    for(i = 0; i < gpset->entries.len; i++) {
      if(gpset->entries.b[i].seq != NULL) {
        gpset->entries.b[i].seq = gpset->seqs.b +
                                     (gpset->entries.b[i].seq - old_seq);
      }
    }
  }
}

// Always adds new path. If newpath could be a duplicate, use gpathhash
// Threadsafe only if resize is false. GPath* not safe to edit until it returns
// Copies newgpath.seq over and wipe new colset
GPath* gpath_set_add_mt(GPathSet *gpset, GPathNew newgpath)
{
  ctx_assert(newgpath.seq != NULL);

  pkey_t pkey;
  GPath *gpath;
  uint8_t *data;
  size_t colset_bytes = (gpset->ncols+7)/8;
  size_t junc_bytes = binary_seq_mem(newgpath.num_juncs);
  size_t nbytes = colset_bytes + junc_bytes;

  if(gpset->can_resize)
  {
    _check_resize(gpset, nbytes);
    pkey = gpath_buf_add(&gpset->entries, (GPath){.seq = NULL, .num_juncs = 0});
    gpath = &gpset->entries.b[pkey];
    data = gpset->seqs.b + gpset->seqs.len;
    gpset->seqs.len += nbytes;
  }
  else
  {
    gpath = gpset->entries.b + __sync_fetch_and_add((volatile size_t*)&gpset->entries.len, 1);
    data = gpset->seqs.b + __sync_fetch_and_add((volatile size_t*)&gpset->seqs.len, nbytes);

    if(gpath >= gpset->entries.b + gpset->entries.size ||
       data+nbytes+SEQ_STORE_PADDING > gpset->seqs.b + gpset->seqs.size)
    {
      gpath_set_print_stats(gpset);
      status("%zu >= %zu; %zu > %zu nbytes: %zu\n",
             gpath - gpset->entries.b,
             gpset->entries.size,
             data+nbytes+SEQ_STORE_PADDING - gpset->seqs.b,
             gpset->seqs.size,
             nbytes);
      die("Out of memory");
    }

    pkey = gpath - gpset->entries.b;
  }

  uint8_t *colset = data;
  gpath->seq = data + colset_bytes;
  gpath->num_juncs = newgpath.num_juncs;
  gpath->orient = newgpath.orient;
  gpath->next = NULL;

  // copy seq and zero colset
  memcpy(gpath->seq, newgpath.seq, junc_bytes);

  if(newgpath.colset)
    memcpy(colset, newgpath.colset, colset_bytes);
  else
    memset(colset, 0, colset_bytes);

  // klen, nseen
  if(gpath_set_has_nseen(gpset))
  {
    gpset->klen_buf.b[pkey] = newgpath.klen;
    __sync_fetch_and_add((volatile size_t*)&gpset->klen_buf.len, 1);
    __sync_fetch_and_add((volatile size_t*)&gpset->nseen_buf.len, gpset->ncols);

    uint8_t *nseen = gpath_set_get_nseen(gpset, gpath);
    ctx_assert(nseen != NULL);

    if(newgpath.nseen)
      memcpy(nseen, newgpath.nseen, gpset->ncols);
    else
      memset(nseen, 0, gpset->ncols);
  }

  return gpath;
}

GPathNew gpath_set_get(const GPathSet *gpset, const GPath *gpath)
{
  GPathNew newgpath = {.seq = gpath->seq,
                       .colset = gpath_get_colset(gpath, gpset->ncols),
                       .nseen = gpath_set_get_nseen(gpset, gpath),
                       .klen = gpath_set_get_klen(gpset, gpath),
                       .num_juncs = gpath->num_juncs,
                       .orient = gpath->orient};
  return newgpath;
}
