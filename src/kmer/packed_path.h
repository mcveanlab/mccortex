#ifndef PACKED_PATH_H_
#define PACKED_PATH_H_

typedef uint64_t PathIndex;
typedef uint16_t PathLen;

#define PATH_NULL UINT64_MAX
#define PATH_LEN_BITS 15
#define MAX_PATHLEN ((1UL<<PATH_LEN_BITS)-1)

//
// Functions for PackedPaths
//
// PackedPath doesn't have a struct, instead is laid out as:
//   {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}
// Where N=round_up(num_of_colours/8), M=round_up(len/4)
// len holds both the length (lower 15 bits), and orientations (top bit)

// Number of bytes a packed path takes up
#define packedpath_mem2(colbytes,nbytes) \
        (sizeof(PathIndex)+(colbytes)+sizeof(PathLen)+(nbytes))

#define packedpath_mem(p,colbytes) \
        packedpath_mem2(colbytes,packedpath_pbytes(p,colbytes))

#define packedpath_pbytes(p,colbytes) \
        packedpath_len_nbytes(packedpath_get_len(p,colbytes))

//
// Previous
//

static inline PathIndex packedpath_get_prev(const uint8_t *ptr) {
  PathIndex idx; memcpy(&idx, ptr, sizeof(PathIndex)); return idx;
}

static inline void packedpath_set_prev(uint8_t *ptr, PathIndex idx) {
  memcpy(ptr, &idx, sizeof(PathIndex));
}

//
// Colour bitset
//
#define packedpath_del_col(p,idx,col) bitset_del((p)+sizeof(PathIndex),col)
#define packedpath_set_col(p,idx,col) bitset_set((p)+sizeof(PathIndex),col)
#define packedpath_has_col(p,idx,col) bitset_get((p)+sizeof(PathIndex),col)

// Get pointer to colset
#define packedpath_get_colset(ptr) ((ptr) + sizeof(PathIndex))

//
// Length and Orientation
//

// Number of bytes needed to store n bases (2 bits per base, 4 per byte)
#define packedpath_len_nbytes(nbases) (((nbases)+3)/4)

#define PP_ORIENTMASK (1U << PATH_LEN_BITS)
#define PP_LENMASK (PP_ORIENTMASK-1)

#define packedpath_len(w) ((w) & ~PP_ORIENTMASK)
#define packedpath_or(w) ((w) >> PATH_LEN_BITS)
#define packedpath_combine_lenorient(len,orient) ((len)|((PathLen)(orient)<<PATH_LEN_BITS))

static inline PathLen packedpath_get_lenword(const uint8_t *ptr, size_t colbytes)
{
  PathLen len; memcpy(&len, ptr+sizeof(PathIndex)+colbytes, sizeof(PathLen));
  return len;
}

static inline void packedpath_set_lenword(uint8_t *ptr, size_t colbytes,
                                          PathLen len_orient)
{
  memcpy(ptr+sizeof(PathIndex)+colbytes, &len_orient, sizeof(PathLen));
}

static inline PathLen packedpath_get_len(const uint8_t *ptr, size_t colbytes) {
  return packedpath_len(packedpath_get_lenword(ptr,colbytes));
}

static inline Orientation packedpath_get_orient(const uint8_t *ptr, size_t colbytes) {
  return packedpath_or(packedpath_get_lenword(ptr,colbytes));
}

static inline void packedpath_set_len(uint8_t *ptr, size_t colbytes, PathLen len)
{
  // Combine with orientation
  len |= packedpath_get_lenword(ptr, colbytes) & PP_ORIENTMASK;
  packedpath_set_lenword(ptr, colbytes, len);
}

static inline void packedpath_set_orient(uint8_t *ptr, size_t colbytes,
                                         Orientation orient)
{
  // Combine with orientation
  PathLen len, len_orient;
  len = packedpath_get_lenword(ptr, colbytes) & PP_LENMASK;
  len_orient = packedpath_combine_lenorient(len, orient);
  packedpath_set_lenword(ptr, colbytes, len_orient);
}

// Get length and orientation together
static inline PathLen packedpack_get_len_orient(const uint8_t *ptr,
                                                size_t colbytes, PathLen *len,
                                                Orientation *orient)
{
  PathLen lenword = packedpath_get_lenword(ptr, colbytes);
  *len = packedpath_len(lenword);
  *orient = packedpath_or(lenword);
  return lenword;
}

static inline void packedpack_set_len_orient(uint8_t *ptr, size_t colbytes,
                                             PathLen len, Orientation orient)
{
  PathLen len_orient = packedpath_combine_lenorient(len, orient);
  packedpath_set_lenword(ptr, colbytes, len_orient);
}

//
// Path
//

// Get pointer to first byte of path
#define packedpath_path(p,colbytes) ((p)+sizeof(PathIndex)+colbytes+sizeof(PathLen))

#endif /* PACKED_PATH_H_ */
