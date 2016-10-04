#ifndef CORTEX_TYPES_H_
#define CORTEX_TYPES_H_

typedef size_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#define COVG_MAX UINT_MAX

#define SAFE_ADD_COVG(a,b) ((uint64_t)(a)+(b) > COVG_MAX ? COVG_MAX : (a)+(b))
#define SAFE_SUM_COVG(a,b) ((a) = SAFE_ADD_COVG((a), (b)))

typedef uint8_t Orientation;
#define FORWARD 0
#define REVERSE 1

typedef uint8_t ReadMateDir;
#define READPAIR_FF 0
#define READPAIR_FR 1
#define READPAIR_RF 2
#define READPAIR_RR 3
// See seq_reader.h to get string representations (MP_DIR_STRS[dir])

#define read_mate_r1(r) ((r)&2)
#define read_mate_r2(r) ((r)&1)

// don't ever use the top bit of hkey, used later for orientation
typedef uint64_t hkey_t;

typedef struct {
  hkey_t orient:1, key:63;
} dBNode;

#endif /* CORTEX_TYPES_H_ */
