#ifndef CORTEX_TYPES_H_
#define CORTEX_TYPES_H_

typedef size_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#define COVG_MAX UINT_MAX

#define FORWARD 0
#define REVERSE 1
typedef uint_fast8_t Orientation;

#define true 1
#define false 0
typedef uint_fast8_t boolean;

typedef uint_fast8_t ReadMateDir;
#define READPAIR_FF 0
#define READPAIR_FR 1
#define READPAIR_RF 2
#define READPAIR_RR 3

#define read_mate_r1(r) ((r)&2)
#define read_mate_r2(r) ((r)&1)

// DEV: define combine READPAIR

#endif /* CORTEX_TYPES_H_ */
