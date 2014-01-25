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

typedef struct dBGraph dBGraph;

#endif /* CORTEX_TYPES_H_ */
