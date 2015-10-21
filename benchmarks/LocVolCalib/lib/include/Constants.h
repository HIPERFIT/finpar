#ifndef   CONSTANTS_H
#define   CONSTANTS_H

#include "real.h"

#define REAL3_CT 4
#define TRANSPOSE_UV 1

#define WARP            (1<<lgWARP)

#define BLOCK_DIM         16
#define logWORKGROUP_SIZE 8
#define WORKGROUP_SIZE    (1<<logWORKGROUP_SIZE)

typedef struct {
    real_t   nu;
    real_t   alpha;
    real_t   beta;

    real_t   dtInv;
    real_t   timeline_i;
    int      t_ind;

    unsigned int NUM_X;
    unsigned int NUM_Y;
    unsigned int NUM_XY;
} RWScalars __attribute__ ((aligned (32)));


#endif // end constants file
