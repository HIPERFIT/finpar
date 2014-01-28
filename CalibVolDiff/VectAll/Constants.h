#ifndef   CONSTANTS_H
#define   CONSTANTS_H

#define WITH_FLOATS 1

#define REAL3_CT 4
#define TRANSPOSE_UV 1

#if (WITH_FLOATS)
    typedef float        REAL;
    typedef unsigned int ULONG;
    #define TILE 16
#else
    #pragma OPENCL EXTENSION cl_khr_fp64: enable
    typedef double         REAL;
    typedef unsigned long  ULONG;
    #define TILE 8
#endif

#define BLOCK_DIM 16
#define logWORKGROUP_SIZE 8
#define    WORKGROUP_SIZE 256 //256   // Put Back to 256 please!

typedef struct {
    REAL   nu;
    REAL   alpha;
    REAL   beta;

    REAL   dtInv;
    REAL   timeline_i;
    int    t_ind;

    unsigned int NUM_X;
    unsigned int NUM_Y;
    unsigned int NUM_XY;

    //REAL   timeline_i;
    //REAL   timeline_ip1;
} RWScalars __attribute__ ((aligned (32)));


#endif // end constants file
