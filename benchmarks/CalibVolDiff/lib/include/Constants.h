#ifndef   CONSTANTS_H
#define   CONSTANTS_H

#define REAL3_CT 4
#define TRANSPOSE_UV 1

#ifdef REAL_IS_DOUBLE

/* Some OpenCL implementations require us to define a pragma if we
   want to use double-precision numbers. */
#ifdef __OPENCL_VERSION__
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

typedef double real_t;
#define REAL_FLAG "REAL_IS_DOUBLE"

#ifdef __OPENCL_VERSION__
#define real4_t double4
#define real3_t double3
#define real2_t double2
#endif

#elif REAL_IS_FLOAT

typedef float real_t;
#define REAL_FLAG "REAL_IS_FLOAT"

#ifdef __OPENCL_VERSION__
#define real4_t float4
#define real3_t float3
#define real2_t float2
#endif

#else

#error "Must set REAL_IS_DOUBLE or REAL_IS_FLOAT."

#endif

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
