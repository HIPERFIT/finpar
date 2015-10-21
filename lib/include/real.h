#ifndef FINPAR_REAL_H
#define FINPAR_REAL_H

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

#endif
