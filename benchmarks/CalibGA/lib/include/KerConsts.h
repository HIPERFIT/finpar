#ifndef KERNEL_CONSTANTS
#define KERNEL_CONSTANTS

/**
 * If set then Sobol's pseudo-rand num gen is used,
 * Otherwise (if GPU_VERSION!=2) `rand' from C is used.
 */
#define WITH_SOBOL 1 

/** 
 * 0 -> Multi-Core Execution
 * 1 -> Only the pricing is done on GPU,
 *        mutation, crossover done on CPU
 * 2 -> Everything on GPU using Sobol's 
 *        independent formula for pseudo-
 *        random numbers.
 */
#define GPU_VERSION 2

/**
 * If WITH_FLOAT == 0 THEN double precision is used, 
 * Otherwise single float precision.
 * 
 * If you modify this flag, DO NOT forget to perform
 *    a dummy change in file `SrcCL/CalibKers.cl',
 *    so that the cached version will be invalidated!
 */
#define WITH_FLOAT 1

/*********************************************/
/*********************************************/
/*********************************************/

#if WITH_FLOAT
    typedef float  REAL;
    #define EPS0   (0.3e-2)
    #define EPS    (0.2e-4) 
    #define PI     (3.14159265358)

#else
    #pragma OPENCL EXTENSION cl_khr_fp64: enable

    typedef double REAL;
    #define EPS0   (1.0e-3)
    #define EPS    (1.0e-5)
    #define PI     (3.1415926535897932384626433832795)

#endif

#define WARP   (1<<lgWARP)

#define LWG_EG 256

// ugaussian_Pinv(KKK)=1.0e~4
#define KKK -3.71901648545568   

#define R     0.03

#endif // KERNEL_CONSTANTS

