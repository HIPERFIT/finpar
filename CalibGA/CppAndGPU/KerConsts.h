#ifndef KERNEL_CONSTANTS
#define KERNEL_CONSTANTS

/**
 * If set than Sobol's pseudo-rand num gen is used,
 * Otherwise (if GPU_VERSION!=2) `rand' from C is used.
 */
#define WITH_SOBOL 1 

/** 
 * 0 -> Multi-Core Execution
 * 1 -> Only the pricing is done on GPU,
 *        mutation, crossover done on CPU
 * 2 -> Everything on GPU using Sobol's 
 *        independent formula for random numbers.
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
    #define EPS0   (0.3e-2) //was (1.0e-2)
    #define EPS    (0.2e-4) // -4 //was (1.0e-5)
    #define PI     (3.14159265358)
    //const REAL PI    = 3.14159265358; // define PI as a const, changes slightly results.
//    const REAL INFTY = 1.0e38;

#else
    #pragma OPENCL EXTENSION cl_khr_fp64: enable

    typedef double REAL;
    #define EPS0   (1.0e-3)
//    #define EPS0   (1.0e-6)
    #define EPS    (1.0e-5)
    #define PI     (3.1415926535897932384626433832795)
    //const REAL INFTY = 1.0e49;   //// don't make INFTY a macro -- wrong results!
#endif

#define WARP   32
#define lgWARP 5

#define LWG_EG 256

// ugaussian_Pinv(KKK)=1.0e~4
#define KKK -3.71901648545568   

//const REAL R     = 0.03;
#define R     0.03
// TODAY should be passed as a compile-time parameter,
// or as a constant structure!

#endif // KERNEL_CONSTANTS

