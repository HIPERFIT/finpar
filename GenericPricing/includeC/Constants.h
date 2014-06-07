#ifndef   CONSTANTS_H
#define   CONSTANTS_H

#include "Optimizations.h"

#define WARP            (1<<lgWARP) 
#define logMAX_CHUNK    8

#if (_OPTIMIZATION_USE_FLOATS)
    typedef float        REAL;
    typedef unsigned int ULONG;
#else
    #pragma OPENCL EXTENSION cl_khr_fp64: enable
    typedef double         REAL;
    typedef unsigned long  ULONG;
#endif

typedef unsigned int   UINT;
typedef unsigned char  UCHAR;

#define MY_UCHAR UCHAR

typedef struct {
    /*** Contract Characteristics ***/
    UINT   contract;        // identifies the dataset
    UINT   num_mcits;       // # of Monte-Carlo iterations
    UINT   sobol_bits;      // # of bits in Sobol-number rep
    UINT   num_under;       // # of underlyings, i.e., paths
    UINT   num_dates;       // # of dates on a path

    /* MD instances data */
    UINT   num_models;      // # of models (each computes a price)
    UINT   num_det_pricers; // # of deterministic pricers (used by payoff)
    UINT   num_cash_flows;  // # of cash flows, i.e., discounts

    /*** GPU execution info ***/
    UINT   num_gpuits;      // # of iterations that can execute on GPU at a time
    UINT   num_mcbigits;    // # of chunked-unrolled iterations

    UINT   sobol_count_ini; // initial number from which the Sobol sequence starts 
    REAL   sob_norm_fact;   // maximal sobol int 1/(1<<sobol_bits)

    UINT   chunk;           // loop strip mining factor    
    UINT   log_chunk;       // log_2 of strip mining factor

    UINT   logBLOCK;
    UINT   BLOCK;

//    UINT   sobol_dim;       // == md_dim * md_nb_path_dates
} LoopROScalars __attribute__ ((aligned (32)));

#endif
