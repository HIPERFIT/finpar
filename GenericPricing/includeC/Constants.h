#ifndef   CONSTANTS_H
#define   CONSTANTS_H

#include "Optimizations.h"

#define logWARP         5  //ro_scal->logBLOCK // 5
#define WARP            32 //ro_scal->BLOCK    // 32

#if (_OPTIMIZATION_USE_FLOATS)
	typedef float        REAL;
	typedef unsigned int ULONG;
	#define TILE _OPTIMIZATION_TILE
#else
	#pragma OPENCL EXTENSION cl_khr_fp64: enable
	typedef double         REAL;
	typedef unsigned long  ULONG;
	#define TILE _OPTIMIZATION_TILE/2
#endif

typedef unsigned int   UINT;
typedef unsigned char  UCHAR;

#define MY_UCHAR UCHAR
//#define MY_UCHAR UINT

#define DEV_ID 0 // the GPU device number to run the code on

typedef struct {
    // Is the function-member of a structure supported in OpenCL?
    void init() { 
        chunk = 64; // must be a power of two!
        sob_norm_fact = 1.0 / (1<<sobol_bits); 
    }

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

	UINT   chunk;           // loop strip mining factor    
    UINT   log_chunk;       // log_2 of strip mining factor

//    UINT   sobol_dim;       // == md_dim * md_nb_path_dates

	/* Brownian Bridge */
	UINT   bb_l;

	/*****************************/
	/* Indexes in RO structures! */
	/*****************************/
	// 1. Brownian Bridge
	UINT bb_li_beg;
	UINT bb_bi_beg;
	UINT bb_ri_beg;
	UINT bb_rw_beg;
	UINT bb_lw_beg;
	UINT bb_sd_beg;

	//2. MD Instances
	UINT inst_c_beg;
	UINT inst_vols_beg;
	UINT inst_drifts_beg;
	UINT inst_trajRO_beg;

	// 3. Model Data Attributes:
	UINT pc_deter_val_beg;
	UINT pc_discounts_beg;

	// 4. BLOCK/CHUNK
	UINT logBLOCK;
	UINT BLOCK;

	

	UINT TILE_FACT;
	UINT TOT_MC_ITER_NUM;

	/** SOBOL RO last_denominator_inv **/
	//REAL   sobol_last_den_inv;
    REAL sob_norm_fact;

	/** Model attributes -- in the general case these should be arrays**/
	//REAL   model_deter_val0;  // promoted to array
	//REAL   model_discounts0;  // promoted to array

} LoopROScalars __attribute__ ((aligned (16)));

#endif
