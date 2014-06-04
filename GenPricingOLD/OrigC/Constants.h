#ifndef   CONSTANTS_H
#define   CONSTANTS_H

/****************************************************/
/**
 TO DO:  1. solve include names from kernels
         2. set TILE automatically somehow
         3. set MY_UCHAR to UCHAR (preferred if supported by GPU) or UINT
         4. in ``discriminate_cost_model'' decide what
            (a) if you want to set an upper bound for CHUNK, i.e.,
                if(ro_scal->logCHUNK > 6) ro_scal->logCHUNK = 6;
            (b) the ratio of global memory you want to use: 3/4 or 2/5 or ... ?
**/
/****************************************************/

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

#define dev_id 0 // the GPU device number to run the code on

typedef struct {
	/*** NUM ITERS ***/
	UINT   mc_iter_num ;
	UINT   mc_bigit_num;
	UINT   sobol_count_ini;

	/* SOBOL RO data */
	UINT   sobol_bit_count;
	UINT   sobol_dim; // == md_dim * md_nb_path_dates

	/* MD data */
	UINT   md_dim;
	UINT   md_nb_path_dates;
	UINT   rng_dim;
	UINT   rng_nb_dates;

	/* MD instances data */
	UINT   inst_num;         // md->nb_instances == num_models, where md :: driver*
	UINT   num_contracts;    // pc->nb_contracts,               where pc :: mlfi_nmc_pricing_code*
	UINT   num_det_pricers;  // pc->nb_deterministic_pricers,   where pc :: mlfi_nmc_pricing_code*
	UINT   num_cash_flows;   // pc->nb_cash_flows,              where pc :: mlfi_nmc_pricing_code*

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

	UINT logCHUNK;
	UINT CHUNK;

	UINT TILE_FACT;
	UINT TOT_MC_ITER_NUM;

	/** SOBOL RO last_denominator_inv **/
	REAL   sobol_last_den_inv;

	/** Model attributes -- in the general case these should be arrays**/
	//REAL   model_deter_val0;  // promoted to array
	//REAL   model_discounts0;  // promoted to array

} LoopROScalars __attribute__ ((aligned (16)));

#endif