/************************************/
/************* MACROS ***************/
/************************************/
#include "Constants.h"

#define underlyings(i,j)   \
		(inst_traj[i*ro_scal->md_dim + j]) 

/********************************************/
/********** SOBOL NUMER GENERATOR ***********/
/********************************************/


void mlfi_genmatrix_uniformGPUind (
		UINT                  			seq_count,
		__constant  LoopROScalars* 		ro_scal,
		__constant  int* 				sobol_v_dir,
		            int* 				sobol_last_num_vec,
		            REAL* 				md_zd
) {
	UINT  j, k, gs, gv_k = 0;

	seq_count += 1;
	gs = seq_count >> 1;
	gs = seq_count ^  gs;

	for( j = 0; j < ro_scal->sobol_dim; j++ )
		sobol_last_num_vec[j] = 0;

	for( k = 0; k < ro_scal->sobol_bit_count; ++k ) {
		if(gs & 1) {
            __constant int* dir_vect
				= sobol_v_dir + k*ro_scal->sobol_dim;
			for( j=0; j < ro_scal->sobol_dim; j++ ) {
				// xor term g_k * v_k to direction i
				sobol_last_num_vec[j] ^= dir_vect[j];
			}
		}
        gs = gs >> 1;
	}
	for( j = 0; j < ro_scal->sobol_dim; j++ ) {
		md_zd[j] = sobol_last_num_vec[j] * ro_scal->sobol_last_den_inv;
	}
}

inline void mlfi_genmatrix_uniformGPUrecOpt(
		UINT                  			f_ind,
		__constant  LoopROScalars* 		ro_scal,
		__constant  int* 				sobol_v_dir,
		            int* 				sobol_last_num_vec,
		            REAL* 				md_zd
) {
	UINT j;

	for(j=0; j < ro_scal->sobol_dim; j++) {
	    sobol_last_num_vec[j] ^= sobol_v_dir[ f_ind*ro_scal->sobol_dim + j];
		md_zd[j]               = sobol_last_num_vec[j] * ro_scal->sobol_last_den_inv;
	}
}


inline void mlfi_genmatrix_uniformGPUrec(
		UINT                  			seq_count,
		__constant  LoopROScalars* 		ro_scal,
		__constant  int* 				sobol_v_dir,
		            int* 				sobol_last_num_vec,
		            REAL* 				md_zd
) {
	UINT ell, j;

    UINT c = seq_count;
    ell = 0;
    while(c & 1) {
        ell++;
        c >>= 1;
    }
	for(j=0; j < ro_scal->sobol_dim; j++) {
	
		sobol_last_num_vec[j] ^= //sobol_last_num_vec[i] ^
				sobol_v_dir[ell*ro_scal->sobol_dim + j];
		md_zd[j] = sobol_last_num_vec[j] * ro_scal->sobol_last_den_inv;
	}
}


/********************************************/
/********** GAUSSIAN DISTRIBUTION ***********/
/********************************************/

#define rat_eval(a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7) \
    (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0)/ \
    (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

inline REAL small_case(REAL q) {
  REAL x = 0.180625 - q * q;
  return q * rat_eval(
                      3.387132872796366608,
                      133.14166789178437745,
                      1971.5909503065514427,
                      13731.693765509461125,
                      45921.953931549871457,
                      67265.770927008700853,
                      33430.575583588128105,
                      2509.0809287301226727,

                      1.0,
                      42.313330701600911252,
                      687.1870074920579083,
                      5394.1960214247511077,
                      21213.794301586595867,
                      39307.89580009271061,
                      28729.085735721942674,
                      5226.495278852854561);
}

inline REAL intermediate(REAL r) {
  REAL x = r - 1.6;
  return rat_eval(
                  1.42343711074968357734,
                  4.6303378461565452959,
                  5.7694972214606914055,
                  3.64784832476320460504,
                  1.27045825245236838258,
                  0.24178072517745061177,
                  0.0227238449892691845833,
                  7.7454501427834140764e-4,

                  1.0,
                  2.05319162663775882187,
                  1.6763848301838038494,
                  0.68976733498510000455,
                  0.14810397642748007459,
                  0.0151986665636164571966,
                  5.475938084995344946e-4,
                  1.05075007164441684324e-9);
}

inline REAL tail(REAL r) {
  REAL x = r - 5.0;
  return rat_eval(
                  6.6579046435011037772,
                  5.4637849111641143699,
                  1.7848265399172913358,
                  0.29656057182850489123,
                  0.026532189526576123093,
                  0.0012426609473880784386,
                  2.71155556874348757815e-5,
                  2.01033439929228813265e-7,

                  1.0,
                  0.59983220655588793769,
                  0.13692988092273580531,
                  0.0148753612908506148525,
                  7.868691311456132591e-4,
                  1.8463183175100546818e-5,
                  1.4215117583164458887e-7,
                  2.04426310338993978564e-15);
}

//#define BRANCH_UNOPT

#ifndef _OPTIMIZATION_BRANCH_OPT
inline void mlfi_ugaussian_Pinv_vector(REAL* p, UINT N, __local UCHAR* IA, UINT logBLOCK) {
	UINT i, UB = N; 
  
	for ( i=0; i < UB; i++ ) {
		REAL dp = p[i] - 0.5;
		if (fabs(dp) <= 0.425) { 
			p[i] = small_case(dp); 
		} else {
			REAL pp = (dp < 0.0) ? p[i] : (1.0 - p[i]);
			REAL r  = sqrt (- log(pp));
			REAL x  = (r <= 5.0) ? intermediate(r) : tail(r);
			p[i]    = (dp < 0.0) ? (0.0 - x) : x; 
		}
	}
}
#else
inline void mlfi_ugaussian_Pinv_vector(   // TILED VERSION!
    REAL* p , UINT SOBOL_DIM, __local ULONG* IA, UINT logBLOCK
) {
    UINT i;
    UCHAR ptr_beg = 0, ptr_end = SOBOL_DIM;

    if(SOBOL_DIM == 1) {
		REAL dp = *p - 0.5;
		if (fabs(dp) <= 0.425) { 
			*p = small_case(dp); 
		} else {
			REAL pp = (dp < 0.0) ? (dp + 0.5) : (0.5 - dp);
			REAL r  = sqrt (- log(pp));
			REAL x  = (r <= 5.0) ? intermediate(r) : tail(r);
			*p      = (dp < 0.0) ? (0.0 - x) : x; 
		}
    } else  {

        for( i = 0; i < SOBOL_DIM; i++ ) {
            p[i] = p[i] - 0.5;
            UINT ind = ( fabs(p[i]) <= 0.425 ) ? ptr_beg++ : --ptr_end;
            IA[ind<<logBLOCK] = i;
        }
        // compute the first branch
        for( i = 0; i < ptr_beg; i++ ) {
            UINT ind = IA[i<<logBLOCK];
            p[ ind ] = small_case( p[ ind ] ); 
        }
        //compute the second branch
        for(i=ptr_end; i<SOBOL_DIM; i++) {
            UINT ind = IA[i<<logBLOCK];
            REAL dp  = p[ ind ];
            REAL pp  = (dp < 0.0) ? (dp + 0.5) : (0.5 - dp);
            REAL r   = sqrt (- log(pp));
            REAL x   = (r <= 5.0) ? intermediate(r) : tail(r);
            p[ ind ] = (dp < 0.0) ? (0.0 - x) : x; 
        }
    }
}  
#endif

/********************************************/
/************* BROWNIAN BRIDGE **************/
/********************************************/
/*
inline void mlfi_brownianbridge_wiener_pathNoTransGPU(
		__constant LoopROScalars*	ro_scal,
		__local    REAL*			md_zd,
		UINT                        block_size
) {
	UINT m, md_dim = ro_scal->md_dim * block_size;
	REAL bb_sd0 = ro_scal->bb_sd0;

	bool cond = (ro_scal->bb_l == ro_scal->md_nb_path_dates) && (ro_scal->bb_l > 0);
	if(cond) {
		for (m=0; m < md_dim; m += block_size) {
			md_zd[m] = md_zd[m] * bb_sd0;
		}
	}
}
*/


// IMPORTANT: PUT THEM BACK, PLEASE!
inline void mlfi_brownianbridge_wiener_pathNoTransGPU(
		__constant  LoopROScalars*	ro_scal,
        __constant  int*             bb_ia,
        __constant  REAL*            bb_coefs,
		            REAL*			 md_zd,
        __local     REAL*            md_z,
		UINT                         block_size
) {

	UINT n, m, md_dim = ro_scal->md_dim * block_size, nb_path_dates = ro_scal->md_nb_path_dates;
    __constant int  *bb_li, *bb_bi, *bb_ri;
    __constant REAL *bb_rw, *bb_lw, *bb_sd;

	//if      ( ro_scal->bb_l != ro_scal->md_nb_path_dates ) return; 
	//else if ( ro_scal->bb_l <= 0 )                         return;
    bool cond = (ro_scal->bb_l == ro_scal->md_nb_path_dates) && (ro_scal->bb_l > 0);    
    if(!cond) return;  //number of rows doesn't match length of bridge

    bb_li = bb_ia    + ro_scal->bb_li_beg;
    bb_bi = bb_ia    + ro_scal->bb_bi_beg;
    bb_ri = bb_ia    + ro_scal->bb_ri_beg;
    bb_rw = bb_coefs + ro_scal->bb_rw_beg;
    bb_lw = bb_coefs + ro_scal->bb_lw_beg;
    bb_sd = bb_coefs + ro_scal->bb_sd_beg;

	for (n=0, m=0; m < md_dim; n++, m+=block_size) {
		UINT i;

        md_z [ (bb_bi[0]-1) * md_dim + m ] = bb_sd[0] * md_zd[n]; 

		for(i=1; i < nb_path_dates; i++) {
			int j = bb_li[i] - 1;
			int k = bb_ri[i] - 1;
			int l = bb_bi[i] - 1; 

			REAL wk = md_z [k*md_dim+m];
			REAL zi = md_zd[i*ro_scal->md_dim+n];
			
			md_z[l*md_dim+m] = (j == -1) ?
					bb_rw[i] * wk + bb_sd[i] * zi :
					bb_rw[i] * wk + bb_sd[i] * zi + bb_lw[i] * md_z[j*md_dim+m];
		}
	}
}

// NOTIFY CASH FLOW!
inline void trajectory_inner( 
		UINT                      model_num,      
        UINT                      block_size, 
		__constant LoopROScalars* ro_scal,  
        __constant REAL*          model_discounts,
        __local    REAL*          local_vhat,
		const UINT                contract_number,  
		const REAL                amount, 
		const UINT                date_index
) {
    UINT index = model_num * ro_scal->num_contracts + contract_number;
    index *= block_size;
    local_vhat[index] += amount * model_discounts[model_num*ro_scal->num_cash_flows+date_index];
           //pc_coefs[ro_scal->pc_discounts_beg + model_num*ro_scal->num_cash_flows+date_index];       
}

void trajectory_contract( // SIMPLE     
		UINT                      model_num, 
		__constant LoopROScalars* ro_scal,
        __constant REAL*          pc_coefs,
        UINT                      block_size, 
                   REAL*          inst_trajWF,
        __local    REAL*          local_vhat
);

__kernel void payoffGPU( 
	__constant LoopROScalars* ro_scal,      // RO SCALARS
	__constant int*           sobol_v_dir,  // RO SOBOL
    // RO Brownian Bridge
    __constant int*           bb_ia,
    __constant REAL*          bb_coefs,
    // RO MD INSTANCE DATA
    __constant REAL*          inst_coefs,
    // RO MODEL ATTRIBUTES
    __constant REAL*          pc_coefs,
    __constant UCHAR*         fix_index,
    // WO (GLOBAL) VHAT
	__global   REAL*          model_vhat,
    // LOCAL ARRAYS
    __local    REAL*          md_z,
	//__local    REAL*          inst_trajWF,
	__local    REAL*          vhat_local
) { 
	UINT k, block_size = get_local_size(0); 
	const UINT ct_init = 0; //ro_scal->sobol_count_ini;

    int  sobol_last_num_vec[15];//UINT[ro_scal->sobol_dim];
    REAL md_zd              [15]; //[ro_scal->sobol_dim];

    { // each local array is ajusted w.r.t. the local offset!
        UINT localID    = get_local_id(0); 
        md_z               = md_z               + localID;
        //inst_trajWF        = inst_trajWF        + localID;
        vhat_local         = vhat_local         + localID;
    }

    // initialize vhat_local to 0
    {   // ro_scal->num_contracts*k*block_size + i*block_size
        UINT step_outer = ro_scal->num_contracts*block_size, i;

        for( k = 0; k < ro_scal->inst_num*step_outer; k+=step_outer) {
            for( i = 0; i < step_outer; i+=block_size) {
	            vhat_local[k+i] = 0.0;    
            }
        }
    }

    // ready to go to the main loop
	UINT lb = (get_global_id (0) << ro_scal->logCHUNK);  
	UINT ub = min(lb + ro_scal->CHUNK, ro_scal->mc_iter_num);

	
	for( k = lb; k < ub; k++ ) { 
	
		////////////////////////////////////////////  
		/// FOR ALL FUNCTION CALLS PASS LOCAL_ID ///  
		/// AS PARAM; START LOOP'S FROM LOCAL_ID ///
		////////////////////////////////////////////

		// 1. random number generation phase -- 10 ms   
#ifndef _OPTIMIZATION_SOBOL_STRENGTH_RED_RECURR
	mlfi_genmatrix_uniformGPUind (
                ct_init+k,    ro_scal, 
                sobol_v_dir,  sobol_last_num_vec,  md_zd
            );
#else
        if( k > lb && k < ub-1) 	
            mlfi_genmatrix_uniformGPUrecOpt (
                fix_index[k-lb], ro_scal, sobol_v_dir, sobol_last_num_vec, md_zd
            ); 
        else if (k==lb)
            mlfi_genmatrix_uniformGPUind (
                ct_init+k,    ro_scal, 
                sobol_v_dir,  sobol_last_num_vec,  md_zd
            );
        else
            mlfi_genmatrix_uniformGPUrec (  
                ct_init+k, ro_scal, sobol_v_dir, sobol_last_num_vec, md_zd
            ); 
#endif
		// 2. generate uniform distribution -- 12-13 ms
		mlfi_ugaussian_Pinv_vector(md_zd, ro_scal->sobol_dim, (__local ULONG*)md_z, ro_scal->logBLOCK);

		// 3. brownian bridge refinement -- 2-3 ms  -- PUT THEM BACK, PLEASE!!!
		mlfi_brownianbridge_wiener_pathNoTransGPU(
                ro_scal, bb_ia, bb_coefs, md_zd, md_z, block_size
            );
		
		{ // 4. expand md_zd -- does not affect execution time
			UINT i;
			UINT LB =  ro_scal->md_dim * block_size;
			UINT UB = (ro_scal->md_dim * ro_scal->md_nb_path_dates - 1) * block_size;
			for  ( i = UB; i >= LB; i -= block_size ) {
				md_z[i] -= md_z[i - LB];
			}
		}
		
		{ // 5. compute trajectory -- about 7ms   
			UINT m, i, j, q, l;
			UINT dim        = ro_scal->md_dim;
			UINT npathdates = ro_scal->md_nb_path_dates;
			UINT dim_sq     = dim*dim, dim_paths = dim*npathdates;
            REAL* trajWF    = md_zd; 
			
			for (m = 0; m < ro_scal->inst_num; m++) {
                __constant REAL* c     = inst_coefs + ( ro_scal->inst_c_beg      + m*dim_sq    );
		        __constant REAL* vols  = inst_coefs + ( ro_scal->inst_vols_beg   + m*dim_paths );
		        __constant REAL* drift = inst_coefs + ( ro_scal->inst_drifts_beg + m*dim_paths );
		        __constant REAL* trajRO= inst_coefs + ( ro_scal->inst_trajRO_beg + m*dim       );

                for (j = 0; j < dim; j++) {
                    REAL accum = trajRO[j];
                    for( i=0; i < npathdates; i++)  {                    
                        REAL temp = 0.0;
						q   = dim*i + j;

						for (l = 0; l <= j; l++)
							temp += c[dim * j + l] * md_z[(q - j + l)*block_size]; //md_z[(q - j + l) << logBLOCK]; 

                        temp = exp(temp * vols[q] + drift[q]);
                        accum *= temp;
                        trajWF[q] = accum;
                    }
                }

                // 6. reduce each instance/model locally on the multiprocessor
                trajectory_contract( m, ro_scal, pc_coefs, block_size, trajWF, vhat_local ); // inst_trajWF
			}
		} // end 5. compute trajectory   

	} // END  for( k ) 

	barrier(CLK_LOCAL_MEM_FENCE);    


    if( get_local_id(0) == 0) {   // ro_scal->num_contracts*k*block_size + i*block_size
        UINT glob_size = get_global_size(0)/block_size;
        UINT i, j;

        for( i=0; i<ro_scal->inst_num; i++ ) {
            UINT ind1 = i*ro_scal->num_contracts*block_size;
            UINT ind2 = i*ro_scal->num_contracts* glob_size;

            for( j = 0; j < ro_scal->num_contracts; j++ ) {
                UINT ind11 = ind1 + j*block_size;
                UINT ind22 = ind2 + j* glob_size;

                for( k=ind11+1; k<ind11+block_size; k++ ) {
                    vhat_local[ind11] += vhat_local[k];
                }
                // the number of GPU groups forms the innermost dimension
                model_vhat[ ind22 + get_group_id(0) ] = vhat_local[ind11];
            }
        }         
    }

}

