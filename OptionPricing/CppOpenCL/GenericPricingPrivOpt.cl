/***************************************************/
/*** Contract Definitions Are Linked at Run time ***/
/*** & The RO-Scalar Data-Structure is included, ***/
/**  for e.g., ContractDefs/SmallContract.cl     ***/
/***************************************************/

#define TYPE REAL
#define FLAG uchar
#include "../../include/Utilities.cl"

/********************************************/
/********** SOBOL NUMER GENERATOR ***********/
/********************************************/

void mlfi_genmatrix_uniformGPUind (
                UINT                                    seq_count,
                __constant  LoopROScalars*              ro_scal,
                __constant  int*                        sobol_v_dir,
                            int*                        sobol_last_num_vec,
                            REAL*                       md_zd
) {
        UINT  j, k, gs, gv_k = 0;

        seq_count += 1;
        gs = seq_count >> 1;
        gs = seq_count ^  gs;

    UINT sob_dim = ro_scal->num_under * ro_scal->num_dates;

        for( j = 0; j < sob_dim; j++ )
                sobol_last_num_vec[j] = 0;
        for( k = 0; k < ro_scal->sobol_bits; ++k ) {
                if(gs & 1) {
            __constant int* dir_vect
                                = sobol_v_dir + k*sob_dim;
                        for( j=0; j < sob_dim; j++ ) {
                                // xor term g_k * v_k to direction i
                                sobol_last_num_vec[j] ^= dir_vect[j];
                        }
                }
        gs = gs >> 1;
        }
        for( j = 0; j < sob_dim; j++ ) {
                md_zd[j] = sobol_last_num_vec[j] * ro_scal->sob_norm_fact;
        }
}

inline void mlfi_genmatrix_uniformGPUrecOpt(
                UINT                                    f_ind,
                __constant  LoopROScalars*              ro_scal,
                __constant  int*                        sobol_v_dir,
                            int*                        sobol_last_num_vec,
                            REAL*                       md_zd
) {
    UINT j;
    UINT sob_dim = ro_scal->num_under * ro_scal->num_dates;
    f_ind *= sob_dim;
        for(j=0; j < sob_dim; j++) {
            sobol_last_num_vec[j] ^= sobol_v_dir[ f_ind + j ]; //f_ind * sob_dim
                md_zd[j]               = sobol_last_num_vec[j] * ro_scal->sob_norm_fact;
        }
}


inline void mlfi_genmatrix_uniformGPUrec(
                UINT                                    seq_count,
                __constant  LoopROScalars*              ro_scal,
                __constant  int*                        sobol_v_dir,
                            int*                        sobol_last_num_vec,
                            REAL*                       md_zd
) {
    UINT ell, j;

    UINT c = seq_count;
    ell = 0;
    while(c & 1) {
        ell++;
        c >>= 1;
    }
    UINT sob_dim = ro_scal->num_under * ro_scal->num_dates;
        for(j=0; j < sob_dim; j++) {
        
                sobol_last_num_vec[j] ^= sobol_v_dir[ell*sob_dim + j];
                md_zd[j] = sobol_last_num_vec[j] * ro_scal->sob_norm_fact; //sobol_last_den_inv;
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

inline void mlfi_ugaussian_Pinv_vector(REAL* p, UINT N, UINT logBLOCK) {
        UINT i, UB = N; 
  
        for ( i=0; i < UB; i++ ) {
                REAL dp = p[i] - 0.5;
                if (fabs(dp) <= 0.425) { 
                        p[i] = small_case(dp); 
                } else {
                        REAL pp = (dp < 0.0) ? dp + 0.5 : (0.5 - dp);
                        REAL r  = sqrt (- log(pp));
                        REAL x  = (r <= 5.0) ? intermediate(r) : tail(r);
                        p[i]    = (dp < 0.0) ? (0.0 - x) : x; 
                }
        }
}


/********************************************/
/************* BROWNIAN BRIDGE **************/
/********************************************/
inline void mlfi_brownianbridge_wiener_pathNoTransGPU(
        __constant  LoopROScalars*   ro_scal,
        __constant  int*             bb_inds,
        __constant  REAL*            bb_data,
        REAL*                        md_zd,
        __local     REAL*            md_z,
        UINT                         block_size
) {

    UINT n, m, md_dim = ro_scal->num_under * block_size, nb_path_dates = ro_scal->num_dates;
    __constant int  *bb_li, *bb_bi, *bb_ri;
    __constant REAL *bb_rw, *bb_lw, *bb_sd;

    bb_bi = bb_inds;
    bb_li = bb_inds +  nb_path_dates;
    bb_ri = bb_inds + (nb_path_dates<<1);

    bb_sd = bb_data;
    bb_lw = bb_data +  nb_path_dates;
    bb_rw = bb_data + (nb_path_dates<<1);

        for (n=0, m=0; m < md_dim; n++, m+=block_size) {
                UINT i;

        md_z [ (bb_bi[0]-1) * md_dim + m ] = bb_sd[0] * md_zd[n]; 

                for(i=1; i < nb_path_dates; i++) {
                        int j = bb_li[i] - 1;
                        int k = bb_ri[i] - 1;
                        int l = bb_bi[i] - 1; 

                        REAL wk = md_z [k*md_dim+m];
                        REAL zi = md_zd[i*ro_scal->num_under+n];
                        
                        md_z[l*md_dim+m] = (j == -1) ?
                                        bb_rw[i] * wk + bb_sd[i] * zi :
                                        bb_rw[i] * wk + bb_sd[i] * zi + bb_lw[i] * md_z[j*md_dim+m];
                }
        }
}



__kernel void payoffGPU( 
    __constant LoopROScalars* ro_scal,      // RO SCALARS
    __constant int*           sobol_v_dir,  // RO SOBOL
    __constant UCHAR*         fix_index,
    // RO Brownian Bridge
    __constant int*           bb_ia,
    __constant REAL*          bb_coefs,
    // RO MODELS DATA
    __constant REAL*          model_coefs, 
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
        md_z            = md_z               + localID;
        //inst_trajWF     = inst_trajWF        + localID;
        vhat_local      = vhat_local         + localID;
    }

    // initialize vhat_local to 0
    for( k = 0; k < ro_scal->num_models*block_size; k+=block_size) {
        vhat_local[k] = 0.0;    
    }

    // ready to go to the main loop
    UINT lb = (get_global_id (0) << ro_scal->log_chunk);  
    UINT ub = min(lb + ro_scal->chunk, ro_scal->num_gpuits);

        
    for( k = lb; k < ub; k++ ) { 
        
        // 1. random number generation phase 
#ifndef _OPTIMIZATION_SOBOL_STRENGTH_RED_RECURR
        mlfi_genmatrix_uniformGPUind (
                ct_init+k,    ro_scal, 
                sobol_v_dir,  sobol_last_num_vec,  md_zd
            );
#else
        if( k > lb && k < ub-1) {
            mlfi_genmatrix_uniformGPUrecOpt (
                fix_index[k-lb], ro_scal, sobol_v_dir, sobol_last_num_vec, md_zd
            ); 
        } else if (k==lb) {
            mlfi_genmatrix_uniformGPUind (
                ct_init+k,    ro_scal, 
                sobol_v_dir,  sobol_last_num_vec,  md_zd
            );
        } else {
            mlfi_genmatrix_uniformGPUrec (  
                ct_init+k, ro_scal, sobol_v_dir, sobol_last_num_vec, md_zd
            ); 
        }
#endif


        // 2. generate uniform distribution
        mlfi_ugaussian_Pinv_vector(md_zd, ro_scal->num_under*ro_scal->num_dates, ro_scal->logBLOCK);

        // 3. brownian bridge refinement 
        mlfi_brownianbridge_wiener_pathNoTransGPU(
                ro_scal, bb_ia, bb_coefs, md_zd, md_z, block_size
            );

        { // 4. expand md_zd 
            UINT i;
            UINT LB =  ro_scal->num_under * block_size;
            UINT UB = (ro_scal->num_under * ro_scal->num_dates - 1) * block_size;
            for  ( i = UB; i >= LB; i -= block_size ) {
                md_z[i] -= md_z[i - LB];
            }
        }

                { // 5. compute trajectory
                        UINT m, i, j, q, l;
                        UINT num_under  = ro_scal->num_under;
            UINT num_dates  = ro_scal->num_dates;
            REAL* trajWF    = md_zd;
            
            __constant REAL *md_c, *md_vols, *md_drifts, *md_starts, *md_discts, *md_detvals;
            { 
                UINT num_mods  = ro_scal->num_models;
                UINT dim_paths = num_under*num_dates;
                UINT offset = 0;
                md_c      = model_coefs + offset; offset += num_mods * num_under * num_under;
                md_vols   = model_coefs + offset; offset += num_mods * dim_paths;
                md_drifts = model_coefs + offset; offset += num_mods * dim_paths;
                md_starts = model_coefs + offset; offset += num_mods * num_under;
                md_discts = model_coefs + offset; offset += num_mods * ro_scal->num_cash_flows;
                md_detvals= model_coefs + offset; 
            }

                        for (m = 0; m < ro_scal->num_models; m++) {

                for (j = 0; j < num_under; j++) {
                    REAL accum = md_starts[j];
                    for( i = 0; i < num_dates; i ++)  {
                        REAL temp = 0.0;
                        q   = num_under*i + j;

                        for (l = 0; l <= j; l++)
                            temp += md_c[num_under * j + l] * md_z[(q - j + l)*block_size]; //md_z[(q - j + l) << logBLOCK]; 

                        temp = exp(temp * md_vols[q] + md_drifts[q]);
                        accum *= temp;
                        trajWF[q] = accum;
                    }
                }

                // 6. reduce each instance/model locally on the multiprocessor
                payoffFunction( m, num_under, 
                                ro_scal->num_cash_flows,  
                                ro_scal->num_det_pricers, 
                                md_discts, 
                                md_detvals, 
                                trajWF, 
                                vhat_local
                        );  

                if ( m+1 < ro_scal->num_models ) { 
                    UINT dim_paths = num_under*num_dates; 
                    md_c      += num_under * num_under;   md_vols   += dim_paths; 
                    md_drifts += dim_paths;               md_starts += num_under;
                    md_discts += ro_scal->num_cash_flows; md_detvals += ro_scal->num_det_pricers;
                }
                        }
                } // end 5. compute trajectory   

        } // END  for( k ) 

        barrier(CLK_LOCAL_MEM_FENCE);    

    segm_scan_reg_block ( vhat_local - get_local_id(0), block_size * ro_scal->num_models, block_size );
    if( get_local_id(0) == 0) {
#if 0
        for( int jj = 0; jj < ro_scal->num_models*block_size; jj+=block_size ) { 
            for( int ii = 0; ii < block_size-1; ii++ ) {
	        vhat_local[jj+block_size-1] += vhat_local[jj+ii];
	    }
	}
#endif
        UINT glob_size = get_global_size(0)/block_size;
        for ( k = 0; k < ro_scal->num_models; k++ ) {
            model_vhat[ k * glob_size + get_group_id(0) ] = 
                vhat_local[ (k + 1)*block_size - 1 ];
        }
    }
}
