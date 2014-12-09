/************************************/
/************* MACROS ***************/
/************************************/

#define TYPE REAL
#define FLAG uchar
#include "../../include/Utilities.cl"

#define DUMMY 1

/************************************/
/************* MACROS ***************/
/************************************/
//#include "Constants.h"
//
//#define underlyings(i,j)   \
//		(inst_traj[(i*ro_scal->md_dim + j)])  
//
//#define DUMMY 1
 
/********************************************/
/********** SOBOL NUMER GENERATOR ***********/
/********************************************/

// TO DO: if we want to overlap kernels, then, since we reuse md_zd for traj_wf multiply also with *ro_scal->inst_num



__kernel void mlfi_genmatrix_uniform2 ( 
        __constant LoopROScalars* 		ro_scal,
        __constant UCHAR*               sobol_fix_ind,
        __global int* 				    sobol_v_dir,
        __global   REAL*                md_zd
) {
    UINT i, j, k;    
    UINT seq_count  = ( get_global_id (0) << ro_scal->logCHUNK );
    
    md_zd = md_zd +  (get_global_id (0) << ro_scal->logCHUNK)*ro_scal->sobol_dim;
    __global REAL* md_zd_tmp = md_zd;
 
    //UINT UB = min(seq_count + (1 << ro_scal->logCHUNK), ro_scal->mc_iter_num);
    UINT UB = (seq_count + (1 << ro_scal->logCHUNK) < ro_scal->mc_iter_num) ?
               seq_count + (1 << ro_scal->logCHUNK) : ro_scal->mc_iter_num  ;
    UB = (UB > seq_count) ? UB - seq_count : 0;

    //seq_count += ro_scal->sobol_count_ini;
 
    if( UB > 0 ) {
        UINT gs;
        // Compute the random number under the INDEPENDENT formulas
        for( i = 0, j=0; j < ro_scal->sobol_dim; j++, i ++ ) {
            UINT accum = 0;
            md_zd_tmp = md_zd;
 
            // FIRST ITER COMPUTED INDEPENDENTLY!
            gs = seq_count+1 + ro_scal->sobol_count_ini; 
            gs = gs ^ (gs>>1);
            for(k=0; k<ro_scal->sobol_bit_count; k++) {
                if(gs & 1) {
                    accum ^= sobol_v_dir[ j*ro_scal->sobol_bit_count + k ];
                }
                gs = gs >> 1;
            }            
            md_zd_tmp[ i ] = ro_scal->sobol_last_den_inv*accum;

            // THE REST OF CHUNK-1 ITERATIONS COMPUTED UNDER RECURRENT FORMULA
            for( k = 1; k < UB-1; k++ ) {
                accum = accum ^ sobol_v_dir[ j*ro_scal->sobol_bit_count + sobol_fix_ind[k] ];
                
                md_zd_tmp += ro_scal->sobol_dim;
                md_zd_tmp[ i ] = ro_scal->sobol_last_den_inv*accum;
            }
            
            // last one
            gs = (seq_count + UB - 1 + ro_scal->sobol_count_ini) >> ro_scal->logCHUNK;
            UINT ell = ro_scal->logCHUNK;
            while(gs & 1) {
                ell++;
                gs >>= 1;
            }
            accum = accum ^ sobol_v_dir[j*ro_scal->sobol_bit_count + ell];
            md_zd_tmp += ro_scal->sobol_dim;
            md_zd_tmp[ i ] = ro_scal->sobol_last_den_inv*accum;
        }
    }   // end IF (seq_count < UB)

}



/********************************************/
/********** INV GAUSSIAN DISTRIB ************/
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


__kernel void mlfi_ugaussian_Pinv_vector1(  // UNTILED VERSION
    __global REAL* p , __constant LoopROScalars* ro_scal, __local REAL* tmp_p 
) {
    UINT i, UB, id = get_global_id(0) << ro_scal->logCHUNK;
    p += id * ro_scal->sobol_dim;

    UB = min(id + (1<<ro_scal->logCHUNK), ro_scal->mc_iter_num);
    UB = (UB > id) ? UB - id : 0;
    UB *= ro_scal->sobol_dim;

    for( i=0; i < UB; i++ ) {
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

/********************************************/
/************* BROWNIAN BRIDGE **************/
/********************************************/

__kernel void mlfi_brownianbridge_wiener_path1(
		__constant LoopROScalars*	ro_scal,
        __constant int*             bb_ia,
        __constant REAL*            bb_coefs,
		__global    REAL*			md_zd,
        __global    REAL*           md_z
) {
    __constant int  *bb_li, *bb_bi, *bb_ri;
    __constant REAL *bb_rw, *bb_lw, *bb_sd;

    ULONG cur_it = get_global_id (0) << ro_scal->logCHUNK;
    UINT m, i, md_dim = ro_scal->md_dim, nb_path_dates = ro_scal->md_nb_path_dates;
    UINT sobol_dim = ro_scal->sobol_dim;

    UINT UB = ( (cur_it + (1<<ro_scal->logCHUNK) )  < ro_scal->mc_iter_num) ? 
              (  cur_it + (1<<ro_scal->logCHUNK) ) : ro_scal->mc_iter_num; 

	if ( ro_scal->bb_l != ro_scal->md_nb_path_dates || ro_scal->bb_l <= 0) return; 

    bb_li = bb_ia    + ro_scal->bb_li_beg;
    bb_bi = bb_ia    + ro_scal->bb_bi_beg;
    bb_ri = bb_ia    + ro_scal->bb_ri_beg;
    bb_rw = bb_coefs + ro_scal->bb_rw_beg;
    bb_lw = bb_coefs + ro_scal->bb_lw_beg;
    bb_sd = bb_coefs + ro_scal->bb_sd_beg;
 
    { 
        UINT tmp = cur_it*sobol_dim;
        md_zd = md_zd + tmp;
        md_z  = md_z  + tmp;
    }

    for( ; cur_it < UB; cur_it++ ) {
        for (m=0; m < md_dim; m++) {
            md_z[ (bb_bi[0]-1) * md_dim + m ] = bb_sd[0] * md_zd[m]; 

    		for(i=1; i < nb_path_dates; i++) {
    			int j = bb_li[i] - 1;
    			int k = bb_ri[i] - 1;
    			int l = bb_bi[i] - 1; 

    			REAL wk = md_z [k*md_dim+m];
    			REAL zi = md_zd[i*md_dim+m];
			
    			md_z[l*md_dim+m] = (j == -1) ?
					bb_rw[i] * wk + bb_sd[i] * zi :
					bb_rw[i] * wk + bb_sd[i] * zi + bb_lw[i] * md_z[j*md_dim+m];
    		}

        }
        md_z  += sobol_dim;
        md_zd += sobol_dim;
    }
}

/**************************************************/
/************* TRAJECTORY COMPUTATION *************/
/**************************************************/

void trajectory_contract(           
		UINT                      model_num,   
		__constant LoopROScalars* ro_scal,
        __constant REAL*          pc_coefs,
        __global    REAL*         inst_traj,
        __global    REAL*         vhat
);


__kernel void mlfi_comp_traj1( 
	__constant LoopROScalars* ro_scal,      // RO SCALARS
    // RO MD INSTANCE DATA
    __constant REAL*          inst_coefs,
    // RO MODEL ATTRIBUTES
    __constant REAL*          pc_coefs,
    __local    REAL*          buff,
    // RO (GLOBAL) MD_Z
	__global   REAL*          md_z,
    // WO (GLOBAL) traj
    __global   REAL*          trajWF
) {    
    UINT dim       = ro_scal->md_dim;
    UINT dim_sq    = dim*dim;
    UINT dim_paths = dim*ro_scal->md_nb_path_dates;
    ULONG cur_it = get_global_id (0) << ro_scal->logCHUNK;

    UINT UB = (cur_it+(1<<ro_scal->logCHUNK) < ro_scal->mc_iter_num) ? cur_it + (1<<ro_scal->logCHUNK) : ro_scal->mc_iter_num; 

    trajWF += cur_it*dim_paths*ro_scal->inst_num;
    md_z   += cur_it*dim_paths;
    
    for( ; cur_it < UB; cur_it++ )    { // ...   
        UINT m, i, j, l;
			 
        { // 4. expand md_zd -- does not affect execution time
            UINT LB = dim;
            for( i = (dim_paths - 1); i >= LB; i -- ) {
                md_z[i] -= md_z[i - LB];
            }
        }

        for (m = 0; m < ro_scal->inst_num; m++) {
            __constant REAL* c     = inst_coefs + ( ro_scal->inst_c_beg      + m*dim_sq    );
            __constant REAL* vols  = inst_coefs + ( ro_scal->inst_vols_beg   + m*dim_paths );
            __constant REAL* drift = inst_coefs + ( ro_scal->inst_drifts_beg + m*dim_paths );
            __constant REAL* trajRO= inst_coefs + ( ro_scal->inst_trajRO_beg + m*dim       );
				
            for (i = 0; i < ro_scal->md_nb_path_dates; i++)  {
                UINT ind = dim*i;
                for (j = 0; j < dim; j++) {
                    REAL temp = 0.0;

                    for (l = 0; l <= j; l++) {
                        temp += c[dim * j + l] * md_z[ (ind + l) ];
                    }
                    temp = exp(temp * vols[ind + j] + drift[ind + j]);

                    trajWF[(ind+j)] = (i > 0) ? trajWF[(ind + j - dim)] * temp : trajRO[ j ] * temp ;
                }
            }

            trajWF += dim_paths; 
        } // end for m

        md_z += dim_paths;

    } // end for cur_it
}




__kernel void mlfi_reduction_step1( 
	__constant LoopROScalars* ro_scal,      // RO SCALARS
    // RO MD INSTANCE DATA
    __constant REAL*          inst_coefs,
    // RO MODEL ATTRIBUTES
    __constant REAL*          pc_coefs,
    // RO (GLOBAL) traj
    __global   REAL*          inst_trajWF,
    // WF (GLOBAL) traj
    __global   REAL*          vhat
) { 
    UINT j;
    ULONG cur_it = get_global_id (0) << ro_scal->logCHUNK;
    
    //UINT UB = min(cur_it + (1<<ro_scal->logCHUNK), ro_scal->mc_iter_num);
    UINT UB = ( cur_it + (1<<ro_scal->logCHUNK) < ro_scal->mc_iter_num) ?
                cur_it + (1<<ro_scal->logCHUNK) : ro_scal->mc_iter_num  ;

    { // adjust the start index
        UINT id = get_global_id (0);
        vhat        += ( (id >> logWARP) << logWARP )*ro_scal->num_contracts*ro_scal->inst_num  + (id & (WARP-1));
        inst_trajWF += cur_it*ro_scal->sobol_dim*ro_scal->inst_num;
    }
    { // init vhat
        UINT step_outer = ro_scal->num_contracts<<logWARP, i;
        for( i = 0; i < ro_scal->inst_num*step_outer; i+=step_outer) {
            for( j = 0; j < step_outer; j+=WARP) {
	            vhat[j+i] = 0.0;
            }
        }
    }

    for( ; cur_it < UB; cur_it++ ) {
        __global REAL* vhat_loc = vhat;
        for (j = 0; j < ro_scal->inst_num; j++) {   
            trajectory_contract( j, ro_scal, pc_coefs, inst_trajWF, vhat_loc ); 
            inst_trajWF += ro_scal->sobol_dim;   
            vhat_loc    += (ro_scal->num_contracts << logWARP); 
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);    

    if( get_local_id(0) == 0) {
#if 1
        UINT step_inner = ro_scal->num_contracts << logWARP;
        UINT step_outer = step_inner * ro_scal->inst_num;
        UINT BD         = (step_outer >> logWARP) * ro_scal->BLOCK;
        UINT i, k;
        for( i = 0; i < BD; i += step_outer ) {
            for( j = 0; j < step_outer; j += step_inner ) {
                for( k = j; k < j+step_inner; k += WARP )       {
                    for( cur_it = k+(i==0? 1 : 0); cur_it < k+WARP; cur_it++ )  {
                        vhat[k] += vhat[i+cur_it];
                    }
                }
            }
        }
#else
        UINT step_j = ro_scal->num_contracts<<logWARP;
        UINT BD_j   = ro_scal->inst_num*step_j;
        for (j = 0; j < BD_j; j+=step_j) { // for all instances
            for(cur_it=0; cur_it < step_j; cur_it+=WARP) { // for all contracts
                UINT i;
                for(i=1; i<WARP; i++) { // reduce WARP its
                    vhat[ j+cur_it ] += vhat[ j+cur_it+i ]; 
                }
            }
        }
#endif
    }

}

  
// NOTIFY CASH FLOW!
inline void trajectory_inner( 
		UINT                      model_num,
		__constant LoopROScalars* ro_scal,  
        __constant REAL*          model_discounts,
        __global   REAL*          vhat,
		const UINT                contract_number,  
		const REAL                amount, 
		const UINT                date_index
) {
    vhat[(contract_number << logWARP)] += amount * model_discounts[model_num*ro_scal->num_cash_flows+date_index];
           //pc_coefs[ro_scal->pc_discounts_beg + model_num*ro_scal->num_cash_flows+date_index];       
}

