#include "CosminSobol.h"
#include "GPU_Trimmed_DS.h"

#include <string.h>
#include <stdio.h>

#include <sys/time.h>
#include <time.h>

#include "TimeHelperChr.h"

#define DEBUG_ON 0

double RAND_SUM = 0.0;

/************************************************/
/************************************************/
/********* Re-Structure/Filter Data BEG *********/
/************************************************/
/************************************************/

void setCHUNK(LoopROScalars* ret) {
	ret->logCHUNK = 4;
	ret->CHUNK    = 1 << ret->logCHUNK;
}

//#define MAX_DATES 12
#define MAX_DATES 10000
void buildROscalars(
		LoopROScalars*               ret,
		UINT                         numiter,
		void*                        data,
		UINT                         num_models,
		mlfi_nmc_model**             model,
		const mlfi_nmc_pricing_code* pc
) {
	driver *md = data;

	setCHUNK(ret);
	ret->logBLOCK = 7;
	ret->BLOCK    = 1 << ret->logBLOCK;
#ifdef _OPTIMIZATION_TILE
        ret->TILE_FACT= TILE;
#else
        ret->TILE_FACT= 16;
#endif


	{ // fill number of iterations
		ret->mc_iter_num     = numiter;
		ret->TOT_MC_ITER_NUM = ret->mc_iter_num;

		UINT cur_iter = 0, num_iter = numiter/ret->CHUNK;
		if(numiter % ret->CHUNK != 0)  num_iter++;
		ret->mc_bigit_num = num_iter;
#if DEBUG_ON
//      printf("Monte-Carlo Iters: %u, BIG ITERS: %u\n\n", ret->mc_iter_num, ret->mc_bigit_num);
#endif
    }

	{ // fill md/rng fields
		assert(md->nb_path_dates == md->rng->nb_dates && "md->nb_path_dates != md->rng->nb_dates !!!\n");
		UINT my_num_dates = (md->nb_path_dates > MAX_DATES) ? MAX_DATES : md->nb_path_dates;
		ret->md_dim           = md->dim;
		ret->md_nb_path_dates = my_num_dates; //md->nb_path_dates;
		ret->rng_dim          = md->rng->dim;
		ret->rng_nb_dates     = my_num_dates; //md->rng->nb_dates;
	}

# if(DEBUG_ON)
//	printf("md_dim: %d, nb_path_dates: %d, rng_dim: %d, rng_nb_dates: %d\n",
//          ret->md_dim, ret->md_nb_path_dates,
//          ret->rng_dim, ret->rng_nb_dates
//        );
#endif
	{ // fill sobol fields
		sobol_state_t* s_state;
		s_state   = (sobol_state_t *) md->rng->sobol->state;

		ret->sobol_count_ini    = s_state->sequence_count;

		ret->sobol_dim          = ret->md_dim * ret->md_nb_path_dates; //md->rng->sobol->dimension;
		assert(md->dim*md->nb_path_dates == md->rng->sobol->dimension &&
				"ERROR: md->dim*md->nb_path_dates != ro_scals->sobol_dim !!!");

		ret->sobol_last_den_inv = (REAL) s_state->last_denominator_inv;
		ret->sobol_bit_count    = SOBOL_BIT_COUNT;
	}

# if(DEBUG_ON)
//    printf("sobol_count_ini: %d , sobol_last_den_inv: %.32f, sobol_bit_count: %d\n",
//        ret->sobol_count_ini, ret->sobol_last_den_inv, ret->sobol_bit_count
//    );
#endif

	// brownian bridge
	ret->bb_l = md->bb->l;
	//ret->bb_sd0   = (REAL) md->bb->sd[0]; // promoted to array

	// model instance number
	ret->inst_num = md->nb_instances;
	assert(ret->inst_num == num_models && "md->nb_instances != num_models is buildROscalars!!!");


	{ // pricer attributes number
		ret->num_contracts   = pc->nb_contracts;
		ret->num_det_pricers = pc->nb_deterministic_pricers;
		ret->num_cash_flows  = pc->nb_cash_flows;
	}

	{ // model
		//ret->model_deter_val0 = (REAL) model->deterministic_values[0];
		//ret->model_discounts0 = (REAL) ((instance*)model->model_data)->discounts[0];
	}

# if(DEBUG_ON)
//	printf("inst_num: %d , bb_l: %d \n", ret->inst_num, ret->bb_l);
//	printf("NUM_CONTRACTS: %d, NUM_DET_PRICERS: %d, NUM_CASH_FLOWS: %d\n",
//			pc->nb_contracts, pc->nb_deterministic_pricers, pc->nb_cash_flows);
    printf("%u // number of Monte-Carlo Iterations\n", ret->mc_iter_num);
    printf("%u       // number of path dates\n%u       // number of underlyings\n", 
            ret->md_nb_path_dates, ret->md_dim);
    printf("%u      // integer bit-length representation for Sobol\n\n", 
            ret->sobol_bit_count );

#endif
	//printf("inst_num: %d , bb_l: %d, bb_sd0: %.32f, model_deter_val0: %.32f, model_discount0: %.32f \n\n\n",
	//		ret->inst_num, ret->bb_l, ret->bb_sd0, ret->model_deter_val0, ret->model_discounts0
	//	);


	{ // Indexes:
		// Brownian Bridge
		ret->bb_li_beg = 0;
		ret->bb_bi_beg = ret->md_nb_path_dates;
		ret->bb_ri_beg = ret->md_nb_path_dates + ret->md_nb_path_dates;
		ret->bb_rw_beg = 0;
		ret->bb_lw_beg = ret->md_nb_path_dates;
		ret->bb_sd_beg = ret->md_nb_path_dates + ret->md_nb_path_dates;

		// MD instance data
		UINT dim_sq        = ret->md_dim*ret->md_dim;
		UINT dim_path      = ret->md_dim*ret->md_nb_path_dates;
		UINT inst_dim_path = ret->inst_num*dim_path;

		ret->inst_c_beg      = 0;
		ret->inst_vols_beg   = ret->inst_num*dim_sq;
		ret->inst_drifts_beg = ret->inst_vols_beg + inst_dim_path;
		ret->inst_trajRO_beg = ret->inst_drifts_beg + inst_dim_path;

		// PC data
		ret->pc_deter_val_beg = 0;
		ret->pc_discounts_beg = ret->inst_num * ret->num_det_pricers;
	}
}

#define WITH_ASSERTS 0

void buildROarrays(
		struct LoopROArrays*   ret,
		LoopROScalars*         ro_scal,
		void*                  data,
		mlfi_nmc_model**       models
) {
	UINT k, m, i, j;
	driver *md = data;

	if(WITH_ASSERTS) {
		assert(md->dim == ro_scal->md_dim && "ERROR: md->dim != ro_scals->md_dim !!!");
		assert(md->nb_instances == ro_scal->inst_num && "ERROR: md->nb_instances != ro_scals->inst_num !!!");
		assert(md->nb_path_dates == ro_scal->md_nb_path_dates && "ERROR: md->nb_path_dates != ro_scals->md_nb_path_dates !!!");
		assert(md->rng->sobol->dimension == ro_scal->sobol_dim && "ERROR: md->rng->sobol->dimension != ro_scals->sobol_dim !!!");
	}

	{ // sobol's v_direction
		sobol_state_t* s_state = (sobol_state_t *) md->rng->sobol->state;

//		if(DEBUG_ON) printf("v_dir: [sobol_bit_count]x[sobol_dim]: \n");
		ret->sobol_v_dir = (int*) malloc(sizeof(int)*SOBOL_BIT_COUNT*md->rng->sobol->dimension);
		for(k=0; k<SOBOL_BIT_COUNT; k++) {
			int* row_w = ret->sobol_v_dir + k*md->rng->sobol->dimension;
			int* row_r = s_state->v_direction[k];

//			if(DEBUG_ON) printf("\t{ ");

			for(i=0; i<ro_scal->sobol_dim; i++) { // md->rng->sobol->dimension
//				if(DEBUG_ON)  {
//					printf("%d", row_r[i]); if(i!=ro_scal->sobol_dim-1) printf(",");
//				}
				row_w[i] = row_r[i];
			}

//			if(DEBUG_ON) printf(" }\n ");
		}
//		if(DEBUG_ON)  printf("\n\n");

        if(DEBUG_ON) 
            printf("// Begin Direction vectors [%u*%u][%u]\n[", 
                    ro_scal->md_dim, ro_scal->md_nb_path_dates, 
                    ro_scal->sobol_bit_count );

		// sobol_v_dir transposed!
		ret->sobol_v_dir_t = (int*) malloc(sizeof(int)*SOBOL_BIT_COUNT*md->rng->sobol->dimension);
        for(i=0; i<ro_scal->sobol_dim; i++) {
            if(DEBUG_ON) printf("\n\t[\n\t\t");

            for(k=0; k<SOBOL_BIT_COUNT; k++) {
				ret->sobol_v_dir_t[i*SOBOL_BIT_COUNT + k] = 
                    ret->sobol_v_dir[k*ro_scal->sobol_dim+i];

                if(DEBUG_ON) {
                    printf("%d", ret->sobol_v_dir_t[i*SOBOL_BIT_COUNT + k]);
                    if (k < SOBOL_BIT_COUNT-1) printf(", ");
                }
			}

            if(DEBUG_ON) {
                printf("\n\t]");
                if (i < ro_scal->sobol_dim-1) printf(",");
            }
		}

        if(DEBUG_ON) 
            printf("\n] // End   Direction Vectors\n\n\n"); 


		ret->sobol_fix_ind = NULL;
	}


	{ // md per instance data
#if 0
		ret->inst_c      = (REAL*) malloc(md->nb_instances*(md->dim*md->dim          )*sizeof(REAL));
		ret->inst_drifts = (REAL*) malloc(md->nb_instances*(md->dim*md->nb_path_dates)*sizeof(REAL));
		ret->inst_vols   = (REAL*) malloc(md->nb_instances*(md->dim*md->nb_path_dates)*sizeof(REAL));
		ret->inst_trajRO = (REAL*) malloc(md->nb_instances*(md->dim                  )*sizeof(REAL));

		ret->model_deter_val = (REAL*) malloc(md->nb_instances * ro_scal->num_det_pricers*sizeof(REAL));
		ret->model_discounts = (REAL*) malloc(md->nb_instances * ro_scal->num_cash_flows *sizeof(REAL));
#endif
		ret->inst_coefs = (REAL*) malloc( (ro_scal->inst_trajRO_beg  + md->nb_instances*md->dim)*sizeof(REAL) );
		ret->pc_coefs   = (REAL*) malloc( (ro_scal->pc_discounts_beg + md->nb_instances*ro_scal->num_cash_flows)*sizeof(REAL) );

		ret->inst_coefs[ro_scal->inst_trajRO_beg] = 0.0; //ret->inst_trajRO[0] = 0.0;

        assert(ro_scal->inst_num == 1 && "ERROR: Number of Instances != 1 !!!");

		for(m=0; m<ro_scal->inst_num; m++) { // md->nb_instances
			instance *inst = md->instances[m];

//          if(DEBUG_ON)  printf("\n\n\nFOR INSTANCE: %d\n\n", m);
            printf("\n\n// The Market / Model Data\n//\t md_c\t md_vols\t md_drifts\t md_starts\t md_detval\t md_disc\n\n");

			{ // inst_c & trajectory
				REAL* row_w_c = ret->inst_coefs + (ro_scal->inst_c_beg + (md->dim*md->dim)*m );
				//REAL* row_w_c = ret->inst_c + (md->dim*md->dim)*m;
				double* row_r_c = inst->c;

//                if(DEBUG_ON)  printf("inst_c: [%d]x[%d] IS: \n", ro_scal->md_dim, ro_scal->md_dim);
                if(DEBUG_ON)  printf("[\t // md_c[%u][%u] (begin)\n", 
                                        ro_scal->md_dim, ro_scal->md_dim);

                int tot_size = ro_scal->md_dim*ro_scal->md_dim;
				for(i=0; i<tot_size; i++) { // md->dim*md->dim
					row_w_c[i] = (REAL) row_r_c[i];
					if(DEBUG_ON) {
						if(i%ro_scal->md_dim == 0) printf("\t[ ");
						printf("%.7f ", row_w_c[i]);
						if(i%ro_scal->md_dim == (ro_scal->md_dim-1)) {
                            if (i == tot_size-1)  printf(" ]\n");
                            else                  printf(" ],\n");
                        } else printf(", ");
					}
				}
            
                if(DEBUG_ON)  printf("]\t // md_c (end)\n\n");
			}

			{ // inst_drifts & inst_vols
				REAL* row_w_drifts = ret->inst_coefs + ( ro_scal->inst_drifts_beg + (md->dim*md->nb_path_dates)*m );
				REAL* row_w_vols   = ret->inst_coefs + ( ro_scal->inst_vols_beg   + (md->dim*md->nb_path_dates)*m );
				//REAL* row_w_drifts = ret->inst_drifts + (md->dim*md->nb_path_dates)*m;
				//REAL* row_w_vols   = ret->inst_vols   + (md->dim*md->nb_path_dates)*m;

				double* row_r_drifts = inst->drift;
				double* row_r_vols   = inst->vols;

//              if(DEBUG_ON)  printf("inst_vols: [%d]x[%d] IS: \n", ro_scal->md_dim, ro_scal->md_dim);
                if(DEBUG_ON)  printf("[\t // md_vols[%u][%u] volatility (begin)\n", 
                                        ro_scal->md_nb_path_dates, ro_scal->md_dim);

                int tot_size = ro_scal->md_dim*ro_scal->md_nb_path_dates;

				for(i=0; i<tot_size; i++) {  // md->dim*md->nb_path_dates
					row_w_vols  [i] = (REAL) row_r_vols  [i];

					if(DEBUG_ON) {
						if(i%ro_scal->md_dim == 0) printf("\t[ ");
						if(DEBUG_ON) printf("%.7f ", row_w_vols[i]);
						if(i%ro_scal->md_dim == (ro_scal->md_dim-1)) {
                            if (i == tot_size-1)  printf(" ]\n");
                            else                  printf(" ],\n");
                        } else printf(", ");
					}
				}

//              if(DEBUG_ON)  printf("inst_drifts: [%d]x[%d] IS: \n", ro_scal->md_dim, ro_scal->md_nb_path_dates);
                if(DEBUG_ON)  printf("]\t // md_vols volatility (end)\n\n");

                if(DEBUG_ON)  printf("[\t // md_drifts[%u][%u] (begin)\n", 
                                        ro_scal->md_nb_path_dates, ro_scal->md_dim);

				for(i=0; i<tot_size; i++) {  // md->dim*md->nb_path_dates
					row_w_drifts[i] = (REAL) row_r_drifts[i];

					if(DEBUG_ON) {
						if(i%ro_scal->md_dim == 0) printf("\t[ ");
						if(DEBUG_ON) printf("%.16f ", row_r_drifts[i]);
						if(i%ro_scal->md_dim == (ro_scal->md_dim-1)) {
                            if (i == tot_size-1)  printf(" ]\n");
                            else                  printf(" ],\n");
                        }
						else                                         printf(", ");
					}

				}

                if(DEBUG_ON)  printf("]\t // md_drifts (end)\n\n");

			}

			{ // fill-in trajectory
				REAL*   traj_w = ret->inst_coefs + (ro_scal->inst_trajRO_beg+m*md->dim);
				double* traj_r = inst->trajectory - md->dim;

//              if(DEBUG_ON)  printf("instRO: [%d] IS: \n\t{", ro_scal->md_dim);
                if(DEBUG_ON)  printf("// md_starts[%u]\n[ ", ro_scal->md_dim);

				for(i=0; i<ro_scal->md_dim; i++) { // md->dim
					traj_w[i] = (REAL) traj_r[i];
					//ret->inst_trajRO[m*md->dim+i] = (REAL) traj_r[i]; //inst->trajectory[i-md->dim];
					if(DEBUG_ON) { printf("%.16f", traj_r[i]); if(i!=ro_scal->md_dim-1) printf(", "); }
				}

                if(DEBUG_ON)  printf(" ]\n\n");
			}


			{ // fill in model_deter_val
				REAL* deter_val_out = ret->pc_coefs + ( ro_scal->pc_deter_val_beg + m*ro_scal->num_det_pricers );
				//REAL*   deter_val_out = ret->model_deter_val + m*ro_scal->num_det_pricers;
				double* deter_val_in  = models[m]->deterministic_values;

				if(DEBUG_ON)  printf("// model deterministic values [%d]\n[ ", ro_scal->num_det_pricers);

				for(i=0; i<ro_scal->num_det_pricers; i++) {
					deter_val_out[i] = deter_val_in[i];
					if(DEBUG_ON) { printf("%.16f ", deter_val_in[i]); if(i != ro_scal->num_det_pricers-1) printf(", "); }
				}

				if(DEBUG_ON) printf(" ]\n\n");
			}

			{ // fill in model_discounts
				REAL*   discounts_out = ret->pc_coefs + ( ro_scal->pc_discounts_beg + m*ro_scal->num_cash_flows );
				//REAL*   discounts_out = ret->model_discounts + m*ro_scal->num_cash_flows;
				instance *inst = (instance*) models[m]->model_data;
				double* discounts_in  = inst->discounts;

//              if(DEBUG_ON)  printf("discounts: [%d] IS: \n\t{", ro_scal->num_cash_flows);
                if(DEBUG_ON)  printf("// model discounts [%u]\n[ ", ro_scal->num_cash_flows);

				for(i=0; i<ro_scal->num_cash_flows; i++) {
					discounts_out[i] = discounts_in[i];
					if(DEBUG_ON) { printf("%.16f ", discounts_in[i]); if(i != ro_scal->num_cash_flows-1) printf(", "); }
				}

				if(DEBUG_ON) printf(" ]\n\n");
			}

		}
	}



	{ // Brownian Bridge arrays:
#if 0
		ret->bb_li = (UINT*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]
		ret->bb_bi = (UINT*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]
		ret->bb_ri = (UINT*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]

		ret->bb_rw = (REAL*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]
		ret->bb_lw = (REAL*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]
		ret->bb_sd = (REAL*)malloc(ro_scal->md_nb_path_dates*sizeof(REAL)); // [ro_scal->md_nb_path_dates]
#endif
        if(DEBUG_ON) printf("// Brownian Bridge (BB) Metadata\n\n//BB indirect accessing:\n[");

		ret->bb_ia    = (UINT*)malloc(3*ro_scal->md_nb_path_dates*sizeof(REAL));
		ret->bb_coefs = (REAL*)malloc(3*ro_scal->md_nb_path_dates*sizeof(REAL));

//      if(DEBUG_ON) printf("bb_bi : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_ia[ro_scal->bb_bi_beg+i] = md->bb->bi[i];

			if(DEBUG_ON) { printf("%d", md->bb->bi[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ],\t//bb_bi[%u]\n", ro_scal->md_nb_path_dates);


//		if(DEBUG_ON) printf("bb_li : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_ia[ro_scal->bb_li_beg+i] = md->bb->li[i];

			if(DEBUG_ON) { printf("%d", md->bb->li[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ],\t//bb_li[%u]\n", ro_scal->md_nb_path_dates);


//        if(DEBUG_ON) printf("bb_ri : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_ia[ro_scal->bb_ri_beg+i] = md->bb->ri[i];

			if(DEBUG_ON) { printf("%d", md->bb->ri[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ] \t//bb_ri[%u]\n]\n\n", ro_scal->md_nb_path_dates);


        if(DEBUG_ON) printf("//BB (real) data:\n[");
//      if(DEBUG_ON) printf("bb_rw : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_coefs[ro_scal->bb_sd_beg+i] = md->bb->sd[i];

			if(DEBUG_ON) { printf("%.16f", md->bb->sd[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ],\t//bb_sd[%u] (standard deviation)\n", ro_scal->md_nb_path_dates);

//      if(DEBUG_ON) printf("bb_lw : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_coefs[ro_scal->bb_lw_beg+i] = md->bb->lw[i];

			if(DEBUG_ON) { printf("%.16f", md->bb->lw[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ],\t//bb_lw[%u]\n", ro_scal->md_nb_path_dates);

//		if(DEBUG_ON) printf("bb_rw : [%d]\n\t{ ", ro_scal->md_nb_path_dates);
        if(DEBUG_ON) printf("\t[ ");
		for(i=0; i<ro_scal->md_nb_path_dates; i++) {
			ret->bb_coefs[ro_scal->bb_rw_beg+i] = md->bb->rw[i];

			if(DEBUG_ON) { printf("%.16f", md->bb->rw[i]); if(i!=ro_scal->md_nb_path_dates-1) printf(", "); }
		}
        if(DEBUG_ON) printf(" ]\t//bb_rw[%u]\n]\n\n", ro_scal->md_nb_path_dates);
	}

}

void buildPrivateArrays(
		struct LoopPrivateArrays* ret,
		LoopROScalars*            ro_scal,
		void*                     data
) {
	UINT k, m, i, j;
	driver *md = data;
	UINT num_iter = ro_scal->mc_bigit_num;

	if(WITH_ASSERTS) {
		assert(md->dim*md->nb_path_dates == ro_scal->sobol_dim &&
					"ERROR: md->dim*md->nb_path_dates != ro_scals->sobol_dim !!!");
		assert(md->rng->sobol->dimension == ro_scal->sobol_dim &&
					"ERROR: md->rng->sobol->dimension != ro_scals->sobol_dim !!!");
	}

	{ // sobol's v_direction
		sobol_state_t* s_state = (sobol_state_t *) md->rng->sobol->state;

		ret->sobol_last_num_vec = (int*) malloc(sizeof(int)*md->rng->sobol->dimension);
		for(i=0; i<ro_scal->sobol_dim; i++) {  // md->rng->sobol->dimension
			ret->sobol_last_num_vec[i] = s_state->last_numerator_vec[i];
		}
	}

	{ // md->z / md->zd arrays
		ret->md_zd = (REAL*) malloc(ro_scal->sobol_dim*sizeof(REAL));
		ret->md_z  = (REAL*) malloc(ro_scal->sobol_dim*sizeof(REAL));
		for(k=0; k<ro_scal->sobol_dim; k++) {
			ret->md_zd[k] = (REAL) md->zd[k];
			ret->md_z [k] = (REAL) md->z [k];
		}
	}

	{ // trajectory!
		ret->inst_trajWF = (REAL*)malloc(md->nb_instances*md->dim*md->nb_path_dates*sizeof(REAL));
		for(k=0; k<ro_scal->inst_num; k++) { //md->nb_instances
			instance *inst = md->instances[k];
			double* inst_traj_r = inst->trajectory;
			REAL  * inst_traj_w = (REAL*) ret->inst_trajWF + k*md->dim*md->nb_path_dates;

			for(i=0; i<ro_scal->md_dim * ro_scal->md_nb_path_dates; i++) { //md->dim*md->nb_path_dates
				inst_traj_w[i] = (REAL) inst_traj_r[i];
			}
		}
	}

	{ // the final accumulator -- this should be a two dimensional vector
	  //     in the general case is a reuction over an array!
	  // model_vhat -- [num_big_iter]x[inst_num]x[num_contracts]
		UINT sz = ro_scal->inst_num * num_iter * ro_scal->num_contracts; //sz_vhat*num_iter*num_models;
//		if(DEBUG_ON) printf("SIZE OF VHAT: %d, numbigiter: %d\n\n\n", sz, num_iter);
		ret->model_vhat = (REAL*) malloc(sz*sizeof(REAL));
		for(i=0; i<sz; i++) {
			ret->model_vhat[i] = (REAL) 0.0;
		}

		ret->vhat_fin = (double*)malloc(ro_scal->inst_num*ro_scal->num_contracts*sizeof(double));
	}
}

void updateState(struct LoopPrivateArrays* ret, UINT seq_count, void* data) {
	UINT k, i;
	driver *md = data;
	sobol_state_t* s_state = (sobol_state_t *) md->rng->sobol->state;

	// 1. update sobol's sequence_count
	s_state->sequence_count = seq_count;
#if 0
	// 2. update all instance trajectories
	for(k=0; k< md->nb_instances; k++) {
		instance *inst = md->instances[k];
		REAL* inst_traj_w = ret->inst_trajWF + k*md->dim*md->nb_path_dates;

		for(i=0; i<md->dim*md->nb_path_dates; i++) {
			inst->trajectory[i] = inst_traj_w[i];
		}
	}
#endif
}

void updateModels(
		mlfi_nmc_model **         models,
		LoopROScalars*            ro_scal,
		struct LoopPrivateArrays* priv_arr
) {
	UINT j, k, l;

	for (j=0; j<ro_scal->inst_num; j++) {
		mlfi_nmc_model * model = models[j];
		instance* inst = (instance*)model->model_data;
		REAL* priv_vhat_row = priv_arr->model_vhat + j*(ro_scal->num_contracts * ro_scal->mc_bigit_num);

		for( k=0; k<ro_scal->num_contracts; k++) {
			REAL* priv_vhat_c = priv_vhat_row + k*ro_scal->mc_bigit_num;
			double  accum = 0.0;
			for( l=0; l<ro_scal->mc_bigit_num; l++ ) {
				accum += priv_vhat_c[l];
			}
			inst->vhat[k] = accum;
			printf("\n\nFINAL RESULT CPU FOR (inst_num,contr_num)=(%d,%d) is: %.16f !\n\n", j, k, inst->vhat[k]/ro_scal->mc_iter_num);
		}
	}
}

/************************************************/
/************************************************/
/********* Re-Structure/Filter Data END *********/
/************************************************/
/************************************************/




/************************************/
/***** Random Number Generator ******/
/************************************/

/**
 * COSMIN
 *
 *	1. WHAT CONDITIONS NEED TO HOLD FOR ||
 *   	(i)   rng->kind == Sobol
 *		(ii)  md->rng->sobel->dimension <= max_poly
 *		(iii) update md->rng->sobel->state->sequence_count at the end!
 *
 *  2. SUMMARIZATION:
 *
 *   READ-ONLY:  dim = md->rng->sobel->dimension
 *               md->rng->sobel->state->v_direction[0 : SOBOL_BIT_COUNT-1][1:dim]
 *               md->rng->sobel->state->last_denominator_inv
 *   WRITE-FIRST (PRIVATIZED):
 *               zd
 *               md->rng->sobel->state->last_numerator_vec[1:dim]
 *
 */
void mlfi_genmatrix_uniformGPUind(
		UINT seq_count,
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	UINT  i, k, gs, gv_k = 0;

	seq_count += 1;
	gs = seq_count >> 1;
	gs = seq_count ^  gs;

	for( i = 0; i < ro_scal->sobol_dim; ++i )
		priv_arr->sobol_last_num_vec[i] = 0;

	for( k = 0; k < ro_scal->sobol_bit_count; ++k ) {
		int* dir_vect = ro_arr->sobol_v_dir + k*ro_scal->sobol_dim;
		UINT gv_mod = gs % 2;
		gs = gs >> 1;

		if(gv_mod != 0) {
			assert(gv_mod == 1 && "ERROR gv_mod should be 1!");

			for( i = 0; i < ro_scal->sobol_dim; ++i ) {
				// xor term g_k * v_k to direction i
				priv_arr->sobol_last_num_vec[i] ^= dir_vect[i];
			}
		}
	}
	for(i = 0; i < ro_scal->sobol_dim; ++i) {
		priv_arr->md_zd[i] = priv_arr->sobol_last_num_vec[i] * ro_scal->sobol_last_den_inv;
		//RAND_SUM += priv_arr->md_zd[i];
	}

}

void mlfi_genmatrix_uniformGPUrec(
		UINT seq_count,
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	UINT ell = 0, i;
	UINT c = seq_count;
	while(c & 1) {
		ell++;
		c >>= 1;
	}

	if(ell > SOBOL_BIT_COUNT) { printf("FAILURE!!!\n\n\n"); exit(0); }

	for(i=0; i<ro_scal->sobol_dim; i++) {

		priv_arr->sobol_last_num_vec[i] = priv_arr->sobol_last_num_vec[i] ^
				ro_arr->sobol_v_dir[ell*ro_scal->sobol_dim + i];
		priv_arr->md_zd[i] =
				priv_arr->sobol_last_num_vec[i] * ro_scal->sobol_last_den_inv;

		//RAND_SUM += priv_arr->md_zd[i];
	}

#if 0
	if(seq_count_orig > 512 && seq_count_orig < 532) {

		double test = (double)(1 << ro_scal->sobol_bit_count);

		printf("RanNum: %d, output: %e, den_inv: %.32e for count: %d %d %.32e\n\n\n",
				priv_arr->sobol_last_num_vec[0],
				priv_arr->md_zd[0],
				ro_scal->sobol_last_den_inv,
				seq_count_orig,
				ro_scal->sobol_bit_count,
				((double)1)/((double)  (1 << ro_scal->sobol_bit_count))
			);
	}
#endif

}


/********************************************/
/********** GAUSSIAN DISTRIBUTION ***********/
/********************************************/

#define rat_eval(a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7) \
    (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0)/ \
    (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

inline static REAL small_case(REAL q) {
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

inline static REAL intermediate(REAL r) {
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

inline static REAL tail(REAL r) {
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


inline void mlfi_ugaussian_Pinv_vectorGPU(REAL* p, UINT N) {
	UINT i;

	for ( i=0; i < N; i++) {
		REAL dp = p[i] - 0.5;
		if ( fabs(dp) <= 0.425 ) {
			p[i] = small_case(dp);
		} else {
			REAL pp = (dp < 0.0) ? p[i] : (1.0 - p[i]);
			REAL r  = sqrt (- log(pp));
			REAL x  = (r <= 5.0) ? intermediate(r) : tail(r);
			p[i]    = (dp < 0.0) ? (0.0 - x) : x;
		}
	}
}


///////////////////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////

/* fill a matrix with random normal deviates */
void mlfi_genmatrix_normalGPU(
		UINT seq_count,
		UINT ind_flag,
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	if(ind_flag)
		mlfi_genmatrix_uniformGPUind( seq_count, ro_scal, ro_arr, priv_arr );
	else
		mlfi_genmatrix_uniformGPUrec( seq_count, ro_scal, ro_arr, priv_arr );

	// Cosmin summary: RW `zd' in [0,rng->nb_dates * rng->dim]
	mlfi_ugaussian_Pinv_vectorGPU(priv_arr->md_zd, ro_scal->rng_nb_dates * ro_scal->rng_dim);

}


///////////////////////////////////////////
////////////// BROWNIAN BRIDGE ////////////
///////////////////////////////////////////

// (rz, cz) == (npathdates, dim)
int mlfi_brownianbridge_wiener_pathNoTransGPU_SIMPLE(
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	UINT m;
	REAL *bb_sd = ro_arr->bb_coefs + ro_scal->bb_sd_beg;
	if( ro_scal->bb_l != ro_scal->md_nb_path_dates || ro_scal->bb_l <= 0 )
		return -1; /* number of rows doesn't match length of bridge */

	for (m=0; m < ro_scal->md_dim; m++) {
		priv_arr->md_zd[m] = priv_arr->md_zd[m] * bb_sd[0]; //ro_scal->bb_sd0;
	}

	return 0;
}

// mlfi_brownianbridge_wiener_path(md->bb, md->zd, z, npathdates, dim, npathdates, dim);
//int mlfi_brownianbridge_wiener_path(const mlfi_brownianbridge *b, enum CBLAS_TRANSPOSE trans, const matrix z, const matrix w,
//		const int rz, const int cz, const int rw, const int cw)
int mlfi_brownianbridge_wiener_pathNoTransGPU(
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	size_t m;
	REAL *bb_sd, *bb_rw, *bb_lw;
	UINT *bb_li, *bb_ri, *bb_bi;

	if      ( ro_scal->bb_l != ro_scal->md_nb_path_dates ) return -1; /* number of rows doesn't match length of bridge */
	else if ( ro_scal->bb_l <= 0 )                         return  0;

	bb_li = ro_arr->bb_ia    + ro_scal->bb_li_beg;
	bb_bi = ro_arr->bb_ia    + ro_scal->bb_bi_beg;
	bb_ri = ro_arr->bb_ia    + ro_scal->bb_ri_beg;

	bb_rw = ro_arr->bb_coefs + ro_scal->bb_rw_beg;
	bb_lw = ro_arr->bb_coefs + ro_scal->bb_lw_beg;
	bb_sd = ro_arr->bb_coefs + ro_scal->bb_sd_beg;

	//printf("Cosmin Bridge, (cz, rz): (%d, %d)\n", cz, rz);
	for (m=0; m < ro_scal->md_dim; m++) {
		size_t i;

		priv_arr->md_z[ (bb_bi[0]-1) * ro_scal->md_dim + m ] = bb_sd[0] * priv_arr->md_zd[m]; // z[m]  needs copying?

		for(i=1; i < ro_scal->md_nb_path_dates; i++) {
			int j = bb_li[i] - 1;
			int k = bb_ri[i] - 1;
			int l = bb_bi[i] - 1;

			double wk = priv_arr->md_z [k*ro_scal->md_dim+m];
			double zi = priv_arr->md_zd[i*ro_scal->md_dim+m];

			priv_arr->md_z[l*ro_scal->md_dim+m] = (j == -1) ?
																	bb_rw[i] * wk + bb_sd[i] * zi :
					bb_lw[i] * priv_arr->md_z[j*ro_scal->md_dim+m] +bb_rw[i] * wk + bb_sd[i] * zi ;

		}
	}

	return 0;
}

/******************************/
/***** new_trajectoryGPU ******/
/******************************/

void new_trajectoryGPU(
		UINT seq_count,
		UINT ind_flag,
		LoopROScalars* ro_scal,
		struct LoopROArrays*  ro_arr,
		struct LoopPrivateArrays* priv_arr
) {
	UINT m, i, j, k, l;
	UINT dim = ro_scal->md_dim, npathdates = ro_scal->md_nb_path_dates;
	UINT dim_sq = dim*dim, dim_paths = dim*npathdates;

	// generate normal distribution
	mlfi_genmatrix_normalGPU(seq_count, ind_flag, ro_scal, ro_arr, priv_arr);

	// brownian bridge refinement
	// mlfi_brownianbridge_wiener_path(md->bb, CblasNoTrans, md->zd, z, npathdates, dim, npathdates, dim);
	mlfi_brownianbridge_wiener_pathNoTransGPU(ro_scal, ro_arr, priv_arr);

	// expand md->z
	for (i = dim * npathdates - 1; i >= dim; i--) {
		priv_arr->md_z[i] -= priv_arr->md_z[i - dim];
		//priv_arr->md_zd[i] -= priv_arr->md_zd[i - dim];
	}

#if 0
		if(seq_count == 1){ // was == 0
			int i;
			printf("{ ");
			for(i=0; i<ro_scal->sobol_dim; i++) {
				printf("%.6f, ", priv_arr->md_z[i]);
			}
			printf(" }\n");
		}
#endif


	//if(seq_count == 1){ printf("\n\n{ "); }

	// compute trajectory
	for (m = 0; m < ro_scal->inst_num; m++) {
		REAL* c     = ro_arr->inst_coefs + ( ro_scal->inst_c_beg      + m*dim_sq    );
		REAL* vols  = ro_arr->inst_coefs + ( ro_scal->inst_vols_beg   + m*dim_paths );
		REAL* drift = ro_arr->inst_coefs + ( ro_scal->inst_drifts_beg + m*dim_paths );
		REAL* trajRO= ro_arr->inst_coefs + ( ro_scal->inst_trajRO_beg + m*dim       );

		REAL* trajWF= priv_arr->inst_trajWF+ m*dim_paths; // * dim

		for (i = 0; i < npathdates; i++)  {
			for (j = 0; j < dim; j++) {
				REAL temp = 0.0;
				k = dim * i + j;

				for (l = 0; l <= j; l++)
					temp += c[dim * j + l] * priv_arr->md_z[dim * i + l]; //priv_arr->md_zd[dim * i + l];

				//if(seq_count == 1){ printf("%.6f, ", temp); }

				temp = exp(temp * vols[k] + drift[k]);

				trajWF[k] = (k >= dim) ? trajWF[k - dim] * temp : trajRO[ k ] * temp ;

				//RAND_SUM += trajWF[k];
			}
		}
		//if(seq_count == 1){ printf(" }\n\n"); }
#if 0
		if(seq_count == 1){ // was == 1
			int i;
			printf("\n\nTrajectory vector: [ ");
			for(i=0; i<ro_scal->sobol_dim; i++) {
				printf("%.6f, ", trajWF[i]);
			}
			printf(" ]\n\n");
		}
#endif

	}


}


/*****************************/
/***** NOTIFY CASH FLOW ******/
/*****************************/

#include "CosContracts.h"


void trajectoryGPU(
		UINT                      model_num,
		UINT                      cur_iter,
		LoopROScalars*            ro_scal,
		struct LoopROArrays*      ro_arr,
		struct LoopPrivateArrays* priv_arr
//		const mlfi_nmc_model*     model
) {
#if   CONTRACT==1
	trajectory_contract1(model_num, cur_iter, ro_scal, ro_arr, priv_arr); // , model
#elif CONTRACT==2
	trajectory_contract2(model_num, cur_iter, ro_scal, ro_arr, priv_arr);
#elif CONTRACT==3
	trajectory_contract3(model_num, cur_iter, ro_scal, ro_arr, priv_arr);
#else
	trajectory_contract1(model_num, cur_iter, ro_scal, ro_arr, priv_arr);
#endif
}


/*****************************/
/******** MAIN LOOP **********/
/*****************************/


//extern "C"
int parallel_predicate_holds(void* data) {
	driver *md = data;
	int dim = md->dim;
	int npathdates = md->nb_path_dates;

	//printf("Predicate: dim=%d, npathdates=%d, md->bb->bi[0]=%d, sobol_dim=%d, max_poly=%d\n\n\n",
	//		dim, npathdates, md->bb->bi[0], md->rng->sobol->dimension, max_poly);

	return 	(dim!=0 || npathdates != 0) &&
			(md->rng->kind == Sobol)    &&
			//(md->bb->bi[0] == 1)        &&
			//(md->nb_path_dates == 1)    &&   // the last two are probably not necessary for ||
			(md->rng->sobol->dimension <= max_poly);
}

//extern "C"
void setupGPU (
		LoopROScalars*     ro_scal,
		struct LoopROArrays *     ro_arr ,
		struct LoopPrivateArrays* priv_arr
);

/** Cosmin:   run_pricing_index(GPU) contains the main loop to parallelize
 * 	new_trajectory is "void new_trajectory(void *data)"
 * 		from nblackscholes.c
 * 	pc->trajectory is "void trajectory(const mlfi_nmc_model* model)"
 * 		from call.c
 * 		// Original iteration space:
 * 		[(EXsimple_price(EUR, 169.8910648662922, "MC BS (C) 1dim, Sobol, 10000000 paths, 937 ms"))]
 *
 * 		// with changed iteration space for reduction:
 * 		[(EXsimple_price(EUR, 169.8910648662875, "MC BS (C) 1dim, Sobol, 10000000 paths, 975 ms"))]
 *		169.8910648662850 -- OK!!!!
 * 		HASKELL:
 * 		                      169.89106486629217 and 169.89106486628876
 * 		                      169.89106486628876 and 169.89106486628876
 *
 *
 *
 */
//extern "C"
void mainLoopGPU (
		mlfi_nmc_model **models, long nmodels,
		mlfi_nmc_pricing_code *pc,
		void (*new_trajectory)(void *),
		void *data,
		long niters, char *name, long * index
) {

	UINT i,j,k;
	UINT seq_count;

	LoopROScalars            ro_scal;
	struct LoopROArrays      ro_arr;
	struct LoopPrivateArrays priv_arr;

	driver *md = data;
	{ // compute seq_count

		sobol_state_t* s_state = (sobol_state_t *) md->rng->sobol->state;
		seq_count = s_state->sequence_count;
	}


	buildROscalars    (&ro_scal,  niters, data, nmodels, models, pc);
	buildROarrays     (&ro_arr, &ro_scal, data, models);
	buildPrivateArrays(&priv_arr,  &ro_scal, data);

//	printf("Checking: sobol_dim: %d, md_dim: %d, nb_path_dates: %d, rng_nb_dates: %d, rng_dim: %d, inst_num: %d\n\n\n",
//	    	    			ro_scal.sobol_dim, ro_scal.md_dim, ro_scal.md_nb_path_dates,
//	    	    			ro_scal.rng_nb_dates, ro_scal.rng_dim, ro_scal.inst_num);

	{
		setupGPU(&ro_scal, &ro_arr, &priv_arr);

		{ // reset vhat!
			setCHUNK(&ro_scal);
			ro_scal.mc_iter_num = ro_scal.TOT_MC_ITER_NUM;
			UINT sz = ro_scal.inst_num * ro_scal.mc_bigit_num * ro_scal.num_contracts;
			for(i=0; i<sz; i++) {
				priv_arr.model_vhat[i] = 0.0;
			}
		}

	}
    printf("CHUNK: %d\n", ro_scal.CHUNK);

	assert(parallel_predicate_holds(data) && "UNSAFE || ERROR!");
	//printf("Before CPU Cosmin num models: %d %d, seq_count %d \n\n\n", nmodels, niters, seq_count);

	UINT cur_iter = 0;

#ifdef TIME_RESOLUTION_MICROSECOND
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);
#else
        mlfi_timeb  t_start, t_end;
        mlfi_ftime(&t_start);
#endif
        unsigned long int elapsed;

	for( i=0; i < ro_scal.mc_iter_num; i+=ro_scal.CHUNK) {
		UINT ind_flag = 1;

		UINT bound = i+ro_scal.CHUNK;
		bound = (bound > ro_scal.mc_iter_num) ? ro_scal.mc_iter_num : bound;

		for(k = i; k<bound; k++) {
			new_trajectoryGPU(seq_count, ind_flag, &ro_scal, &ro_arr, &priv_arr);
            // assert(seq_count == k && "seq_count invar violated");
			++seq_count;
			//printf("Traj: %.6f ", priv_arr.inst_trajWF[0]);
			for (j = 0; j < nmodels; j++) {
				trajectoryGPU(j, cur_iter, &ro_scal, &ro_arr, &priv_arr); //, models[j]);
			}

			ind_flag = 0;
		}

		cur_iter++;
	}

#ifdef TIME_RESOLUTION_MICROSECOND
        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
#else
        mlfi_ftime(&t_end);
        elapsed = mlfi_diff_time(t_end,t_start);
#endif
    //printf("RAND_SUM: %.16f \n\n\n", RAND_SUM);
    printf("CPU Run Time: %lu !!! \n\n\n", elapsed);
    //printf("Fast branch %llu, slow branches: %llu \n\n\n", fast_branch, slow_branch);

	//seq_count += niters;
	updateState(&priv_arr, seq_count, data);                          // UPDATE IMPLEM!
	updateModels(models, &ro_scal, &priv_arr);
	//updateModels(models, nmodels, cur_iter, (UINT)niters, &priv_arr); // UPDATE IMPLEM!
#if 0
	free(ro_arr.sobol_v_dir);
	free(ro_arr.bb_ia      );
	free(ro_arr.bb_coefs   );
	free(ro_arr.inst_coefs );
	free(ro_arr.pc_coefs   );
	free(priv_arr.sobol_last_num_vec);
	free(priv_arr.md_zd);
	free(priv_arr.md_z );
	free(priv_arr.inst_trajWF);
	free(priv_arr.model_vhat );
	free(priv_arr.vhat_fin);
#endif
	*index = niters;
}


// ToDo:	0. WRITE THE GPU VERSION!!!
// 			1. Test the new way of computing the least-significant bit for sobol.
// 			2. Write multiple-kernel version in which the vector is actually created.
//			3. Implement the reduction in parallel at GPU level, please!


