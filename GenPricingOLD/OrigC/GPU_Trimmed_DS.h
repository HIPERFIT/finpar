/*********************************************************/
/*************** Compact Structures **********************/
/*********************************************************/
#ifndef GPU_TRIMMED_DS_H
#define GPU_TRIMMED_DS_H

#include "Constants.h"


struct LoopROArrays {
	/* sobol direction vector */
	int*    sobol_v_dir;    // [sobol_bit_count][sobol_dim]
	int*    sobol_v_dir_t;  // [sobol_bit_count][sobol_dim]
	UCHAR*  sobol_fix_ind;  // [CHUNK-1]

	/* Brownian Bridge data */
	UINT*  bb_ia;      // [3]x[ro_scal->md_nb_path_dates]
	//UINT*  bb_li;      // [ro_scal->md_nb_path_dates]
	//UINT*  bb_bi;      // [ro_scal->md_nb_path_dates]
	//UINT*  bb_ri;      // [ro_scal->md_nb_path_dates]

	REAL*  bb_coefs;   //[3]x[ro_scal->md_nb_path_dates]
	//REAL*  bb_rw;      // [ro_scal->md_nb_path_dates]
	//REAL*  bb_lw;      // [ro_scal->md_nb_path_dates]
	//REAL*  bb_sd;      // [ro_scal->md_nb_path_dates]

	/* MD instance data */
	REAL*  inst_coefs;  // inst_num*md_dim^2 + 2*inst_num*md_dim*md_nb_path_dates + inst_num*md_dim
	//REAL*  inst_c;      // [inst_num][md_dim^2]
	//REAL*  inst_vols;   // [inst_num][md_dim*md_nb_path_dates]
	//REAL*  inst_drifts; // [inst_num][md_dim*md_nb_path_dates]
	//REAL*  inst_trajRO; // [inst_num][md_dim]

	/** Model Instance attributes **/
	REAL*   pc_coefs;         // inst_num * pc->nb_deterministic_pricers + inst_num*pc->nb_cash_flows
	//REAL*   model_deter_val;  // [inst_num][pc->nb_deterministic_pricers] where pc :: mlfi_nmc_pricing_code*
	//REAL*   model_discounts;  // [inst_num][pc->nb_cash_flows]            where pc :: mlfi_nmc_pricing_code*

};

struct LoopPrivateArrays {
	/* sobol's last_numerator_vector (changeable state) */
	int*   sobol_last_num_vec; // [sobol_dim] == [md_dim*md_nb_path_dates]

	/* md->z and md->zd WF-arrays -- use the same space for both! */
	REAL*  md_zd; // [md_dim*md_nb_path_dates] == [sobol_dim]
	REAL*  md_z;  // [md_dim*md_nb_path_dates]

	/* the trajectory WF data */
	REAL*  inst_trajWF; // [inst_num][md_dim*md_nb_path_dates]

	/** accumulator -- this should be an array **/
	REAL*    model_vhat;  // [inst_num]x[num_contracts]x[num_big_iter]
	double*  vhat_fin;    // [inst_num]x[num_contracts]

	/***** MULTI-KERNEL VERSION *******/
	REAL*  glb_md_zd;       // [mc_iter_num]x[sobol_dim]
	REAL*  glb_md_z;		// [mc_iter_num]x[sobol_dim]
	REAL*  glb_instTrajWF;  // [mc_iter_num]x[inst_num]x[sobol_dim]
	REAL*  glb_vhat;        // [mc_iter_num/CHUNK]x[inst_num]x[num_contracts]
};


#endif





#if 0

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
void mlfi_genmatrix_uniformGPU(mlfi_rng *rng, matrix zd, unsigned int seq_count, int flag_indep) {
	sobol_state_t * s_state;
	unsigned int    dimension, i;
	assert(rng->kind == Sobol && "ERROR INVARIANT rng->kind != Sobol!");

	s_state   = (sobol_state_t *) rng->sobol->state;
	dimension = rng->sobol->dimension;


	// qrng_get(rng->sobol, zd);
	// sobol_getGPU(void * state, unsigned int dimension, double * v, int flag)
	// sobol_getGPU(rng->sobol->state, rng->sobol->dimension, zd, flag);

	if(flag_indep) {
		unsigned int k, gs, gv_k = 0;

		seq_count += 1;
		gs = seq_count >> 1;
		gs = seq_count ^  gs;

		for( i = 0; i < dimension; ++i ) { s_state->last_numerator_vec[i] = 0; }

		for( k = 0; k < SOBOL_BIT_COUNT; ++k ) {
			// compute g_k * v_k
			unsigned int gv_mod = gs % 2;
			gs = gs >> 1;

			if(gv_mod != 0) {
				assert(gv_mod == 1 && "ERROR");

				for( i = 0; i < dimension; ++i ) {
					// xor term g_k * v_k to direction i
					const int direction_i = s_state->v_direction[k][i];
					const int old_numerator_i = s_state->last_numerator_vec[i];
					const int new_numerator_i = old_numerator_i ^ direction_i;
					s_state->last_numerator_vec[i] = new_numerator_i;
					//last_num_vec[i] = last_num_vec[i] ^ s_state->v_direction[k][i]; // gv_mod*
				}
			}
		}

		for(i = 0; i < dimension; ++i) {
			zd[i] = s_state->last_numerator_vec[i] * s_state->last_denominator_inv;
		}

	} else {

		/* Find the position of the least-significant zero in count. */
		int ell = 0;
		int c = s_state->sequence_count;

		while(1) {
			++ell;
			if((c % 2) == 1) c /= 2;
			else break;
		}

		/* Check for exhaustion. */
		//if(ell > SOBOL_BIT_COUNT) failure("bit count exceeded");

		for(i=0; i<dimension; i++) {
			const int direction_i     = s_state->v_direction[ell-1][i];
			const int old_numerator_i = s_state->last_numerator_vec[i];
			const int new_numerator_i = old_numerator_i ^ direction_i; // old_numerator xor direction_i
			s_state->last_numerator_vec[i] = new_numerator_i;
			zd[i] = new_numerator_i * s_state->last_denominator_inv;
		}
	}

	//s_state->sequence_count++;  // do not increment here but later!
}



/* fill a matrix with random normal deviates */
void mlfi_genmatrix_normalGPU(mlfi_rng *rng, matrix zd, UINT seq_count) {
	printf("md->rng->sobol->dim, rng->nb_dates , rng->dim = 1, 1, 1: %d %d %d\n",
			rng->sobol->dimension, rng->nb_dates, rng->dim);

	// Cosmin summary: see above

	mlfi_genmatrix_uniformGPU(
			rng,
			zd,
			((sobol_state_t *)rng->sobol->state)->sequence_count++,
			0
		);
	// Cosmin summary: RW `zd' in [0,rng->nb_dates * rng->dim]
	mlfi_ugaussian_Pinv_vector(zd, rng->nb_dates * rng->dim);
}



///////////////////////////////////////////
////////////// BROWNIAN BRIDGE ////////////
///////////////////////////////////////////

// Cosmin READ-ONLY:   b->sd[0], z[0..cz-1], b->l, rz, cz
//        WRITE-FIRST: w[0..cz-1]
int mlfi_brownianbridge_wiener_pathTransposeGPU(
			const mlfi_brownianbridge *b,
			const matrix z,  const matrix w,
			const int    rz, const int    cz
) {
	size_t m;
	//printf("B->l=1: %d, rz=1: %d\n", b->l, rz);
	if (b->l != rz || b->l <= 0) return -1; /* number of rows doesn't match length of bridge */

	//printf("Cosmin Bridge, (cz, rz)=(1,1): (%d, %d), (b->bi[0]-1) * cz = 0 = %d\n",
	//		cz, rz, ((b->bi[0]-1) * cz));
	for (m=0; m < cz; m++) {
		//w[(b->bi[0]-1) * cz + m] = b->sd[0] * z[m]; /* TODO: eliminate int mul */
		w[m] = b->sd[0] * z[m];
	}
	return 0;
}


/******************************/
/***** new_trajectoryGPU ******/
/******************************/

void new_trajectoryGPU(void *data)
{
	int i, j, k, l, m;

	driver *md = data;
	int dim = md->dim;
	int npathdates = md->nb_path_dates;

	double *z = md->z;

	//if (dim==0 || npathdates == 0) return;
	/* so we can function as a 'deterministic' model */

	/* This is shared by all instances */
	{
		/*   READ-ONLY:  dim = md->rng->sobel->dimension
		 *               md->rng->sobel->state->v_direction[0 : SOBOL_BIT_COUNT-1][1:dim]
		 *               md->rng->sobel->state->last_denominator_inv
		 *   WRITE-FIRST (PRIVATIZED):
		 *               md->zd[0,rng->nb_dates * rng->dim]
		 *               md->rng->sobel->state->last_numerator_vec[1:dim]
		 */
		mlfi_genmatrix_normalGPU(md->rng, md->zd, seq_count);


		// Cosmin READ-ONLY:   md->bb->sd[0], md->zd[0..cz-1], npathdates, dim
		//        WRITE-FIRST: ms->z[0..cz-1]
		mlfi_brownianbridge_wiener_pathNoTransGPU(
				md->bb, md->zd, z, npathdates, dim
			);
	}

	for (i = dim * npathdates - 1; i >= dim; i--) z[i] -= z[i-dim];

	//driver* drv = md;
	//printf("Cosmin model params -- dim: %d, nb_path_dates: %d, nb instances: %d \n",
	//		drv->dim, drv->nb_path_dates, drv->nb_instances);


	/** Cosmin summary:
	 * RO: md->nb_instances, md->instances[0:md->nb_instances-1]
	 *     md->instances[0:md->nb_instances-1]->c[dim*dim-1]
	 *     md->zd[dim*npathdates-1]
	 *     md->instances[0:md->nb_instances-1]->vols[dim*npathdates-1]
	 *     md->instances[0:md->nb_instances-1]->drifts[dim*npathdates-1]
	 * WF: md->instances[0:md->nb_instances-1]->trajectory[dim*npathdates-1]
	 * RO: md->instances[0:md->nb_instances-1]->trajectory[-dim, -1]
	 *
	 */
	//printf("Cosmin (nb_instances, npathdates, dim) == (1,1,1): (%d, %d, %d)\n", md->nb_instances, npathdates, dim);
	/* Instance-specific processing */
	for (m = 0; m < md->nb_instances; m++) {
		instance *inst = md->instances[m];
		for (i = 0; i < npathdates; i++)  {  // hm -- strange here, but it is correct like this!
			for (j = 0; j < dim; j++) { // Cosmin this one is || for trajectory
				double temp = 0.;
				k = dim * i + j;
				for (l = 0; l <= j; l++) temp += inst->c[dim * j + l] * z[dim * i + l];
				temp = exp(temp * inst->vols[k] + inst->drift[k]);

				//printf("Cosmin BUG i %d, j %d, k %d, dim %d, k-dim %d\n", i,j,k,dim,k-dim);
				// COSMIN BUG ... = trajectory[-1] for the first iteration!

				inst->trajectory[k] = inst->trajectory[k - dim] * temp;
				// COSMIN: is in fact inst->trajectory[k] read in call.c via
				//    model->underlyings?
			}
		}
	}
}

static void trajectoryGPU(const mlfi_nmc_model* model) {
	const double underl_val = model->underlyings[0][0];
	const mlfi_contract_index contract_number = 0;
	const double amount = fmax(0., ((underl_val - 4000.) * model->deterministic_values[0]));
	const mlfi_pay_date_index date_index = 0;

	instance *inst = (instance*) model->model_data;
	inst->vhat[contract_number] += amount * inst->discounts[date_index];

	return;

}

void mainLoopGPU (
		mlfi_nmc_model **models, long nmodels,
		mlfi_nmc_pricing_code *pc,
		void (*new_trajectory)(void *),
		void *data,
		long niters, char *name, long * index ) {

	UINT i,j,k;

	assert(parallel_predicate_holds(data) && "UNSAFE || ERROR!");
	//printf("Cosmin num models: %d %d\n\n\n", nmodels, niters);


	for (i = 0; i < niters; i++) { // make it blocked!

		//new_trajectoryGPU(data, i);
		new_trajectoryGPU(data);

		for (j = 0; j < nmodels; j++) {
			trajectoryGPU(models[j]);
		}
	}

	*index = niters;
}


#endif


#if 0
static void trajectoryGPU_GEN1(
		UINT               model_num,
		UINT               cur_iter,
		LoopROScalars*     ro_scal,
		struct LoopPrivateArrays* priv_arr,
		const mlfi_nmc_model*     model
) {
	//const double underl_val = model->underlyings[0][0];
	const REAL underl_val = priv_arr->inst_trajWF[model_num];

	const REAL amount = fmax(0., ((underl_val - 4000.) * ro_scal->model_deter_val0));
	//const double amount = fmax(0., ((underl_val - 4000.) * model->deterministic_values[0]));

	//const mlfi_contract_index contract_number = 0;
	//const mlfi_pay_date_index date_index = 0;
	//instance *inst = (instance*) model->model_data;
	//inst->vhat[contract_number] += amount * inst->discounts[date_index];

	// more parameters are needed here to implement access to vhat accumulator
	//      but we keep it simple for testing purposes!
	priv_arr->model_vhat[cur_iter] += amount * ro_scal->model_discounts0;
	//inst->vhat[0] += amount * ro_scal->model_discounts0;
	return;

}
#endif
