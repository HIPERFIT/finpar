#include "CosminGPU.h"
#include "TimeHelperChr.h"

//int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1);

/**********************************/
/******* HELPER FUNCTIONS  ********/
/**********************************/

UINT logNextPow2(UINT n) {
	UINT ret    = 1;
	UINT logret = 0;
	while(n>ret) {
		ret     = ret << 1;
		logret += 1;
	}
	return logret;
}

UINT prevMultipleOf(UINT n, UINT mul) {
	if(n % mul == 0) return n;
	UINT res = 0;
	while(res < n) res += mul;
	return res - mul;
}

static UINT getWorkSize(UINT num_iter, UINT CHUNK, UINT BLOCK) {
	UINT tmp_rem = num_iter % CHUNK;
	UINT res = (tmp_rem == 0) ? (num_iter / CHUNK) : ((num_iter / CHUNK) + 1);

	//printf("IN GETWORKSIZE: NUM_ITER: %d, part num_iter/CHUNK: %d\n\n", num_iter, res);

	if(res % BLOCK == 0) return res;

	res = (res / BLOCK) * BLOCK + BLOCK;
	return res;
}


#define logMAX_CHUNK 8
enum GPU_KERNEL { PRIV, VECT };

GPU_KERNEL priv_or_vect_kernel(LoopROScalars* ro_scal, cl_device_id device) {
#if   (_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV == 1)
	return VECT;
#elif (_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV == 2)
	return PRIV;
#else
	cl_ulong mem_size, THRESHOLD;
	clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
	THRESHOLD = (mem_size >> 10) / 3;
	bool is_priv = (ro_scal->sobol_dim+ro_scal->num_contracts <= THRESHOLD);
	printf("\n\nCOST MODEL THRESHOLD: %u, PRIV_KERNEL: %d\n", (UINT)THRESHOLD, (int)is_priv);
	return (is_priv) ? PRIV : VECT;
#endif
}

GPU_KERNEL discriminate_cost_model(LoopROScalars* ro_scal, cl_device_id device) {
	bool is_priv;
	UINT THRESHOLD;

	{
		cl_ulong local_mem_size;
		clGetDeviceInfo(device,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(local_mem_size),&local_mem_size,NULL);
#ifdef _OPTIMIZATION_TILE
                ro_scal->TILE_FACT = TILE;
#else
		ro_scal->TILE_FACT = local_mem_size/(256*3*sizeof(REAL));
#endif
		printf("LOC MEM SIZE: %lu, TILE_FACT: %u, TILE: %u \n\n\n", local_mem_size, ro_scal->TILE_FACT, TILE);

		assert(ro_scal->TILE_FACT == TILE && "SET PROPER TILE in Constants.h");
	}

	ro_scal->logBLOCK = 7;
	ro_scal->BLOCK = 1 << ro_scal->logBLOCK;

	is_priv = (priv_or_vect_kernel(ro_scal, device) == PRIV);

    if(!is_priv) {
    	// get device memory size info
        cl_ulong glob_mem_size;
        clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(glob_mem_size),&glob_mem_size,NULL);
        //printf("GLOB MEM SIZE: %lu\n\n\n", glob_mem_size);

    	// check if there is enough memory
    	UINT perit_mem_size = ro_scal->sobol_dim*(ro_scal->inst_num+1) + ro_scal->num_contracts;
    	perit_mem_size *= sizeof(REAL);

    	UINT max_mc_iter_num = ( (glob_mem_size * 3)/4 ) / perit_mem_size;

    	UINT BLOCK = ro_scal->BLOCK;
    	if(max_mc_iter_num % BLOCK != 0)
    		max_mc_iter_num = (max_mc_iter_num / BLOCK)*BLOCK + BLOCK;

    	max_mc_iter_num = prevMultipleOf(max_mc_iter_num, ro_scal->CHUNK);

    	ro_scal->mc_iter_num = MIN(max_mc_iter_num, ro_scal->mc_iter_num);
    }

	{
		cl_uint compute_units, tot_cuda_cores;
		clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
							sizeof(compute_units), &compute_units, NULL);
		cl_uint comp_capabil_major, comp_capabil_minor;
		clGetDeviceInfo(device, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
							sizeof(cl_uint), &comp_capabil_major, NULL);
		clGetDeviceInfo(device, CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,
							sizeof(cl_uint), &comp_capabil_minor, NULL);
		tot_cuda_cores = compute_units *
				ConvertSMVer2CoresCopy(comp_capabil_major, comp_capabil_minor);

		//printf("Compute units: %d, tot_core_num: %d\n", compute_units, tot_cuda_cores);

		UINT div = (1 << logNextPow2(tot_cuda_cores)) * 64;
		ro_scal->logCHUNK = logNextPow2(ro_scal->mc_iter_num / div );  // was div*2
		ro_scal->logCHUNK = MIN(ro_scal->logCHUNK, logMAX_CHUNK);
		if(ro_scal->logCHUNK <= 3) ro_scal->logCHUNK++;

		// TO DO: move here the computation of the FIX_IND ARRAY!
	}

	ro_scal->CHUNK = 1 << ro_scal->logCHUNK;

	printf("\n\nCOST MODEL COMPUTED: MC_ITERS: %u, MC_TOT_ITERS: %u, CHUNK: %u, BLOCK: %u\n\n",
			ro_scal->mc_iter_num, ro_scal->TOT_MC_ITER_NUM, ro_scal->CHUNK, ro_scal->BLOCK);

	return (is_priv) ? PRIV : VECT;
}

void computeSobolFixIndex(LoopROArrays * ro_arr, UINT CHUNK) {
	// Given CHUNK, the most-significant zero of iterations 1 ... CHUNK-1 is fixed.
	if(ro_arr->sobol_fix_ind == NULL)
		ro_arr->sobol_fix_ind = (UCHAR*)malloc(sizeof(UCHAR)*(1<<logMAX_CHUNK));
	for(int k=1; k<CHUNK-1; k++) {
		UINT gs  = k;
		UINT ell = 0;
		while(gs & 1) {
			ell++;
		    gs >>= 1;
		}
		ro_arr->sobol_fix_ind[k] = ell;
		//printf("SOB_FI[%d]=%d,   ", k, ell);
	}
}

/*****************************/
/***** PRIVATIZED KERNEL *****/
/*****************************/

#define CREATE_AND_ENQUEUE_CL_BUFFER(X) \
		ocl_arr.X = clCreateBuffer( cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, \
															cur_size, ro_arr->X, &ciErr2 ); \
		ciErr1 |= ciErr2;


static void oclAllocArrays_PrivKernel (
		oclLoopArrays&     ocl_arr,
		cl_context         cxGPUContext,
		cl_command_queue&  cqCommandQueue,
		cl_int&            ciErr1,
		LoopROScalars*     ro_scal,
		LoopROArrays *     ro_arr,
		LoopPrivateArrays* priv_arr
) {
	cl_int ciErr2;
	unsigned long int cur_size;
	/*********************/
	/*** 1. RO scalars ***/
	/*********************/
	{
		cur_size = sizeof(LoopROScalars);
		ocl_arr.ro_scals = clCreateBuffer(
					cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
				);
		ciErr1 |= ciErr2;
		ciErr1 |= clEnqueueWriteBuffer(cqCommandQueue, ocl_arr.ro_scals, CL_TRUE, 0,
					cur_size, ro_scal, 0, NULL, NULL);
	}

	/***************************/
	/*** 2. SOBOL DIR VECTOR ***/
	/***************************/
	cur_size = (UINT)ro_scal->sobol_bit_count*ro_scal->sobol_dim*sizeof(UINT);
	CREATE_AND_ENQUEUE_CL_BUFFER(sobol_v_dir); // sobol_v_dir

	/**********************************/
	/*** 3. BROWNIAN BRIDGE Data RO ***/
	/**********************************/
	{
		cur_size = 3*ro_scal->md_nb_path_dates*sizeof(UINT);
		CREATE_AND_ENQUEUE_CL_BUFFER(bb_ia);

		cur_size = 3*ro_scal->md_nb_path_dates*sizeof(REAL);
		CREATE_AND_ENQUEUE_CL_BUFFER(bb_coefs);
	}

	/******************************/
	/*** 4. MD INSTANCE Data RO ***/
	/******************************/
	{
		cur_size = (ro_scal->inst_trajRO_beg + ro_scal->inst_num*ro_scal->md_dim) * sizeof(REAL);
		CREATE_AND_ENQUEUE_CL_BUFFER(inst_coefs);
	}

	/******************************/
	/*** 5. MODEL ATTRIBUTES RO ***/
	/******************************/
	{
		cur_size = (ro_scal->pc_discounts_beg + ro_scal->inst_num*ro_scal->num_cash_flows) * sizeof(REAL);
		CREATE_AND_ENQUEUE_CL_BUFFER(pc_coefs);

		{
			ocl_arr.sobol_fix_ind = clCreateBuffer( cxGPUContext, CL_MEM_READ_ONLY |
							CL_MEM_COPY_HOST_PTR, ro_scal->CHUNK, ro_arr->sobol_fix_ind, &ciErr2 );
			ciErr1 |= ciErr2;
		}
	}


	/*******************************/
	/*** FINALLY, GLOBAL VHAT WO ***/
	/*******************************/
	{
		//red_sz = getGlobalRedArrSize( ro_scal->mc_iter_num, ro_scal->CHUNK*ro_scal->BLOCK );
		cur_size = getWorkSize(ro_scal->mc_iter_num, ro_scal->CHUNK, ro_scal->BLOCK) / ro_scal->BLOCK;
		cur_size *= ro_scal->inst_num*ro_scal->num_contracts*sizeof(REAL);
		priv_arr->glb_vhat = (REAL*)malloc(cur_size);
		ocl_arr.model_vhat = clCreateBuffer(
				cxGPUContext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
				cur_size, priv_arr->glb_vhat, &ciErr2
			);
		ciErr1 |= ciErr2;
	}

	// 2. the rest: priv_arr->sobol_last_num_vec/md_zd/inst_trajWF
	//      are implemented as __local (private) hence do not require
	//      global memory allocation.

	oclCheckError(ciErr1, CL_SUCCESS);
}


void runGPU_PRIV(
		LoopROScalars*     	ro_scal,
		LoopPrivateArrays* 	priv_arr,
		oclLoopArrays&		dev_arr,
		cl_command_queue& 	cqCommandQueue,
		cl_program&       	cpProgram,
	    size_t*				globalWorkSize,
	    size_t*				localWorkSize
) {
	cl_kernel ckGenPricing = NULL; // OpenCL kernel
	cl_int    ciErr1, ciErr2;

	int  counter = 0, priv_sz = 0; // ro_scals
	const int PRIV_MULT   = sizeof(REAL)*ro_scal->BLOCK;


	ckGenPricing = clCreateKernel(cpProgram, "payoffGPU", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


	// 1. RO SCALARS struct  //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals   );

	// 2. RO SOBOL ARR //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_v_dir);

	// 3. RO BROWNIAN BRIDGE //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_ia);
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_coefs);


	// 4. RO MD INSTANCE DATA //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.inst_coefs);

	// 5. MODEL ATTRIBUTES RO //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.pc_coefs);

	// 6. fix indexes to optimize sobol rand num gen. //
	ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_fix_ind);

    ////////////////////////
    // 7. GLOBAL WO Array //
    ////////////////////////
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.model_vhat );

    // 8. LOCAL (PRIVATE ARRAYS) MD_Z
    priv_sz = PRIV_MULT*ro_scal->sobol_dim;
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, priv_sz, NULL);

    // 9. LOCAL (PRIVATE ARRAY) VHAT
    priv_sz = PRIV_MULT * ro_scal->inst_num * ro_scal->num_contracts;
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, priv_sz, NULL); // model_vhat (the private one)

    oclCheckError(ciErr1, CL_SUCCESS);

    // 10. ENQUEUE KERNEL!!!  //
    ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing, 1, NULL,
    								globalWorkSize, localWorkSize, 0, NULL, NULL);
#if 0
    printf("AFTER ENQUEUEING KERNEL code: %d %d %lu %lu %lu %lu %lu\n\n",
    		ciErr1,  CL_INVALID_WORK_GROUP_SIZE,
    		globalWorkSize[0], localWorkSize[0],
    		sizeof(size_t), sizeof(UINT), globalWorkSize[0] % localWorkSize[0]);
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    // 11. FINALLY, WRITE BACK!!! //
    priv_sz = (globalWorkSize[0]/ro_scal->BLOCK)*ro_scal->inst_num*ro_scal->num_contracts*sizeof(REAL);
    ciErr1 |= clEnqueueReadBuffer(cqCommandQueue, dev_arr.model_vhat, CL_TRUE, 0,
    		priv_sz, priv_arr->glb_vhat, 0, NULL, NULL);
    oclCheckError(ciErr1, CL_SUCCESS);

    clReleaseKernel(ckGenPricing);
}

void reduceVHAT_PrivCPU(LoopROScalars* ro_scal, REAL* vhat, double* cpu_vhat) {
	// final reduction on CPU!!!
	UINT SZ = getWorkSize(ro_scal->mc_iter_num, ro_scal->CHUNK, ro_scal->BLOCK) / ro_scal->BLOCK;

	for(int ii = 0; ii<ro_scal->inst_num; ii++) {
		for(int jj = 0; jj<ro_scal->num_contracts; jj++) {
			double final_res = 0.0;

			// this forms an inner dimension!
			for( UINT k=0; k<SZ; k++ ) {
				final_res += vhat[k];
			}
			vhat += SZ;
			cpu_vhat[ii*ro_scal->num_contracts + jj] = final_res;

	    	printf("FINAL RESULT GPU for(inst,contract_num): (%d,%d) IS: %.16f !",
	    			ii, jj, final_res/ro_scal->mc_iter_num);
		}
	}
}

/*****************************/
/***** VECTORIZED KERNEL *****/
/*****************************/

static void oclAllocArrays_VectKernel (
		oclLoopArrays&     ocl_arr,
		cl_context         cxGPUContext,
		cl_command_queue&  cqCommandQueue,
		cl_int&            ciErr1,
		LoopROScalars*     ro_scal,
		LoopROArrays *     ro_arr,
		LoopPrivateArrays* priv_arr
) {
	const UINT CHUNK = ro_scal->CHUNK;
	const UINT BLOCK = ro_scal->BLOCK;
	cl_int ciErr2;
	unsigned long int cur_size;
	/*********************/
	/*** 1. RO scalars ***/
	/*********************/
	cur_size = sizeof(LoopROScalars);
	ocl_arr.ro_scals = clCreateBuffer(
				cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
			);
	ciErr1 |= ciErr2;
	ciErr1 |= clEnqueueWriteBuffer(cqCommandQueue, ocl_arr.ro_scals, CL_TRUE, 0,
				cur_size, ro_scal, 0, NULL, NULL);

	/***************************/
	/*** 2. SOBOL DIR VECTOR ***/
	/***************************/
	cur_size = ro_scal->sobol_bit_count*ro_scal->sobol_dim*sizeof(UINT);
	ocl_arr.sobol_v_dir   = clCreateBuffer( cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
													cur_size, ro_arr->sobol_v_dir_t, &ciErr2 );
	ocl_arr.sobol_fix_ind = clCreateBuffer( cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
													CHUNK, ro_arr->sobol_fix_ind, &ciErr2 );
	ciErr1 |= ciErr2;

	/**********************************/
	/*** 3. BROWNIAN BRIDGE Data RO ***/
	/**********************************/
	cur_size = 3*ro_scal->md_nb_path_dates*sizeof(UINT);
	CREATE_AND_ENQUEUE_CL_BUFFER(bb_ia);

	cur_size = 3*ro_scal->md_nb_path_dates*sizeof(REAL);
	CREATE_AND_ENQUEUE_CL_BUFFER(bb_coefs);

	/******************************/
	/*** 4. MD INSTANCE Data RO ***/
	/******************************/
	cur_size = (ro_scal->inst_trajRO_beg + ro_scal->inst_num*ro_scal->md_dim) * sizeof(REAL);
	CREATE_AND_ENQUEUE_CL_BUFFER(inst_coefs);

	/******************************/
	/*** 5. MODEL ATTRIBUTES RO ***/
	/******************************/
	cur_size = (ro_scal->pc_discounts_beg + ro_scal->inst_num*ro_scal->num_cash_flows) * sizeof(REAL);
	CREATE_AND_ENQUEUE_CL_BUFFER(pc_coefs);  // pc_coefs

	/*******************************/
	/*** 6. PRIV: MD_Z and MD_ZD ***/
	/*******************************/
	cur_size = getWorkSize(ro_scal->mc_iter_num, CHUNK, BLOCK) *
					ro_scal->sobol_dim*CHUNK*sizeof(REAL);

#if 1
	ocl_arr.md_zd = clCreateBuffer(
			cxGPUContext, CL_MEM_READ_WRITE,
			cur_size*ro_scal->inst_num, NULL, &ciErr2
		);
#else
	priv_arr->glb_md_zd = (REAL*)malloc(cur_size);
	ocl_arr.md_zd = clCreateBuffer(
			cxGPUContext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
			cur_size*ro_scal->inst_num, priv_arr->glb_md_zd, &ciErr2
		);
#endif
	ciErr1 |= ciErr2;

	ocl_arr.md_z = clCreateBuffer(
			cxGPUContext, CL_MEM_READ_WRITE, //CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
			cur_size, NULL, &ciErr2
		);
	ciErr1 |= ciErr2;

	/*******************************/
	/*** 7. PRIV: INST_TRAJ      ***/
	/*******************************/
	ocl_arr.instTrajWF = ocl_arr.md_zd; // just reuse md_zd

	oclCheckError(ciErr1, CL_SUCCESS);

	/*******************************/
	/*** FINALLY, GLOBAL VHAT WO ***/
	/*******************************/
	cur_size = getWorkSize( ro_scal->mc_iter_num, CHUNK, BLOCK );
	cur_size*= ro_scal->inst_num*ro_scal->num_contracts*sizeof(REAL);
#if 0
	printf("GLB_VHAT SIZE: %u, %u \n\n",
			(UINT)(cur_size/sizeof(REAL)),
			(UINT)(ro_scal->mc_iter_num*ro_scal->inst_num*ro_scal->num_contracts/CHUNK)
		);
#endif
	priv_arr->glb_vhat = (REAL*)malloc( cur_size );
	//assert(priv_arr->glb_vhat && "MALLOC FAILED!!!");
	ocl_arr.model_vhat = clCreateBuffer(
			cxGPUContext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
			cur_size, priv_arr->glb_vhat, &ciErr2
		);
	ciErr1 |= ciErr2;

	oclCheckError(ciErr1, CL_SUCCESS);
}


void runGPU_VECT(
		LoopROScalars*     	ro_scal,
		LoopPrivateArrays* 	priv_arr,
		oclLoopArrays&		dev_arr,
		cl_command_queue& 	cqCommandQueue,
		cl_program&       	cpProgram,
	    size_t*				globalWorkSize,
	    size_t*				localWorkSize
) {
	cl_kernel ckGenPricing_sobol = NULL, ckGenPricing_invg   = NULL, ckGenPricing_brown = NULL;
	cl_kernel ckGenPricing_traj  = NULL, ckGenPricing_finred = NULL;
	cl_int    ciErr1, ciErr2;

	int  counter = 0;


	{ // SOBOL KERNEL!
		counter = 0;
		ckGenPricing_sobol = clCreateKernel(cpProgram, "mlfi_genmatrix_uniform2", &ciErr1);
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals   );
		ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_fix_ind);

		ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_v_dir);
		ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd );

        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing_sobol, 1, NULL,
        								globalWorkSize, localWorkSize, 0, NULL, NULL);
        oclCheckError(ciErr1, CL_SUCCESS);
	}

	{ // INV GAUSSIAN KERNEL!
		//shrLog("Call Inv Gausssian kernel on GPU...\n\n");
		counter = 0;
		ckGenPricing_invg = clCreateKernel(cpProgram, "mlfi_ugaussian_Pinv_vector1", &ciErr1);
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clSetKernelArg(ckGenPricing_invg, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd );

		ciErr1 |= clSetKernelArg(ckGenPricing_invg, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);

        //UINT priv_sz = (sizeof(REAL) + sizeof(char)) * TILE * BLOCK;
		UINT priv_sz = sizeof(REAL) * ro_scal->TILE_FACT * ro_scal->BLOCK;
        ciErr1 |= clSetKernelArg(ckGenPricing_invg, counter++, priv_sz, NULL); // local_space!
        oclCheckError(ciErr1, CL_SUCCESS);

		const unsigned long int inv_ker_sz = globalWorkSize[0]; //* CHUNK * ro_scal->sobol_dim;
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing_invg, 1, NULL,
        		&inv_ker_sz, localWorkSize, 0, NULL, NULL);
        oclCheckError(ciErr1, CL_SUCCESS);

	}

	{ // BROWNIAN BRIDGE KERNEL!
		//shrLog("Call Brownian Bridge kernel on GPU...\n\n");
		counter = 0;
		ckGenPricing_brown = clCreateKernel(cpProgram, "mlfi_brownianbridge_wiener_path1", &ciErr1);
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);

		ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_ia);
		ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_coefs);

		ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd );
		oclCheckError(ciErr1, CL_SUCCESS);
		ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.md_z );
		oclCheckError(ciErr1, CL_SUCCESS);

        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing_brown, 1, NULL,
        								globalWorkSize, localWorkSize, 0, NULL, NULL);
        oclCheckError(ciErr1, CL_SUCCESS);
	}


	{ // compute trajectory
		//shrLog("Call Compute Trajectory kernel on GPU...\n\n");
		bool optim_kernel = (2*ro_scal->md_dim <= TILE);
		counter = 0;
		ckGenPricing_traj = (optim_kernel) ?
				clCreateKernel(cpProgram, "mlfi_comp_traj1",      &ciErr1) :
				clCreateKernel(cpProgram, "mlfi_comp_traj_unopt", &ciErr1) ;
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);

		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.inst_coefs);
		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.pc_coefs);

		UINT priv_sz = sizeof(REAL) * ro_scal->TILE_FACT * ro_scal->BLOCK;
		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, priv_sz, NULL); // local_space!

		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.md_z );
		ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.instTrajWF );

		//ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.model_vhat );
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing_traj, 1, NULL,
					globalWorkSize, localWorkSize, 0, NULL, NULL);
		oclCheckError(ciErr1, CL_SUCCESS);
	}

	{ // final reduction!
		//shrLog("Call Compute Reduction kernel on GPU...\n\n");
		counter = 0;
		ckGenPricing_finred = clCreateKernel(cpProgram, "mlfi_reduction_step1", &ciErr1);
		oclCheckError(ciErr1, CL_SUCCESS);

		ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);


		ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.inst_coefs);
		ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.pc_coefs);

		ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.instTrajWF );
		ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.model_vhat );

		ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing_finred, 1, NULL,
					globalWorkSize, localWorkSize, 0, NULL, NULL);
		oclCheckError(ciErr1, CL_SUCCESS);

        ciErr1 |= clEnqueueReadBuffer(cqCommandQueue, dev_arr.model_vhat, CL_TRUE, 0,
        		globalWorkSize[0]*ro_scal->num_contracts*ro_scal->inst_num*sizeof(REAL),
        		priv_arr->glb_vhat, 0, NULL, NULL);
        oclCheckError(ciErr1, CL_SUCCESS);
	}

    clReleaseKernel(ckGenPricing_sobol );
    clReleaseKernel(ckGenPricing_invg  );
    clReleaseKernel(ckGenPricing_brown );
    clReleaseKernel(ckGenPricing_traj  );
    clReleaseKernel(ckGenPricing_finred);
}

void reduceVHAT_VectCPU(LoopROScalars* ro_scal, REAL* vhat, double* cpu_vhat) {
	UINT WKSZ = getWorkSize(ro_scal->mc_iter_num, ro_scal->CHUNK, ro_scal->BLOCK)/ro_scal->BLOCK;
	UINT step = ro_scal->BLOCK * ro_scal->inst_num * ro_scal->num_contracts;

	for(int b = 0; b < WKSZ; b++) {
		for(int i=0; i<ro_scal->inst_num*ro_scal->num_contracts; i+=ro_scal->num_contracts) {
			for(int j=0; j<ro_scal->num_contracts; j++) {
				cpu_vhat[i+j] += vhat[(i+j)<<logWARP];
			}
		}
		vhat += step;
	}
}

/***********************************************/
/***** BUILD CONTEXT/COMMAND QUEUE/PROGRAM *****/
/***********************************************/

GPU_KERNEL build_for_GPU(
		LoopROScalars*   	ro_scal,
		cl_context& 		cxGPUContext,
		cl_command_queue& 	cqCommandQueue,
		cl_uint& 			nDevice,
		cl_device_id*&     	cdDevices,
		cl_program&       	cpProgram
) {
	cl_platform_id 	cpPlatform;  // OpenCL platform
	cl_int 			ciErr1, ciErr2;

    // 2. get the platforms' ids
    shrLog("Get platforms...\n");
    ciErr1 = oclGetPlatformID(&cpPlatform);
    oclCheckError(ciErr1, CL_SUCCESS);

    // 3. get devices
    shrLog("Get devices...\n");
    ciErr1 = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 0, NULL, &nDevice);
    oclCheckError(ciErr1, CL_SUCCESS);
    cdDevices = (cl_device_id *)malloc(nDevice * sizeof(cl_device_id) );
    ciErr1 = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, nDevice, cdDevices, NULL);
    oclCheckError(ciErr1, CL_SUCCESS);

    // 4. create context
    shrLog("Create context...\n");
    cxGPUContext = clCreateContext(0, nDevice, cdDevices, NULL, NULL, &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);

    // 5. create command queues for all available devices
    shrLog("clCreateCommandQueue\n");

    //for (cl_uint i = 0; i < nDevice; i++) {
    	cqCommandQueue = clCreateCommandQueue(cxGPUContext, cdDevices[dev_id], 0, &ciErr1);
    	oclCheckErrorEX(ciErr1, CL_SUCCESS, NULL);
    //}

    for (cl_uint i = 0; i < nDevice; i++) {
#if (!OpenCL_ONLY)
    	oclPrintDevInfo(LOGBOTH, cdDevices[i]);
#else
    	oclPrintDevInfo(0, cdDevices[i]);
#endif
	}

    shrLog("\nUsing %d GPU(s)...\n\n", nDevice);

    // 6. build program
    GPU_KERNEL kernel_type = priv_or_vect_kernel(ro_scal, cdDevices[dev_id]);
    const char *clSourcefile = (kernel_type == PRIV) ?
        		"payoffPrivGPU.cl" : "payoffVectGPU.cl"; // kernel file

#if (!OpenCL_ONLY)
    shrSetLogFileName ("oclGenericPrincing.txt");
#endif

    shrLog("Create and build program from %s...\n", clSourcefile);
    size_t szKernelLength; // Byte size of kernel code
    char *cSourcePath = (char*)clSourcefile; //shrFindFilePath(clSourcefile, argv[0]);

#if (!OpenCL_ONLY)
    shrCheckError(cSourcePath != NULL, shrTRUE);
#else
    assert(cSourcePath != NULL && "UNDEFINED cSourcePath");
#endif

    char * kernel_code = oclLoadProgSource(cSourcePath, "// My comment\n", &szKernelLength);

    oclCheckError(kernel_code != NULL, shrTRUE);
    cpProgram = clCreateProgramWithSource(cxGPUContext, 1, (const char **)&kernel_code, &szKernelLength, &ciErr1);
    //printf("Before building error: %d, code: %d, success: %d, num_dev: %d\n\n",
    //		ciErr1, CL_BUILD_PROGRAM_FAILURE, CL_SUCCESS, nDevice);
    //
    /*
    // valid options: http://www.khronos.org/registry/cl/sdk/1.2/docs/man/xhtml/clBuildProgram.html
    // -cl-fast-relaxed-math: enables all optimizations -> 8-10% speed increase in contract 1, but segmentation fault in contract 2
    // -cl-unsafe-math-optimizations -> 8-10% speed increase in contract 1, but segmentation fault in contract 2
    // -cl-no-signed-zeros, -cl-strict-aliasing, -cl-mad-enable -> segmentation fault in contract 2
    const char *options = "";
    ciErr1 |= clBuildProgram(
      cpProgram, // program
      0,         // num_devices
      NULL,      // device list
      options,   // options
      NULL,      // callback
      NULL       // user data
    );
    */
    ciErr1 |= clBuildProgram(cpProgram, 0, NULL, NULL, NULL, NULL);

//oclLogBuildInfo(cpProgram, oclGetFirstDev(cxGPUContext));
//oclLogPtx(cpProgram, oclGetFirstDev(cxGPUContext), "payoffGPU.ptx"/*"GenericPricing.ptx"*/);

    // 7. check errors!
    if (ciErr1 != CL_SUCCESS) {
        // write out standard error, Build Log and PTX, then cleanup and exit
#if (!OpenCL_ONLY)
        shrLogEx(LOGBOTH | ERRORMSG, (double)ciErr1, STDERROR);
#endif
        printf("BUILDING COSMIN ERROR: %d", ciErr1);
        oclLogBuildInfo(cpProgram, oclGetFirstDev(cxGPUContext));
        oclLogPtx(cpProgram, oclGetFirstDev(cxGPUContext), "payoffGPU.ptx"/*"GenericPricing.ptx"*/);

        if (ciErr1 == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(cpProgram, cdDevices[dev_id], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(cpProgram, cdDevices[dev_id], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            printf("%s\n", log);
        }

        oclCheckError(ciErr1, CL_SUCCESS);
    }

    return kernel_type;
}

extern "C" void setupGPU (
		LoopROScalars*     ro_scal,
		LoopROArrays *     ro_arr ,
		LoopPrivateArrays* priv_arr
) {
    cl_context cxGPUContext;                        // OpenCL context
    cl_command_queue cqCommandQueue[16]; 			// OpenCL command que
    cl_uint nDevice;                                // OpenCL device count
    cl_device_id* cdDevices;                        // OpenCL device list
    cl_program cpProgram;                           // OpenCL program
    cl_kernel  ckGenPricing = NULL;                 // OpenCL kernel
    cl_int ciErr1, ciErr2;                          // Error code var

    size_t globalWorkSize[1]; 	// = {MT_RNG_COUNT};      // 1D var for Total # of work items
    size_t localWorkSize [1]; 	//  = {128};               // 1D var for # of work items in the work group
    size_t global_arr_sz; // = getGlobalRedArrSize( numiter, CHUNK*BLOCK );



    GPU_KERNEL kernel_type = build_for_GPU(
    		ro_scal, cxGPUContext, cqCommandQueue[dev_id],
    		nDevice, cdDevices, cpProgram
    	);

    // 0. Allocate memory
    oclLoopArrays dev_arr;

#ifdef TIME_RESOLUTION_MICROSECOND
    struct timeval t_start, t_end, t_diff;
    gettimeofday(&t_start, NULL);
#else
    mlfi_timeb  t_start, t_end;
    mlfi_ftime(&t_start);
#endif
    unsigned long int elapsed;

	{
    	GPU_KERNEL kernel_type1 = discriminate_cost_model(ro_scal, cdDevices[dev_id]);
    	assert(kernel_type1 == kernel_type && "INCONSISTANT KERNEL TYPES!");
    	computeSobolFixIndex(ro_arr, ro_scal->CHUNK);

    	globalWorkSize[0]    = getWorkSize( ro_scal->mc_iter_num, ro_scal->CHUNK, ro_scal->BLOCK );
    	localWorkSize [0]    = ro_scal->BLOCK;
	}

    if(kernel_type == PRIV) { // CALL THE PRIVATE KERNEL

    	oclAllocArrays_PrivKernel(
    			dev_arr, cxGPUContext, cqCommandQueue[dev_id],
    			ciErr1, ro_scal, ro_arr, priv_arr
    		);

        runGPU_PRIV (
        		ro_scal, priv_arr, dev_arr, cqCommandQueue[dev_id], cpProgram,
        		globalWorkSize, localWorkSize
        	);

#ifdef TIME_RESOLUTION_MICROSECOND
        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
#else
        mlfi_ftime(&t_end);
        elapsed = mlfi_diff_time(t_end,t_start);
#endif

        printf("\n\nGPU Run Time: %lu !\n\n", elapsed);

        reduceVHAT_PrivCPU(ro_scal, priv_arr->glb_vhat, priv_arr->vhat_fin);

    } else { // VECT VERSION!

    	oclAllocArrays_VectKernel(
    			dev_arr, cxGPUContext, cqCommandQueue[dev_id],
    			ciErr1, ro_scal, ro_arr, priv_arr
    		);

        UINT UB = ro_scal->mc_iter_num, cur_iter = 0, sob_ini_count = ro_scal->sobol_count_ini;

    	for(UINT cur_iter = 0; cur_iter < ro_scal->TOT_MC_ITER_NUM; cur_iter+=ro_scal->mc_iter_num) {
    		if(cur_iter + ro_scal->mc_iter_num > ro_scal->TOT_MC_ITER_NUM) {
    			ro_scal->mc_iter_num = ro_scal->TOT_MC_ITER_NUM - cur_iter;
    			globalWorkSize[0] = getWorkSize( ro_scal->mc_iter_num,ro_scal->CHUNK,ro_scal->BLOCK);
    		}
    		if(cur_iter != 0) {
    			ro_scal->sobol_count_ini = cur_iter + sob_ini_count;

    			clReleaseMemObject(dev_arr.ro_scals);
    			unsigned long cur_size = sizeof(LoopROScalars);
    			cl_int ciErr2;
    			dev_arr.ro_scals = clCreateBuffer(
    						cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
    					);
    			ciErr2 = clEnqueueWriteBuffer(cqCommandQueue[dev_id], dev_arr.ro_scals, CL_TRUE, 0,
    						cur_size, ro_scal, 0, NULL, NULL);
    			oclCheckError(ciErr2, CL_SUCCESS);
    		}

    		runGPU_VECT (
        		ro_scal, priv_arr, dev_arr, cqCommandQueue[dev_id], cpProgram,
        		globalWorkSize, localWorkSize
        	);

    		reduceVHAT_VectCPU(ro_scal, priv_arr->glb_vhat, priv_arr->vhat_fin);

    	}

    	ro_scal->sobol_count_ini = sob_ini_count;

#ifdef TIME_RESOLUTION_MICROSECOND
        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
#else
        mlfi_ftime(&t_end);
        elapsed = mlfi_diff_time(t_end,t_start);
#endif

        printf("\n\nGPU Run Time: %lu !\n\n", elapsed);

    	for(int i=0; i<ro_scal->inst_num*ro_scal->num_contracts; i+=ro_scal->num_contracts) {
    		for(int j=0; j<ro_scal->num_contracts; j++) {
    			printf("FINAL RESULT GPU for(inst,contract_num): (%d,%d) IS: %.16f \n\n\n",
    				i/ro_scal->num_contracts, j, priv_arr->vhat_fin[i+j]/ro_scal->TOT_MC_ITER_NUM);
    		}
    	}
    }


    { // free allocated space
        shrLog("Release CPU buffers and OpenCL objects...\n");
        // clFinish(cqCommandQueue[0]);
        // clReleaseKernel(ckGenPricing);
        clReleaseProgram(cpProgram);
        {
            clReleaseMemObject(dev_arr.ro_scals);
            clReleaseMemObject(dev_arr.sobol_v_dir);
            clReleaseMemObject(dev_arr.sobol_fix_ind);
            clReleaseMemObject(dev_arr.bb_ia);
            clReleaseMemObject(dev_arr.bb_coefs);
            clReleaseMemObject(dev_arr.inst_coefs);
            clReleaseMemObject(dev_arr.pc_coefs);
            clReleaseMemObject(dev_arr.model_vhat);

            free(priv_arr->glb_vhat);

            clReleaseCommandQueue(cqCommandQueue[dev_id]);
        }
        free(cdDevices);
        clReleaseContext(cxGPUContext);
    }

    //exit(0);
}
