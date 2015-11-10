#ifndef GPU_HELPERS
#define GPU_HELPERS

#include "Optimizations.h"
#include "SDK_stub.h"

#if (_OPTIMIZATION_USE_FLOATS)
    #define TILE    GPU_TILE
#else
    #define TILE    GPU_TILE/2
#endif

/****************************/
/*** GPU-Allocated Arrays ***/
/****************************/
struct oclLoopArrays {

    /* RO scalars */
    cl_mem  ro_scals;

    /* sobol direction vector; RO */
    cl_mem  sobol_dirvcts;  // UINT[sobol_bit_count][num_under*num_dates];
    cl_mem  sobol_fix_ind;  // UCHAR[chunk-1]

    /* Brownian Bridge data RO */
    cl_mem  bb_inds;        // UINT[3][num_dates]
    cl_mem  bb_data;        // REAL[3][num_dates]

    /* MD instance data; RO */       
    cl_mem  md_c;             // Includes all arrays bellow:
//    cl_mem  md_c            // RO; REAL[num_models][num_under^2]
//    cl_mem  md_vols;        // RO; REAL[num_models][num_under*num_dates]
//    cl_mem  md_drifts;      // RO; REAL[num_models][num_under*num_dates]
//    cl_mem  md_starts;      // RO; REAL[num_models][num_under]
//    cl_mem  md_discts;      // REAL[num_models][num_cash_flows ]
//    cl_mem  md_detvals;     // REAL[num_models][num_det_pricers]

    /* FOR VECTORIZED VERSION */
    cl_mem  md_zd;    // RW: REAL[num_gpuits]x[num_under*num_dates]
    cl_mem  md_z;  // RW: REAL[num_gpuits]x[num_models]x[num_under*num_dates]

    /* global array to be reduced on CPU; WO */
    cl_mem  model_vhat;  // WO; ~ [num_models]x[num_contracts]x[num_iter/(chunk*BLOCK)]

    /*** local tiled space ***/
    cl_mem 	local_space; 

    void cleanupPRIV() {
        clReleaseMemObject(ro_scals);
        clReleaseMemObject(sobol_dirvcts);
        clReleaseMemObject(sobol_fix_ind);
        clReleaseMemObject(bb_inds);
        clReleaseMemObject(bb_data);
        clReleaseMemObject(md_c);
        clReleaseMemObject(model_vhat);
    }
    void cleanupVECT() {
        cleanupPRIV();
        clReleaseMemObject(md_z );
        clReleaseMemObject(md_zd);
    }
};

/*********************************/
/***  Math-Utility Functions   ***/
/*********************************/

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
	while(res <= n) res += mul;
	return res - mul;
}

/*********************************/
/*** GPU-Glue Helper Functions ***/
/*********************************/

enum GPU_KERNEL { PRIV, VECT };

GPU_KERNEL 
priv_or_vect_kernel(const LoopROScalars& ro_scal) {
#ifndef _OPTIMIZATION_MEM_COALES_ON
    return VECT;
#endif
#if   (_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV == 1)
    return VECT;
#elif (_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV == 2)
    return PRIV;
#else
    UINT threshold = 16;
    bool is_priv = (ro_scal.num_under*ro_scal.num_dates+ro_scal.num_models <= threshold);
    fprintf(stderr, "\nCOST MODEL THRESHOLD: %u, PRIV_KERNEL: %d\n", threshold, (int)is_priv);
    return (is_priv) ? PRIV : VECT;
#endif
}

const char* 
makeGPUprogPreamble(
        const LoopROScalars&    ro_scal,
        const GPU_KERNEL   &    kernel_type
) {
    const char* preamble =  (ro_scal.contract == 1 && kernel_type == VECT) ?
                                "#define VECT_KERNEL\n#include \"ContractDefs/SmallContract.cl\"\n\n" :
                            (ro_scal.contract == 1 && kernel_type == PRIV) ?
                                "#include \"ContractDefs/SmallContract.cl\"\n\n" :
                            (ro_scal.contract == 2 && kernel_type == VECT) ?
                                "#define VECT_KERNEL\n#include \"ContractDefs/MediumContract.cl\"\n\n" :
                            (ro_scal.contract == 2 && kernel_type == PRIV) ?
                                "#include \"ContractDefs/MediumContract.cl\"\n\n" :
                            (ro_scal.contract == 3 && kernel_type == VECT) ?
                                "#define VECT_KERNEL\n#include \"ContractDefs/LargeContract.cl\"\n\n" :
                            (ro_scal.contract == 3 && kernel_type == PRIV) ?
                                "#include \"ContractDefs/LargeContract.cl\"\n\n" : NULL;
    assert(preamble && "Invalid Contract Number (should be 1, 2 or 3) OR Invalid kernel type (PRIV or VECT)\n");
    return preamble;
}


UINT 
getWorkSize(const UINT& num_iter, const UINT& chunk, const UINT& block) {
	UINT tmp_rem = num_iter % chunk;
	UINT res     = (tmp_rem == 0) ? 
                         (num_iter / chunk)      : 
                        ((num_iter / chunk) + 1) ;

	if(res % block == 0) return res;

	res = (res / block) * block + block;
	return res;
}

int
discriminate_cost_model( LoopROScalars& ro_scal, cl_device_id device, const GPU_KERNEL kernel) {
    bool is_priv = (kernel == PRIV);
    cl_uint tot_cuda_cores = 0;

    UINT size_vct = ro_scal.num_under * ro_scal.num_dates;

    ro_scal.logBLOCK = 7;
    ro_scal.BLOCK = 1 << ro_scal.logBLOCK;

    // If not-enough GPU resources => sequentialize 
    //  the loop into big chunks that execute on GPU.
    if(!is_priv) { 
        // get device memory size info 
        cl_ulong glob_mem_size;
        clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(glob_mem_size),&glob_mem_size,NULL);
        { // check global mem matches the user-declared one
            const size_t USER_MEM = ((size_t)GPU_GLB_MEM) << 30;
            const size_t ub       = static_cast<size_t>( USER_MEM * 101.0 / 100.0 );
            const size_t lb       = static_cast<size_t>( USER_MEM *  99.0 / 100.0 );
            if ( !(glob_mem_size >= lb && glob_mem_size <= ub) ) {
                fprintf(stderr, "WARNING! Querried GPU global memory %lu DIFFERS from what user declared: [%ld,%ld]!\n", 
			(unsigned long)glob_mem_size, lb, ub);
                glob_mem_size = static_cast<cl_ulong>(USER_MEM);
                fprintf(stderr, "Using user-declared size for global memory: %lu\n", (unsigned long)glob_mem_size);
            }
        } 

        // check if there is enough memory
    	UINT perit_mem_size = size_vct*( ro_scal.num_models + 1 ) + ro_scal.num_models;
    	perit_mem_size *= sizeof(REAL);

        UINT max_mc_iter_num = ( (glob_mem_size * 3)/4 ) / perit_mem_size;

        UINT BLOCK = ro_scal.BLOCK;
        if(max_mc_iter_num % BLOCK != 0)
            max_mc_iter_num = (max_mc_iter_num / BLOCK)*BLOCK + BLOCK;

        max_mc_iter_num = prevMultipleOf(max_mc_iter_num, ro_scal.chunk);

        ro_scal.num_gpuits = MIN(max_mc_iter_num, ro_scal.num_mcits);

//        fprintf(stderr, "# of total iters: %u, # of gpu iters: %u, chunk: %u\n", 
//                        ro_scal.num_mcits, ro_scal.num_gpuits, ro_scal.chunk );
    } else {
        ro_scal.num_gpuits = ro_scal.num_mcits;
    }

    {
#ifdef __APPLE__
      tot_cuda_cores = GPU_CORES;   // MEMO: No correctness check... mael 2015-03-12
#else
        cl_uint compute_units;
        clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(compute_units), &compute_units, NULL);
        cl_uint comp_capabil_major, comp_capabil_minor;
        clGetDeviceInfo(device, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                            sizeof(cl_uint), &comp_capabil_major, NULL);
        clGetDeviceInfo(device, CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,
                            sizeof(cl_uint), &comp_capabil_minor, NULL);
        tot_cuda_cores = compute_units *
            ConvertSMVer2CoresCopy(comp_capabil_major, comp_capabil_minor);

        if( tot_cuda_cores != GPU_CORES ) {
            fprintf(stderr, "WARNING! Querried GPU number of cores %u DIFFERS from what user declared: %d!\n", 
                             tot_cuda_cores, GPU_CORES );
            tot_cuda_cores = GPU_CORES;
            fprintf(stderr, "Using user-declared GPU number of cores: %u\n", tot_cuda_cores);
        }
#endif
        UINT div = (1 << logNextPow2(tot_cuda_cores)) * 64;
        ro_scal.log_chunk = logNextPow2(ro_scal.num_gpuits / div );  // was div*2
        ro_scal.log_chunk = MIN(ro_scal.log_chunk, logMAX_CHUNK);
        if(ro_scal.log_chunk <= 3) ro_scal.log_chunk++;
    }

    ro_scal.chunk = 1 << ro_scal.log_chunk;

//    fprintf(stderr, "\n\nCOST MODEL COMPUTED: MC_ITERS: %u, MC_TOT_ITERS: %u, CHUNK: %u, BLOCK: %u\n\n",
//                    ro_scal.num_gpuits, ro_scal.num_mcits, ro_scal.chunk, ro_scal.BLOCK);

    return tot_cuda_cores;
}

/*********************************/
/*** Allocating GPU Buffers &  ***/
/***running the Private Kernels***/
/*********************************/

REAL*
oclAllocArrays_PrivKernel (
        oclLoopArrays&          ocl_arrs,
        cl_context              cxGPUContext,
        cl_command_queue&       cqCommandQueue,
        cl_int&                 ciErr,
        const LoopROScalars&    ro_scals,
        const SobolArrays  &    sob_arrs,
        const BrowBridgeArrays&  bb_arrs,
        const ModelArrays  &     md_arrs
) {
    cl_int ciErr2;
    size_t cur_size;
    REAL*  glb_vhat = NULL;
    ciErr = 0;
#if 0
    { // 1. RO scalars
        cur_size = sizeof(LoopROScalars);
        ocl_arrs.ro_scals = clCreateBuffer(
                cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
            );

        ciErr  = ciErr2;
        ciErr |= clEnqueueWriteBuffer(cqCommandQueue, ocl_arrs.ro_scals, CL_TRUE, 0,
                                        cur_size, &ro_scals, 0, NULL, NULL);
    }
#endif

    { // 2. Sobol RO Arrays
        const SobolArrays& cpu_arrs = sob_arrs;

        cur_size = (UINT)   ro_scals.num_under  * ro_scals.num_dates *
                            ro_scals.sobol_bits * sizeof(UINT);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(sobol_dirvcts); 
        
        cur_size+= ro_scals.chunk * sizeof(UCHAR);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(sobol_fix_ind); 
    }

    { // 3. BROWNIAN BRIDGE RO indirect arrays and data arrays
        const BrowBridgeArrays& cpu_arrs = bb_arrs;

        cur_size = 3*ro_scals.num_dates*sizeof(UINT);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(bb_inds);

        cur_size = 3*ro_scals.num_dates*sizeof(REAL);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(bb_data);
    }

    { // 4. Various RO Model Data Arrays
        const ModelArrays& cpu_arrs = md_arrs;
        
        cur_size  = ro_scals.num_under * ( ro_scals.num_under + 2*ro_scals.num_dates + 1);
        cur_size += ro_scals.num_cash_flows + ro_scals.num_det_pricers;
        cur_size *= ro_scals.num_models * sizeof(REAL);

        CREATE_AND_ENQUEUE_RO_CL_BUFFER(md_c);
    }

    { // 5. Finally, the WO Partial Result
        //red_sz = getGlobalRedArrSize( ro_scal->mc_iter_num, ro_scal->CHUNK*ro_scal->BLOCK );
        cur_size = getWorkSize(ro_scals.num_gpuits, ro_scals.chunk, ro_scals.BLOCK) / ro_scals.BLOCK;
        cur_size *= ro_scals.num_models * sizeof(REAL);

        glb_vhat = (REAL*)malloc(cur_size);
        ocl_arrs.model_vhat = clCreateBuffer(
                cxGPUContext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
                cur_size, glb_vhat, &ciErr2
            );
        ciErr |= ciErr2;
    }

    // 2. the rest: priv_arr->sobol_last_num_vec/md_zd/inst_trajWF
    //      are implemented as __local (private) hence do not require
    //      global memory allocation.
    
    oclCheckError(ciErr, CL_SUCCESS);
    
    return glb_vhat;
}


void
runGPU_PRIV(
        const LoopROScalars&    ro_scals,
        REAL*                   glb_vhat,
        oclLoopArrays&          dev_arr,
        cl_command_queue&       cqCommandQueue,
        cl_program&             cpProgram,
        size_t*                 globalWorkSize,
        size_t*                 localWorkSize
) {
    cl_kernel ckGenPricing = NULL; // OpenCL kernel
    cl_int    ciErr1;

    size_t  counter = 0, priv_sz = 0; // ro_scals
    const int PRIV_MULT   = sizeof(REAL)*ro_scals.BLOCK;

    ckGenPricing = clCreateKernel(cpProgram, "payoffGPU", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);

    // 1. RO SCALARS struct
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals   );

    // 2. RO SOBOL ARR //
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_dirvcts);
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_fix_ind);

    // 3. RO BROWNIAN BRIDGE //
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_inds);
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_data);

    // 4. RO MD INSTANCE DATA //
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.md_c);

    // 5. GLOBAL WO Array //
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, sizeof(cl_mem), (void*)&dev_arr.model_vhat );

    // 6. LOCAL (PRIVATE ARRAYS) MD_Z
    priv_sz = PRIV_MULT*ro_scals.num_under*ro_scals.num_dates;
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, priv_sz, NULL);

    // 7. LOCAL (PRIVATE ARRAY) VHAT
    priv_sz = PRIV_MULT * ro_scals.num_models;
    ciErr1 |= clSetKernelArg(ckGenPricing, counter++, priv_sz, NULL); // model_vhat (the private one)

    //oclCheckError(ciErr1, CL_SUCCESS);
    //clFinish(cqCommandQueue);

#if 0
    struct timeval t_start, t_end, t_diff; unsigned long long elapsed;
    gettimeofday(&t_start, NULL);
#endif
    // 8. ENQUEUE KERNEL!!!  //
    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, ckGenPricing, 1, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL);
    //oclCheckError(ciErr1, CL_SUCCESS);

    clFinish(cqCommandQueue);
#if 0
    gettimeofday(&t_end, NULL);
    timeval_subtract(&t_diff, &t_end, &t_start);
    elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
    printf("KERNEL TIME IS: %llu\n\n", elapsed);

    gettimeofday(&t_start, NULL);
#endif
    // 9. FINALLY, WRITE BACK!!! //
    priv_sz = (globalWorkSize[0]/ro_scals.BLOCK)*ro_scals.num_models*sizeof(REAL);

    ciErr1 |= clEnqueueReadBuffer(cqCommandQueue, dev_arr.model_vhat, CL_TRUE, 0,
                                    priv_sz, glb_vhat, 0, NULL, NULL);
    oclCheckError(ciErr1, CL_SUCCESS);

    clReleaseKernel(ckGenPricing);

    clFinish(cqCommandQueue);
#if 0
    gettimeofday(&t_end, NULL);
    timeval_subtract(&t_diff, &t_end, &t_start);
    elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
    printf("COPY-OUT and KERNEL RELEASE TIME IS: %llu\n\n", elapsed);
#endif
}


void reduceVHAT_CPU(const LoopROScalars& ro_scals, REAL* vhat, double* cpu_vhat ) {
    // final reduction on CPU!!!
    UINT SZ = getWorkSize(ro_scals.num_gpuits, ro_scals.chunk, ro_scals.BLOCK) / ro_scals.BLOCK;

    for( UINT ii = 0; ii<ro_scals.num_models; ii++) {
        int    offset    = ii * SZ;
        double final_res = 0.0;

        // this forms an inner dimension!
        for( UINT k = 0; k < SZ; k ++ ) {
            final_res += vhat[offset+k];
        }
        
        cpu_vhat[ii] += final_res;
    }
}


/*****************************/
/***** VECTORIZED KERNEL *****/
/*****************************/

inline REAL* 
oclAllocArrays_VectKernel (
        oclLoopArrays&          ocl_arrs,
        cl_context              cxGPUContext,
        cl_command_queue&       cqCommandQueue,
        cl_int&                 ciErr,
        const LoopROScalars&    ro_scals,
        const SobolArrays  &    sob_arrs,
        const BrowBridgeArrays&  bb_arrs,
        const ModelArrays  &     md_arrs
) {
    cl_int ciErr2;
    size_t cur_size;
    REAL*  glb_vhat;
    ciErr = 0;
#if 0
    { // 1. RO scalars
        cur_size = sizeof(LoopROScalars);
        ocl_arrs.ro_scals = clCreateBuffer(
                                cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
                            );
        ciErr  = ciErr2;
        ciErr |= clEnqueueWriteBuffer(  cqCommandQueue, ocl_arrs.ro_scals, CL_TRUE, 0,
                                        cur_size, &ro_scals, 0, NULL, NULL);
    }
#endif
    { // 2. SOBOL DIR VECTOR
        const SobolArrays& cpu_arrs = sob_arrs;

        cur_size = (UINT)   ro_scals.num_under  * ro_scals.num_dates *
                            ro_scals.sobol_bits * sizeof(UINT);
        //CREATE_AND_ENQUEUE_RO_CL_BUFFER(sobol_dirvctsT);
        ocl_arrs.sobol_dirvcts = 
            clCreateBuffer( cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            cur_size, sob_arrs.sobol_dirvctsT, &ciErr2 );
        ciErr |= ciErr2;

        cur_size = ro_scals.chunk * sizeof(UCHAR);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(sobol_fix_ind); 
    }

    { // 3. BROWNIAN BRIDGE RO indirect arrays and data arrays
        const BrowBridgeArrays& cpu_arrs = bb_arrs;

        cur_size = 3*ro_scals.num_dates*sizeof(UINT);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(bb_inds);

        cur_size = 3*ro_scals.num_dates*sizeof(REAL);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(bb_data);
    }

    { // 4. Various RO Model Data Arrays
        const ModelArrays& cpu_arrs = md_arrs;
        
        cur_size  = ro_scals.num_under * ( ro_scals.num_under + 2*ro_scals.num_dates + 1);
        cur_size += ro_scals.num_cash_flows + ro_scals.num_det_pricers;
        cur_size *= ro_scals.num_models * sizeof(REAL);

        CREATE_AND_ENQUEUE_RO_CL_BUFFER(md_c);
    }

    { // 5. PRIV: md_zd   is used for sob_vct, md_zd and trajectory vector 
      //          md_z    is filled up by Brownian Bridge & used in trajectory..
        cur_size =  getWorkSize( ro_scals.num_gpuits, ro_scals.chunk, ro_scals.BLOCK ) *
                    ro_scals.num_models * ro_scals.num_under * ro_scals.num_dates * 
                    ro_scals.chunk * sizeof(REAL);

        ocl_arrs.md_zd = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, 
                                        cur_size,     NULL,    &ciErr2 );
        ciErr |= ciErr2;

        ocl_arrs.md_z  = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE,
                                        cur_size,     NULL,    &ciErr2 );
        ciErr |= ciErr2;
    }

    oclCheckError(ciErr, CL_SUCCESS);

    { // 7. FINALLY, GLOBAL VHAT WO
        cur_size = getWorkSize( ro_scals.num_gpuits, ro_scals.chunk, ro_scals.BLOCK );
        cur_size*= ro_scals.num_models * sizeof(REAL);

        glb_vhat = (REAL*)malloc( cur_size );

        ocl_arrs.model_vhat = clCreateBuffer(
                cxGPUContext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, 
                cur_size, glb_vhat, &ciErr2
            );
        ciErr |= ciErr2;
        oclCheckError(ciErr, CL_SUCCESS);
    }

    return glb_vhat;
}

void runGPU_VECT(
        const LoopROScalars&    ro_scals,
        REAL*                   glb_vhat,
        oclLoopArrays&          dev_arr,
        cl_command_queue&       cqCommandQueue,
        cl_program&             cpProgram,
        size_t*                 globalWorkSize,
        size_t*                 localWorkSize
) {
    cl_kernel ckGenPricing_sobol = NULL, ckGenPricing_invg   = NULL, ckGenPricing_brown = NULL;
    cl_kernel ckGenPricing_traj  = NULL, ckGenPricing_finred = NULL;
    cl_int    ciErr1;

    int  counter = 0;

    { // SOBOL KERNEL!
        counter = 0;
        ckGenPricing_sobol = clCreateKernel(cpProgram, "mlfi_genmatrix_uniform2", &ciErr1);

        ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals     );
        ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_fix_ind);

        ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.sobol_dirvcts);
        ciErr1 |= clSetKernelArg(ckGenPricing_sobol, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd        );

        ciErr1 |= clEnqueueNDRangeKernel( cqCommandQueue, ckGenPricing_sobol, 1, NULL,
                                          globalWorkSize, localWorkSize, 0, NULL, NULL );
        oclCheckError(ciErr1, CL_SUCCESS);
    }

    { // INV GAUSSIAN KERNEL!
        counter = 0;
        ckGenPricing_invg = clCreateKernel(cpProgram, "mlfi_ugaussian_Pinv_vector1", &ciErr1);

        ciErr1 |= clSetKernelArg(ckGenPricing_invg, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);
        ciErr1 |= clSetKernelArg(ckGenPricing_invg, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd   );

        ciErr1 |= clEnqueueNDRangeKernel(   cqCommandQueue, ckGenPricing_invg, 1, NULL,
                                            globalWorkSize, localWorkSize, 0, NULL, NULL  );
        oclCheckError(ciErr1, CL_SUCCESS);
    }

    { // BROWNIAN BRIDGE KERNEL!
        //shrLog("Call Brownian Bridge kernel on GPU...\n\n");
        counter = 0;
        ckGenPricing_brown = clCreateKernel(cpProgram, "mlfi_brownianbridge_wiener_path1", &ciErr1);

        ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);
        ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_inds );
        ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.bb_data );

        ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd   );
        ciErr1 |= clSetKernelArg(ckGenPricing_brown, counter++, sizeof(cl_mem), (void*)&dev_arr.md_z    );

        ciErr1 |= clEnqueueNDRangeKernel(   cqCommandQueue, ckGenPricing_brown, 1, NULL,
                                            globalWorkSize, localWorkSize, 0, NULL, NULL );
        oclCheckError(ciErr1, CL_SUCCESS);
    }

    { // compute trajectory
        bool optim_kernel = (2*ro_scals.num_under <= TILE);
        counter = 0;

#ifdef _OPTIMIZATION_MEM_COALES_ON
        ckGenPricing_traj = (optim_kernel) ?
                clCreateKernel(cpProgram, "mlfi_comp_traj1",      &ciErr1) :
                clCreateKernel(cpProgram, "mlfi_comp_traj_unopt", &ciErr1) ;
#else
        ckGenPricing_traj = clCreateKernel(cpProgram, "mlfi_comp_traj1", &ciErr1);
#endif

        ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals);
        ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.md_c    );

#ifdef _OPTIMIZATION_MEM_COALES_ON
        if( optim_kernel ) {
            assert( ( 2*ro_scals.num_under <= TILE) && "Illegal dataset: num_under > 4" );
            UINT priv_sz = TILE * ro_scals.BLOCK * sizeof(REAL);
            ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, priv_sz, NULL); // local_space!
        }
#endif

        ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.md_z );
        ciErr1 |= clSetKernelArg(ckGenPricing_traj, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd);

        ciErr1 |= clEnqueueNDRangeKernel(   cqCommandQueue, ckGenPricing_traj, 1, NULL,
                                            globalWorkSize, localWorkSize, 0, NULL, NULL );

        
        oclCheckError(ciErr1, CL_SUCCESS);
    }

    { // final reduction!
        //shrLog("Call Compute Reduction kernel on GPU...\n\n");
        counter = 0;
        ckGenPricing_finred = clCreateKernel(cpProgram, "mlfi_reduction_step1", &ciErr1);
        oclCheckError(ciErr1, CL_SUCCESS);

        ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.ro_scals  );
        ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.md_c      );

        ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.md_zd      );
        ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, sizeof(cl_mem), (void*)&dev_arr.model_vhat );

        assert( (ro_scals.num_models <= TILE) && "Illegal dataset: num_models > 8" );
        UINT priv_sz = TILE * ro_scals.BLOCK * sizeof(REAL);
        ciErr1 |= clSetKernelArg(ckGenPricing_finred, counter++, priv_sz, NULL); // local_space!
        
        ciErr1 |= clEnqueueNDRangeKernel(   cqCommandQueue, ckGenPricing_finred, 1, NULL,
                                            globalWorkSize, localWorkSize, 0, NULL, NULL );
        oclCheckError(ciErr1, CL_SUCCESS);

        ciErr1 |= clEnqueueReadBuffer(  cqCommandQueue, dev_arr.model_vhat, CL_TRUE, 0,
                                        globalWorkSize[0] * ro_scals.num_models * sizeof(REAL),
                                        glb_vhat, 0, NULL, NULL);
        oclCheckError(ciErr1, CL_SUCCESS);
    }

    clReleaseKernel(ckGenPricing_sobol );
    clReleaseKernel(ckGenPricing_invg  );
    clReleaseKernel(ckGenPricing_brown );
    clReleaseKernel(ckGenPricing_traj  );
    clReleaseKernel(ckGenPricing_finred);
}

#endif // GPU_HELPERS

