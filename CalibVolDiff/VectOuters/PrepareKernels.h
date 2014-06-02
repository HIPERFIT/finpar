#ifndef PREPARE_KERNELS
#define PREPARE_KERNELS

/************************************************************/
/***************  HELPER DATA STRUCTURES   ******************/
/************************************************************/

// GPU BUFFERS
struct oclNordeaArrays {
    /* RO scalars */
    cl_mem  ro_scals;

    /* RO arrays */
    cl_mem timeline; // [NUM_T]

    cl_mem myX;   // [NUM_X]
    cl_mem myDx;  // [NUM_X * 3]
    cl_mem myDxx; // [NUM_X * 3]

    cl_mem myY;   // [NUM_Y]
    cl_mem myDy;  // [NUM_Y * 3]
    cl_mem myDyy; // [NUM_Y * 3]

    /* TRIDAG helpers */
    cl_mem a;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem b;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem c;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem y;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem u;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem v;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem tmp4; // [4 * OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    cl_mem tmp2; // [2 * OUTER_LOOP_COUNT * NUM_X * NUM_Y]

    cl_mem res_arr; // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]

    /* WRITE ONLY OUTPUT stored in "res"*/
    //cl_mem res;
};

// CPU BUFFERS
struct NordeaArrays {
    /* RO arrays */
    REAL* timeline; // [NUM_T]

    REAL* myX;   // [NUM_X]
    REAL* myDx;  // [NUM_X * 3]
    REAL* myDxx; // [NUM_X * 3]

    REAL* myY;   // [NUM_Y]
    REAL* myDy;  // [NUM_Y * 3]
    REAL* myDyy; // [NUM_Y * 3]

    /* TRIDAG helpers are created on GPU only! no need for them to exist in main memory*/
    REAL* a;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* b;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* c;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* y;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* u;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* v;   // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]
    REAL* tmp; // [4 * OUTER_LOOP_COUNT * NUM_X * NUM_Y]

    REAL* res_arr; // [OUTER_LOOP_COUNT * NUM_X * NUM_Y]

    /* NOT USED FOR NOW! WRITE ONLY OUTPUT stored in "res"*/
    //REAL* res;
};

// GPU KERNELS
typedef struct {
    cl_kernel ckMatTranspUpdate;
    cl_kernel ckMatTransposeU;
    cl_kernel ckMatTransposeV;

    // TRIDAG KERNELS
    cl_kernel ckNordeaKernelX;  // OUTER_LOOP_COUNT x NUM_Y
    cl_kernel ckNordeaKernelY;	// OUTER_LOOP_COUNT x NUM_X

    size_t FORM_32X [3];
    size_t FORM_32Y [3];
    size_t FORM_Y32 [3];
} GPUkernels __attribute__ ((aligned (32)));

/******************************************************/
/******* compute global/local grid dimensions *********/
/******************************************************/

static bool is_pow2(unsigned int x) {
    while(x > 1) {
        if( x & 1 ) return false;
        x = x >> 1;
    }
    return true;
}

void initGPUgrid ( GPUkernels& kernels ) {
    // fill transposition layouts!

    bool is_safe =  ( (OUTER_LOOP_COUNT*NUM_X) % NUM_Y == 0          ) &&
    				//( (OUTER_LOOP_COUNT*NUM_Y) % NUM_X == 0          ) &&
      				( (OUTER_LOOP_COUNT*NUM_X) % WORKGROUP_SIZE == 0 ) &&
      				( (OUTER_LOOP_COUNT*NUM_Y) % WORKGROUP_SIZE == 0 ) &&
       				( NUM_X % 32 == 0 && NUM_Y % 32 == 0);
    assert( is_safe && "CANNOT EXECUTE ON GPU: sizes do not match!");

    kernels.FORM_32X[0] = 32;    kernels.FORM_32X[1] = NUM_X; kernels.FORM_32X[2] = (OUTER_LOOP_COUNT*NUM_Y)/32;
    kernels.FORM_32Y[0] = 32;    kernels.FORM_32Y[1] = NUM_Y; kernels.FORM_32Y[2] = (OUTER_LOOP_COUNT*NUM_X)/32;
    kernels.FORM_Y32[0] = NUM_Y; kernels.FORM_Y32[1] = 32;    kernels.FORM_Y32[2] = (OUTER_LOOP_COUNT*NUM_X)/32;
}

/************************************************************/
/*************** KERNELS OTHER THAN TRIDAG ******************
/************************************************************/

///////////////////////////
/// NORDEA_ALL_KERNEL_X ///
///////////////////////////

cl_kernel make_NordeaKernelX(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs
) {
    //size_t*             globalWorkSize,
    //size_t*             localWorkSize
    cl_kernel ckAllX = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckAllX = clCreateKernel(cpProgram, "nordea_kernel_x", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    // 1. RO SCALARS //
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.ro_scals   );

    // 2. RO ARRAYS //
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myX);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDx);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDxx);

    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myY);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDy);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDyy);

    // 3. GPU-ONLY GLOBAL ARRAYS //
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.b);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.y);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);

    // 4. OUT ARRAY!
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.res_arr );

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckAllX;
}

void run_NordeaKernelX(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    size_t  globalWorkSize = NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckNordeaKernelX, 1, NULL,
                                       &globalWorkSize, &localWorkSize, 0, NULL, NULL);
    ciErr1 |= clFinish(cqCommandQueue);

    oclCheckError(ciErr1, CL_SUCCESS);
}

///////////////////////////
/// NORDEA_ALL_KERNEL_Y ///
///////////////////////////

cl_kernel make_NordeaKernelY(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs
) {
    cl_kernel ckAllY = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckAllY = clCreateKernel(cpProgram, "nordea_kernel_y", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    // 1. RO SCALARS //
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.ro_scals   );

    // 2. RO ARRAYS //
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDy);
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDyy);

    // 3. GPU-ONLY GLOBAL ARRAYS //
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a);
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.b);
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c);

    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.y); // y <- u^T
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u); // u <- v^T
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckAllY;
}

void run_NordeaKernelY(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    size_t  globalWorkSize = NUM_X * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckNordeaKernelY, 1, NULL,
                                       &globalWorkSize, &localWorkSize, 0, NULL, NULL);

    ciErr1 |= clFinish(cqCommandQueue);

    oclCheckError(ciErr1, CL_SUCCESS);
}

/***************************************************/
/*********  Segmented Matrix Transposition  ********/
/***************************************************/

cl_kernel make_transposeGPU(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        cl_mem              inp_arr,
        cl_mem              out_arr,
        size_t*             size_x,
        size_t*             size_y,
        size_t*             size_z,
        const char*         kernel_name
//        size_t*             globalWorkSize,
//        size_t*             localWorkSize
) {
    cl_kernel ckMatTransp = NULL;
    cl_int    ciErr1, ciErr2;
    unsigned int  counter = 0;

    ckMatTransp = clCreateKernel(cpProgram, kernel_name, &ciErr1); // "transposeUpdateScalars"
    oclCheckError(ciErr1, CL_SUCCESS);

    ciErr1 |= clSetKernelArg(ckMatTransp, counter++, sizeof(cl_mem), (void *) &out_arr); // input array
    ciErr1 |= clSetKernelArg(ckMatTransp, counter++, sizeof(cl_mem), (void *) &inp_arr); // output array
    //ciErr1 |= clSetKernelArg(ckMatTransp, counter++, sizeof(int), &offset);
    ciErr1 |= clSetKernelArg(ckMatTransp, counter++, sizeof(int), size_x);              // size of 0 dim
    ciErr1 |= clSetKernelArg(ckMatTransp, counter++, sizeof(int), size_y);              // size of 1 dim

    // shared memory space
    ciErr1 |= clSetKernelArg(ckMatTransp, counter++, (BLOCK_DIM + 1) * BLOCK_DIM * sizeof(REAL), 0 );

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckMatTransp;
}

void run_transposeGPU (
        cl_command_queue&   cqCommandQueue,
        cl_kernel           kernel,
        size_t*             globalWorkSize

) {
    cl_int  ciErr1;
    size_t   localWorkSize [3]  = {BLOCK_DIM, BLOCK_DIM, 1};

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernel, 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );
    ciErr1 |= clFinish(cqCommandQueue);
    oclCheckError(ciErr1, CL_SUCCESS);
}

cl_kernel make_transposeGPU_WithUpdate(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        cl_mem              inp_arr,
        cl_mem              out_arr,
        size_t*             size_x,
        size_t*             size_y,
        size_t*             size_z,
        oclNordeaArrays&    ocl_arrs
//        size_t*             globalWorkSize,
//        size_t*             localWorkSize
) {
    cl_int    ciErr1;
    cl_kernel kernel = make_transposeGPU(
            cqCommandQueue, cpProgram, inp_arr,
            out_arr, size_x, size_y, size_z, "transposeUpdateScalars"
        );

    ciErr1 = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &ocl_arrs.ro_scals);
    ciErr1 = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &ocl_arrs.timeline);
    oclCheckError(ciErr1, CL_SUCCESS);

    return kernel;
}

void run_transposeGPU_WithUpdate (
        cl_command_queue&   cqCommandQueue,
        cl_kernel           kernel,
        size_t*             globalWorkSize
) {
    cl_int  ciErr1;
    size_t   localWorkSize [3]  = {BLOCK_DIM, BLOCK_DIM, 1};
    //size_t*  localWorkSize  = kernels.LWS_YXZ;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernel, 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );
    ciErr1 |= clFinish(cqCommandQueue);
    oclCheckError(ciErr1, CL_SUCCESS);
}

/**********************************************************/
/**************** MAKE ALL KERNELS ************************/
/**********************************************************/

void make_kernels (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs,
        GPUkernels&         kernels
) {
    initGPUgrid ( kernels );

    // First Kernel
    kernels.ckNordeaKernelX = make_NordeaKernelX (
            cqCommandQueue, cpProgram, ro_scal, ocl_arrs
        );

    // Second Kernel
    kernels.ckNordeaKernelY = make_NordeaKernelY (
            cqCommandQueue, cpProgram, ro_scal, ocl_arrs
        );

    { // matrix transpose!!!
    	// transpose U into Y
    	kernels.ckMatTransposeU = make_transposeGPU (
    			cqCommandQueue, cpProgram, ocl_arrs.u, ocl_arrs.y,
    			&kernels.FORM_32X[0], &kernels.FORM_32X[1], &kernels.FORM_32X[2], "transpose"
    			//&NUM_X, &NUM_Y, &OUTER_LOOP_COUNT, "transpose"
        	);

    	// transpose V into U
    	kernels.ckMatTransposeV = make_transposeGPU (
    			cqCommandQueue, cpProgram, ocl_arrs.v, ocl_arrs.u,
    			&kernels.FORM_Y32[0], &kernels.FORM_Y32[1], &kernels.FORM_Y32[2], "transpose"
    			//&NUM_X, &NUM_Y, &OUTER_LOOP_COUNT, "transpose"
        	);


    	// Transpose into res_arr and update induction variables!
    	kernels.ckMatTranspUpdate = make_transposeGPU_WithUpdate (
    			cqCommandQueue, cpProgram, ocl_arrs.u, ocl_arrs.res_arr,    // ocl_arrs.y
    			&kernels.FORM_32Y[0], &kernels.FORM_32Y[1], &kernels.FORM_32Y[2], ocl_arrs
    			//&NUM_Y, &NUM_X, &OUTER_LOOP_COUNT, ocl_arrs
        	);
    }
}


/**********************************************************/
/********************* RUN KERNELS ************************/
/**********************************************************/

void run_GPUkernels_one_time_iteration (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    size_t globalWorkSize[3];

    run_NordeaKernelX( cqCommandQueue, kernels );
    run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeU, kernels.FORM_32X ); // y <- transpose(u)
    run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeV, kernels.FORM_Y32 ); // u <- transpose(v)
    run_NordeaKernelY( cqCommandQueue, kernels );
    run_transposeGPU_WithUpdate( cqCommandQueue, kernels.ckMatTranspUpdate, kernels.FORM_32Y);
}

/**********************************************************/
/******************* RELEASE KERNELS **********************/
/**********************************************************/

void release_all_kernels ( GPUkernels& kernels ) {
    clReleaseKernel(kernels.ckNordeaKernelX);

    clReleaseKernel(kernels.ckMatTransposeU);
    clReleaseKernel(kernels.ckMatTransposeV);

    clReleaseKernel(kernels.ckNordeaKernelY);

    clReleaseKernel(kernels.ckMatTranspUpdate);
}

void release_all_GPU_resources(
        cl_command_queue&   cqCommandQueue,
        cl_context&         cxGPUContext,
        cl_program          cpProgram,
        cl_device_id*       cdDevices,
        oclNordeaArrays&    ocl_arrs,
        GPUkernels&         kernels
) {

    shrLog(stdlog, "Release Kernels...\n");
    release_all_kernels( kernels );

    shrLog(stdlog, "Release CPU buffers and OpenCL objects...\n");
    clReleaseProgram(cpProgram);

    clReleaseMemObject(ocl_arrs.ro_scals);
    clReleaseMemObject(ocl_arrs.timeline);

    clReleaseMemObject(ocl_arrs.a);
    clReleaseMemObject(ocl_arrs.b);
    clReleaseMemObject(ocl_arrs.c);

    clReleaseMemObject(ocl_arrs.y);
    clReleaseMemObject(ocl_arrs.u);
    clReleaseMemObject(ocl_arrs.v);

    clReleaseMemObject(ocl_arrs.res_arr);

    clReleaseMemObject(ocl_arrs.myX);
    clReleaseMemObject(ocl_arrs.myDx);
    clReleaseMemObject(ocl_arrs.myDxx);

    clReleaseMemObject(ocl_arrs.myY);
    clReleaseMemObject(ocl_arrs.myDy);
    clReleaseMemObject(ocl_arrs.myDyy);

    clReleaseCommandQueue(cqCommandQueue);
    free(cdDevices);
    clReleaseContext(cxGPUContext);
}


/*****************************************************************/
/************ DEBUGGING ******************************************/
/*****************************************************************/


void run_trimmed_GPUkernels_one_time_iteration (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        NordeaArrays&       cpu_arrs,
        oclNordeaArrays&    ocl_arrs,
        int                 time_ind,
        const REAL          alpha,
        const REAL          beta,
        const REAL          nu
) {
	run_NordeaKernelX( cqCommandQueue, kernels );
	run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeU, kernels.FORM_32X ); // y <- transpose(u)
    run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeV, kernels.FORM_Y32 ); // u <- transpose(v)
    run_NordeaKernelY( cqCommandQueue, kernels );
    run_transposeGPU_WithUpdate( cqCommandQueue, kernels.ckMatTranspUpdate, kernels.FORM_32Y);

    cl_int ciErr;
    const unsigned int ARR_SIZE = NUM_X * NUM_Y * OUTER_LOOP_COUNT * sizeof(REAL);
#if 1
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.a, CL_TRUE, 0, ARR_SIZE, cpu_arrs.a, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.b, CL_TRUE, 0, ARR_SIZE, cpu_arrs.b, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.c, CL_TRUE, 0, ARR_SIZE, cpu_arrs.c, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.u, CL_TRUE, 0, ARR_SIZE, cpu_arrs.u, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.v, CL_TRUE, 0, ARR_SIZE, cpu_arrs.v, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.y, CL_TRUE, 0, ARR_SIZE, cpu_arrs.y, 0, NULL, NULL );
    ciErr = clEnqueueReadBuffer( cqCommandQueue, ocl_arrs.res_arr, CL_TRUE, 0, ARR_SIZE, cpu_arrs.res_arr, 0, NULL, NULL );
#endif

    oclCheckError(ciErr, CL_SUCCESS);

    // now the rest of the code!
    unsigned int i, j, k;

    REAL dtInv = 1.0/(cpu_arrs.timeline[time_ind+1]-cpu_arrs.timeline[time_ind]);

    REAL   *res_arr = cpu_arrs.res_arr, *a = cpu_arrs.a,
           *b = cpu_arrs.b, *c = cpu_arrs.c, *u = cpu_arrs.u,
           *v = cpu_arrs.v, *y = cpu_arrs.y;

    unsigned int cos_test_ind = 3*NUM_XY + 23*NUM_X + 211, cos_test_ind_rev = 3*NUM_XY + 211*NUM_Y + 23;


    { // transpose u: u is in form: [i/32, j, i%32] 
      //   and needs to be brought to form [i/32, i%32, j], a.k.a. [i,j]
    	const unsigned int ID = OUTER_LOOP_COUNT * NUM_X / 32;
    	for( k=0; k<ID; ++k ) {  // segmented transpose followed by copy out!
    	    const unsigned int glb_ind = k*NUM_Y*32;

    	    for(j=0; j<NUM_Y; j++) {
    	        for(int i=0; i<32; i++) {
    	           res_arr[glb_ind + i*NUM_Y+j] = u[glb_ind + j*32+i]; //y[glb_ind + i*NUM_Y+j];
    	        }
    	    }
    	}
    }
    // segmented transpose followed by copy out!
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {  
        const unsigned int glb_ind = k*NUM_X*NUM_Y;

        for(i=0; i<NUM_X; i++) {
            for(int j=0; j<NUM_Y; j++) {
            	y[glb_ind + j*NUM_X + i] = res_arr[glb_ind + i*NUM_Y + j];
            }
        }

        for(i=0; i<NUM_X*NUM_Y; i++) res_arr[glb_ind + i] = y[glb_ind + i];
    }
}

#endif // include "PrepareKernels"
