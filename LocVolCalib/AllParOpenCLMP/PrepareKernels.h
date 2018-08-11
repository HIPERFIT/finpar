/***********************************************************/
/*************** Helper DATA STRUCTURES ********************/
/***********************************************************/

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
    // PRE/POST TRIDAG KERNELS
    cl_kernel ckPreTridagX;  // NUM_X x NUM_Y x OUTER_LOOP_COUNT
    cl_kernel ckPreTridagY;  // NUM_Y x NUM_X x OUTER_LOOP_COUNT
    cl_kernel ckMatTranspUpdate;

    cl_kernel ckMatTransposeU;
    cl_kernel ckMatTransposeV;

    // TRIDAG KERNELS
    cl_kernel ckTridagMatMult  [2];     // scan with matrix multiplication:
    cl_kernel ckConcludeMatMult[2];

    cl_kernel ckPreludeFwdFunComp[2];
    cl_kernel ckTridagFwdFunComp[2];   // scan with linear function composition
    cl_kernel ckPostFwdFunComp[2];

    cl_kernel ckPreludeBwdFunComp[2];
    cl_kernel ckPostBwdFunComp[2];

    cl_kernel ckTridagAllOpt[2];

    cl_kernel ckNordeaKernelX;
    cl_kernel ckNordeaKernelY;

    size_t GWS_XYZ [3];
    size_t LWS_XYZ [3];

    size_t GWS_YXZ [3];
    size_t LWS_YXZ [3];

    // NOT USED
    //cl_kernel ckMatTransp;   // NUM_Y x NUM_X x OUTER_LOOP_COUNT
} GPUkernels __attribute__ ((aligned (32)));

/******************************************************/
/******* compute global/local grid dimensions *********/
/******************************************************/

static void local_block_size(
        size_t size_x, size_t size_y, size_t size_z,
        size_t& loc_x, size_t& loc_y, size_t& loc_z
) {
    bool ok = is_pow2(size_x) && is_pow2(size_y);
    assert( ok && "Global Size X or Y is not a power of 2!");
    assert( size_x * size_y >= WORKGROUP_SIZE && "NUM_X * NUM_Y < WORKGROUP_SIZE!!!");

    loc_x = WORKGROUP_SIZE; loc_y = 1; loc_z = 1;
    while(size_x / loc_x == 0) {
        loc_x /= 2;
        loc_y *= 2;
    }

    ok = is_pow2(loc_x) && is_pow2(loc_y) && (loc_x*loc_y == WORKGROUP_SIZE) && (loc_x <= size_x) && (loc_y <= size_y);
    assert( ok && "Local Block Size is not a power of 2!");
}


void initGPUgrid ( GPUkernels& kernels ) {
    { // find a convenient local block size
        size_t loc_x, loc_y, loc_z;

        kernels.GWS_XYZ[0] = NUM_X; kernels.GWS_XYZ[1] = NUM_Y; kernels.GWS_XYZ[2] = OUTER_LOOP_COUNT;
        kernels.GWS_YXZ[0] = NUM_Y; kernels.GWS_YXZ[1] = NUM_X; kernels.GWS_YXZ[2] = OUTER_LOOP_COUNT;

        local_block_size(NUM_X, NUM_Y, OUTER_LOOP_COUNT, loc_x, loc_y, loc_z);
        kernels.LWS_XYZ[0] = loc_x; kernels.LWS_XYZ[1] = loc_y; kernels.LWS_XYZ[2] = loc_z;

        local_block_size(NUM_Y, NUM_X, OUTER_LOOP_COUNT, loc_y, loc_x, loc_z);
        kernels.LWS_YXZ[0] = loc_y; kernels.LWS_YXZ[1] = loc_x; kernels.LWS_YXZ[2] = loc_z;
    }

}

/************************************************************/
/*************** KERNELS OTHER THAN TRIDAG ******************
/************************************************************/

cl_kernel make_prepare_tridag_X(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs,
        size_t*             globalWorkSize,
        size_t*             localWorkSize
) {
    cl_kernel ckPreTridagX = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPreTridagX = clCreateKernel(cpProgram, "prepare_tridag_x", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    { // 1. RO SCALARS //
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.ro_scals   );
    }

    { // 2. RO ARRAYS //
        //ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.timeline);

        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myX);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDx);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDxx);

        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myY);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDy);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDyy);
    }

    { // 3. GPU-ONLY GLOBAL ARRAYS //
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.b);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.y);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp4);
    }

#if 1
    { // 4. LOCAL (PRIVATE ARRAYS)
        const int PRIV_MULT = sizeof(REAL) * localWorkSize[0] * localWorkSize[1] * localWorkSize[2];
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, PRIV_MULT, NULL);
    }
#endif
    { // 5. OUT ARRAY!
        ciErr1 |= clSetKernelArg(ckPreTridagX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.res_arr );
    }
#if 0
    { // Finally, enqueue kernel
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPreTridagX, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPreTridagX;
}

void run_prepare_tridag_X(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPreTridagX, 3, NULL,
                                         kernels.GWS_XYZ, kernels.LWS_XYZ, 0, NULL, NULL);
    oclCheckError(ciErr1, CL_SUCCESS);
}


///////////////////////////
/// NORDEA_ALL_KERNEL_X ///
///////////////////////////

cl_kernel make_NordeaKernelX(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs
) {
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
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u);
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);

    // 4. LOCAL (PRIVATE ARRAYS)
    const int PRIV_MULT = sizeof(REAL) * 8 * WORKGROUP_SIZE;
    ciErr1 |= clSetKernelArg(ckAllX, counter++, PRIV_MULT, NULL);

    // 5. OUT ARRAY!
    ciErr1 |= clSetKernelArg(ckAllX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.res_arr );

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckAllX;
}

void run_NordeaKernelX(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckNordeaKernelX, 1, NULL,
                                       &globalWorkSize, &localWorkSize, 0, NULL, NULL);


    oclCheckError(ciErr1, CL_SUCCESS);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

cl_kernel make_prepare_tridag_Y(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        size_t*             globalWorkSize,
        size_t*             localWorkSize
) {
    cl_kernel ckPreTridagY = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPreTridagY = clCreateKernel(cpProgram, "prepare_tridag_y", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    { // 1. RO SCALARS // here we do not have to reload ro_scals!
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.ro_scals   );
    }

    { // 2. RO ARRAYS //
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDy);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.myDyy);
    }

    { // 3. GPU-ONLY GLOBAL ARRAYS //
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.b);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.y);
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp4);
    }
#if 1
    { // 4. LOCAL (PRIVATE ARRAYS)
        const int PRIV_MULT = sizeof(REAL) * localWorkSize[0] * localWorkSize[1] * localWorkSize[2];
        ciErr1 |= clSetKernelArg(ckPreTridagY, counter++, PRIV_MULT, NULL);
    }
#endif
#if 0
    { // Finally, enqueue kernel
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPreTridagY, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPreTridagY;
}

void run_prepare_tridag_Y(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPreTridagY, 3, NULL,
                                         &kernels.GWS_YXZ[0], &kernels.LWS_YXZ[0], 0, NULL, NULL);
    oclCheckError(ciErr1, CL_SUCCESS);
}



///////////////////////////
/// NORDEA_ALL_KERNEL_X ///
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
#if TRANSPOSE_UV
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.y); // y <- u^T
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u); // u <- v^T
#else
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.u);
    ciErr1 |= clSetKernelArg(ckAllY, counter++, sizeof(cl_mem), (void*)&ocl_arrs.v);
#endif

    // 4. LOCAL (PRIVATE ARRAYS)
    const int PRIV_MULT = sizeof(REAL) * 8 * WORKGROUP_SIZE;
    ciErr1 |= clSetKernelArg(ckAllY, counter++, PRIV_MULT, NULL);

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckAllY;
}

void run_NordeaKernelY(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
    cl_int    ciErr1;
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckNordeaKernelY, 1, NULL,
                                       &globalWorkSize, &localWorkSize, 0, NULL, NULL);


    oclCheckError(ciErr1, CL_SUCCESS);
}


/************************************************************/
/******************** TRIDAG KERNELS ************************/
/************************************************************/



/**
 * Parallel Prefix with matrix multiplication
 */
cl_kernel make_MatMultScan(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        const unsigned int* NUM_P
) {
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;
    cl_kernel ckTridagMatMultX = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckTridagMatMultX = clCreateKernel(cpProgram, "scanMatMultIncl", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    { // 1. GPU-ONLY GLOBAL ARRAYS /(4 * NUM_X * NUM_Y * OUTER_LOOP_COUNT)
        ciErr1 |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp4);
    }

    { // 2. LOCAL: 8 * sizeof(REAL) * worksize
        const int PRIV_MULT = sizeof(REAL) * 8 * localWorkSize;
        ciErr1 |= clSetKernelArg(ckTridagMatMultX, counter++, PRIV_MULT, NULL);
    }

    { // 1. GPU-ONLY GLOBAL ARRAYS //
        ciErr1 |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(uint), (void *)NUM_P);
    }

#if 0
    { // Finally, the array size:
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckTridagMatMultX, 1, NULL,
                                         &globalWorkSize, &localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckTridagMatMultX;
}

void run_MatMultScan(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int    ciErr1;
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckTridagMatMult[X], 1, NULL,
                                        &globalWorkSize, &localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}




/**
 * Parallel Prefix with linear function composition
 */
cl_kernel make_FwdFunCompScan(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        const unsigned int* NUM_P
) {
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;
    cl_kernel ckTridagFwdFunComp = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0;
    unsigned long int cur_size;

    // find kernel
    ckTridagFwdFunComp = clCreateKernel(cpProgram, "scanLinFunCompIncl", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);

    // 1. GPU-ONLY GLOBAL ARRAYS /(4 * NUM_X * NUM_Y * OUTER_LOOP_COUNT)
    ciErr1 |= clSetKernelArg(ckTridagFwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp2);

    // 2. LOCAL: 4 * sizeof(REAL) * worksize
    const int PRIV_MULT = sizeof(REAL) * 4 * localWorkSize;
    ciErr1 |= clSetKernelArg(ckTridagFwdFunComp, counter++, PRIV_MULT, NULL);

    // 1. GPU-ONLY GLOBAL ARRAYS //
    ciErr1 |= clSetKernelArg(ckTridagFwdFunComp, counter++, sizeof(uint), (void *)NUM_P);

    // Finally, the array size:
#if 0
    ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckTridagFwdFunComp, 1, NULL,
                                     &globalWorkSize, &localWorkSize, 0, NULL, NULL);
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckTridagFwdFunComp;
}

void run_FwdFunCompScan(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckTridagFwdFunComp[X], 1, NULL,
                                        &globalWorkSize, &localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}

/**
 * Mapping back the results of parallel prefix with matrix multiplication
 */
cl_kernel make_conclude_matmult (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              d,
        cl_mem              y,
        cl_mem              u
) {
    cl_kernel ckConcludeMatMult = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckConcludeMatMult = clCreateKernel(cpProgram, "conclude_matmult", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckConcludeMatMult, counter++, sizeof(cl_mem), (void*)&ocl_arrs.b  ); // b
    ciErr1 |= clSetKernelArg(ckConcludeMatMult, counter++, sizeof(cl_mem), (void*)&d           );
    ciErr1 |= clSetKernelArg(ckConcludeMatMult, counter++, sizeof(cl_mem), (void*)&y           );
    ciErr1 |= clSetKernelArg(ckConcludeMatMult, counter++, sizeof(cl_mem), (void*)&u           );
    ciErr1 |= clSetKernelArg(ckConcludeMatMult, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp4);


#if 0
    { // Finally, enqueue kernel!
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckConcludeMatMult, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckConcludeMatMult;
}


void run_conclude_matmult(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckConcludeMatMult[X], 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}


/**
 * Preparing for parallel prefix with linear function composition
 * (i.e., filling in ocl_arrs.tmp)
 */
cl_kernel make_preludeFwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              d,
        cl_mem              u
) {
    cl_kernel ckPreludeFwdFunComp = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPreludeFwdFunComp = clCreateKernel(cpProgram, "prelude_fwd_fun_comp", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckPreludeFwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a  ); // b
    ciErr1 |= clSetKernelArg(ckPreludeFwdFunComp, counter++, sizeof(cl_mem), (void*)&d           );
    ciErr1 |= clSetKernelArg(ckPreludeFwdFunComp, counter++, sizeof(cl_mem), (void*)&u           );
    ciErr1 |= clSetKernelArg(ckPreludeFwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp2);


#if 0
    { // Finally, enqueue kernel!
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPreludeFwdFunComp, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPreludeFwdFunComp;
}

void run_preludeFwdMapFunComp(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPreludeFwdFunComp[X], 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}


/**
 * Mapping back the results of parallel prefix with linear function composition
 */
cl_kernel make_postFwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              y
) {
    cl_kernel ckPostFwdFunComp = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPostFwdFunComp = clCreateKernel(cpProgram, "post_fwd_fun_comp", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckPostFwdFunComp, counter++, sizeof(cl_mem), (void*)&y           );
    ciErr1 |= clSetKernelArg(ckPostFwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp2);

#if 0
    { // Finally, enqueue kernel!
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPostFwdFunComp, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPostFwdFunComp;
}

void run_postFwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPostFwdFunComp[X], 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}



/**
 * Preparing for forward-parallel prefix with linear function composition
 * (i.e., filling in ocl_arrs.tmp)
 */
cl_kernel make_preludeBwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              y,
        cl_mem              u
) {
    cl_kernel ckPreludeBwdFunComp = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPreludeBwdFunComp = clCreateKernel(cpProgram, "prelude_bwd_fun_comp", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckPreludeBwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c  ); // b
    ciErr1 |= clSetKernelArg(ckPreludeBwdFunComp, counter++, sizeof(cl_mem), (void*)&y           );
    ciErr1 |= clSetKernelArg(ckPreludeBwdFunComp, counter++, sizeof(cl_mem), (void*)&u           );
    ciErr1 |= clSetKernelArg(ckPreludeBwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp2);

#if 0
    { // Finally, enqueue kernel!
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPreludeBwdFunComp, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPreludeBwdFunComp;
}

void run_preludeBwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPreludeBwdFunComp[X], 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}

/**
 * Mapping back the results of parallel prefix with linear function composition
 */
cl_kernel make_postBwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              y
) {
    cl_kernel ckPostBwdFunComp = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckPostBwdFunComp = clCreateKernel(cpProgram, "post_bwd_fun_comp", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckPostBwdFunComp, counter++, sizeof(cl_mem), (void*)&y                );
    ciErr1 |= clSetKernelArg(ckPostBwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp2    );
    //ciErr1 |= clSetKernelArg(ckPostBwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.ro_scals);
    //ciErr1 |= clSetKernelArg(ckPostBwdFunComp, counter++, sizeof(cl_mem), (void*)&ocl_arrs.timeline);

#if 0
    { // Finally, enqueue kernel!
        ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckPostBwdFunComp, 3, NULL,
                                         globalWorkSize, localWorkSize, 0, NULL, NULL);
    }
#endif
    oclCheckError(ciErr1, CL_SUCCESS);

    return ckPostBwdFunComp;
}

void run_postBwdMapFunComp (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID MatMultScan Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckPostBwdFunComp[X], 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );

    oclCheckError(ciErr1, CL_SUCCESS);
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * Mapping back the results of parallel prefix with matrix multiplication
 */
cl_kernel make_TRIDAG_ALL_OPT (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        cl_mem              d,
        cl_mem              y,
        cl_mem              u,
        const unsigned int* NUM_P
) {
    cl_kernel ckTridagAll = NULL; // OpenCL kernel
    cl_int    ciErr1, ciErr2;

    unsigned int  counter = 0; // ro_scals

    /////
    ckTridagAll = clCreateKernel(cpProgram, "tridag_inlined", &ciErr1);
    oclCheckError(ciErr1, CL_SUCCESS);


    unsigned long int cur_size;

    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&ocl_arrs.a  ); // a
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&ocl_arrs.c  ); // c
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&d           );
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&y           );
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&u           );
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(cl_mem), (void*)&ocl_arrs.tmp4);
    ciErr1 |= clSetKernelArg(ckTridagAll, counter++, sizeof(uint)  , (void *)NUM_P        );

    {
        const int PRIV_MULT = sizeof(REAL) * 8 * WORKGROUP_SIZE;
        ciErr1 |= clSetKernelArg(ckTridagAll, counter++, PRIV_MULT, NULL);
    }

    oclCheckError(ciErr1, CL_SUCCESS);

    return ckTridagAll;
}


void run_TRIDAG_ALL_OPT (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        int                 X
) {
    cl_int  ciErr1;

    size_t  globalWorkSize = NUM_X * NUM_Y * OUTER_LOOP_COUNT;
    size_t  localWorkSize  = WORKGROUP_SIZE;

    //size_t*  globalWorkSize = (X==0) ? kernels.GWS_XYZ : kernels.GWS_YXZ;
    //size_t*  localWorkSize  = (X==0) ? kernels.LWS_XYZ : kernels.LWS_YXZ;

    assert( (X==0 || X==1) && "INVALID TRIDAG_ALL_OPT Kernel!");

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckTridagAllOpt[X], 1, NULL,
                                        &globalWorkSize, &localWorkSize, 0, NULL, NULL );

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
        const unsigned int* size_x,
        const unsigned int* size_y,
        const unsigned int* size_z,
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

#if 0
    // enqueue kernel
    ciErr1 |= clEnqueueNDRangeKernel(cqCommandQueue, ckMatTransp, 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL);
#endif
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
    oclCheckError(ciErr1, CL_SUCCESS);
}

cl_kernel make_transposeGPU_WithUpdate(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        cl_mem              inp_arr,
        cl_mem              out_arr,
        const unsigned int* size_x,
        const unsigned int* size_y,
        const unsigned int* size_z,
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
        GPUkernels&         kernels
) {
    cl_int  ciErr1;
    size_t*  globalWorkSize     = kernels.GWS_YXZ;
    size_t   localWorkSize [3]  = {BLOCK_DIM, BLOCK_DIM, 1};
    //size_t*  localWorkSize  = kernels.LWS_YXZ;

    ciErr1 = clEnqueueNDRangeKernel(cqCommandQueue, kernels.ckMatTranspUpdate, 3, NULL,
                                        globalWorkSize, localWorkSize, 0, NULL, NULL );
    oclCheckError(ciErr1, CL_SUCCESS);
}

/**********************************************************/
/**************** MAKE ALL KERNELS ************************/
/**********************************************************/

void make_TRIDAG(
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        oclNordeaArrays&    ocl_arrs,
        GPUkernels&         kernels,
        const int           X
) {
    assert( (X==0 || X==1) && "INVALID argument X to make_TRIDAG Kernels!");

    cl_mem d; cl_mem y; cl_mem u;
    const unsigned int* NUM_P = NULL;
    //size_t* globalWorkSize = NULL;
    //size_t* localWorkSize  = NULL;

    if(X == 0) {
        // first tridag run: d <- u; y <- u; u <- y NUM_P <- &NUM_X
        d = ocl_arrs.u; y = ocl_arrs.u; u = ocl_arrs.y; NUM_P = &NUM_X;
        //globalWorkSize = kernels.GWS_XYZ;
        //localWorkSize  = kernels.LWS_XYZ;
    } else { // X == 1
        // second tridag run: d <- v; y <- y; u <- u NUM_P <- &NUM_Y
        d = ocl_arrs.v; y = ocl_arrs.v; u = ocl_arrs.y; NUM_P = &NUM_Y; //y = ocl_arrs.y; u = ocl_arrs.u;
        //globalWorkSize = kernels.GWS_XYZ;
        //localWorkSize  = kernels.LWS_XYZ;
    }


    // scan with matrix multiply!
    kernels.ckTridagMatMult[X] = make_MatMultScan (
            cqCommandQueue, cpProgram, ocl_arrs, NUM_P
        );

    // map the results of the parallel prefix with matrix multiplication!
    kernels.ckConcludeMatMult[X] = make_conclude_matmult (
            cqCommandQueue, cpProgram, ocl_arrs, d, y, u
        );

    // prepare for forward scan with function composition!
    kernels.ckPreludeFwdFunComp[X] = make_preludeFwdMapFunComp (
            cqCommandQueue, cpProgram, ocl_arrs, d, u
        );

    // (forward) scan with function composition!
    kernels.ckTridagFwdFunComp[X] = make_FwdFunCompScan (
            cqCommandQueue, cpProgram, ocl_arrs, NUM_P
        );

    // post_fwd_fun_comp
    kernels.ckPostFwdFunComp[X] = make_postFwdMapFunComp (
            cqCommandQueue, cpProgram, ocl_arrs, y
        );

    // prepare for backward scan with linear-function composition
    kernels.ckPreludeBwdFunComp[X] = make_preludeBwdMapFunComp (
            cqCommandQueue, cpProgram, ocl_arrs, y, u
        );

    // (backward) scan with function composition!
    // have already been computed in "kernels.ckTridagFwdFunComp[X]" needs to be run twice!!!
    //kernels.ckTridagFwdFunComp2[X] = make_FwdFunCompScan (
    //        cqCommandQueue, cpProgram, ocl_arrs, NUM_P
    //    );

    // post_bwd_fun_comp
    kernels.ckPostBwdFunComp[X] = make_postBwdMapFunComp (
            cqCommandQueue, cpProgram, ocl_arrs, y
        );

    /////////////////////////
    // Finally, TRIDAG ALL //
    /////////////////////////
    kernels.ckTridagAllOpt[X] = make_TRIDAG_ALL_OPT (
            cqCommandQueue, cpProgram, ocl_arrs, d, y, u, NUM_P
        );
}


void make_kernels (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        RWScalars&          ro_scal,
        oclNordeaArrays&    ocl_arrs,
        GPUkernels&         kernels
) {
    initGPUgrid ( kernels );

    kernels.ckPreTridagX = make_prepare_tridag_X (
            cqCommandQueue,
            cpProgram,
            ro_scal,
            ocl_arrs,
            kernels.GWS_XYZ,
            kernels.LWS_XYZ
        );

    make_TRIDAG( cqCommandQueue, cpProgram, ocl_arrs, kernels, 0 );

    kernels.ckPreTridagY = make_prepare_tridag_Y (
                cqCommandQueue,
                cpProgram,
                ocl_arrs,
                kernels.GWS_YXZ,
                kernels.LWS_YXZ
            );

    make_TRIDAG( cqCommandQueue, cpProgram, ocl_arrs, kernels, 1 );

#if TRANSPOSE_UV
    kernels.ckMatTranspUpdate = make_transposeGPU_WithUpdate (
            cqCommandQueue, cpProgram, ocl_arrs.u, ocl_arrs.res_arr,    // ocl_arrs.y
            &NUM_Y, &NUM_X, &OUTER_LOOP_COUNT, ocl_arrs
        );
#else
    // transpose y into res_arr
    kernels.ckMatTranspUpdate = make_transposeGPU_WithUpdate (
            cqCommandQueue, cpProgram, ocl_arrs.v, ocl_arrs.res_arr,    // ocl_arrs.y
            &NUM_Y, &NUM_X, &OUTER_LOOP_COUNT, ocl_arrs
        );
#endif

    //////////////////////////////////
    // Finally, Finally, TRIDAG ALL //
    //////////////////////////////////
    kernels.ckNordeaKernelX = make_NordeaKernelX (
            cqCommandQueue, cpProgram, ro_scal, ocl_arrs
        );

    kernels.ckNordeaKernelY = make_NordeaKernelY (
            cqCommandQueue, cpProgram, ro_scal, ocl_arrs
        );

    // transpose U into Y
    kernels.ckMatTransposeU = make_transposeGPU (
            cqCommandQueue, cpProgram, ocl_arrs.u, ocl_arrs.y,
            &NUM_X, &NUM_Y, &OUTER_LOOP_COUNT, "transpose"
        );

    // transpose V into U
    kernels.ckMatTransposeV = make_transposeGPU (
            cqCommandQueue, cpProgram, ocl_arrs.v, ocl_arrs.u,
            &NUM_X, &NUM_Y, &OUTER_LOOP_COUNT, "transpose"
        );

}


/**********************************************************/
/********************* RUN KERNELS ************************/
/**********************************************************/

void run_TRIDAG(
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels,
        const int           X
) {
    assert( (X==0 || X==1) && "INVALID argument X to make_TRIDAG Kernels!");

#ifdef TRIDAG_ALL_OPT_ON
    run_TRIDAG_ALL_OPT( cqCommandQueue, kernels, X );
#else
    // scan with matrix multiply!
    run_MatMultScan( cqCommandQueue, kernels, X );

    // map the results of the parallel prefix with matrix multiplication!
    run_conclude_matmult( cqCommandQueue, kernels, X );


    // prepare for forward scan with function composition!
    run_preludeFwdMapFunComp( cqCommandQueue, kernels, X );

    // (forward) scan with function composition!
    run_FwdFunCompScan( cqCommandQueue, kernels, X);

    // post_fwd_fun_comp
    run_postFwdMapFunComp( cqCommandQueue, kernels, X );

    // prepare for backward scan with linear-function composition
    run_preludeBwdMapFunComp( cqCommandQueue, kernels, X );

    // (backward) scan with function composition!
    run_FwdFunCompScan( cqCommandQueue, kernels, X );

    // post_bwd_fun_comp
    run_postBwdMapFunComp( cqCommandQueue, kernels, X );
#endif
}

void run_GPUkernels_one_time_iteration (
        cl_command_queue&   cqCommandQueue,
        GPUkernels&         kernels
) {
#ifdef MOST_OPTIMISED_ON
    run_NordeaKernelX( cqCommandQueue, kernels );

#if TRANSPOSE_UV
    run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeU, kernels.GWS_XYZ ); // y <- transpose(u)
    run_transposeGPU ( cqCommandQueue, kernels.ckMatTransposeV, kernels.GWS_XYZ ); // u <- transpose(v)
#endif

    run_NordeaKernelY( cqCommandQueue, kernels );

#else

    run_prepare_tridag_X( cqCommandQueue, kernels );

    run_TRIDAG( cqCommandQueue, kernels, 0);

    run_prepare_tridag_Y( cqCommandQueue, kernels );

    run_TRIDAG( cqCommandQueue, kernels, 1);
#endif

    run_transposeGPU_WithUpdate( cqCommandQueue, kernels );
}

/**********************************************************/
/******************* RELEASE KERNELS **********************/
/**********************************************************/

void release_all_kernels ( GPUkernels& kernels ) {
    clReleaseKernel(kernels.ckPreTridagX);
    clReleaseKernel(kernels.ckPreTridagY);
    clReleaseKernel(kernels.ckMatTranspUpdate);

    clReleaseKernel(kernels.ckTridagMatMult[0]);
    clReleaseKernel(kernels.ckTridagMatMult[1]);

    clReleaseKernel(kernels.ckConcludeMatMult[0]);
    clReleaseKernel(kernels.ckConcludeMatMult[1]);

    clReleaseKernel(kernels.ckPreludeFwdFunComp[0]);
    clReleaseKernel(kernels.ckPreludeFwdFunComp[1]);

    clReleaseKernel(kernels.ckTridagFwdFunComp[0]);
    clReleaseKernel(kernels.ckTridagFwdFunComp[1]);

    clReleaseKernel(kernels.ckPostFwdFunComp[0]);
    clReleaseKernel(kernels.ckPostFwdFunComp[1]);

    clReleaseKernel(kernels.ckPreludeBwdFunComp[0]);
    clReleaseKernel(kernels.ckPreludeBwdFunComp[1]);

    clReleaseKernel(kernels.ckPostBwdFunComp[0]);
    clReleaseKernel(kernels.ckPostBwdFunComp[1]);

    clReleaseKernel(kernels.ckTridagAllOpt[0]);
    clReleaseKernel(kernels.ckTridagAllOpt[1]);

    clReleaseKernel(kernels.ckNordeaKernelX);
    clReleaseKernel(kernels.ckNordeaKernelY);

    clReleaseKernel(kernels.ckMatTransposeU);
    clReleaseKernel(kernels.ckMatTransposeV);
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

    clReleaseMemObject(ocl_arrs.tmp2);
    clReleaseMemObject(ocl_arrs.tmp4);
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


/**********************************************************/
/*******************  DEBUG TESTING  **********************/
/**********************************************************/

void testMatMultScan (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        cl_context&         cxGPUContext
) {
    //const unsigned NUM_X = 64, NUM_Y = 1;

    cl_mem cl_tmp; cl_int ciErr;
    REAL* tmp = new REAL[NUM_X*NUM_Y*4];

    for(int i=0; i<NUM_X*NUM_Y*4; i+=4) {
        REAL val = (REAL)(i/4 + 1);
        tmp[i]   = val;   //1.0 + val;
        tmp[i+1] = val+1; //0.0;
        tmp[i+2] = val+2; //1.0 + val;
        tmp[i+3] = val+3; //0.0;
    }

    { // 1. enqueue buffer
        unsigned long int cur_size = 4 * NUM_X * NUM_Y * sizeof(REAL);

        cl_tmp = clCreateBuffer( cxGPUContext,
                                  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                  cur_size, tmp, &ciErr
                               );
        //ciErr |= clEnqueueWriteBuffer(cqCommandQueue, cl_tmp, CL_TRUE, 0,
        //                              cur_size, tmp, 0, NULL, NULL);

        oclCheckError(ciErr, CL_SUCCESS);
    }
    { // 2. make kernel and execute
        size_t  globalWorkSize = NUM_X * NUM_Y;
        size_t  localWorkSize  = 256; //NUM_X;
        cl_kernel ckTridagMatMultX = NULL; // OpenCL kernel

        unsigned int  counter = 0; // ro_scals
        unsigned long int cur_size;

        ckTridagMatMultX = clCreateKernel(cpProgram, "scanMatMultIncl", &ciErr);
        oclCheckError(ciErr, CL_SUCCESS);

        // 1. GPU-ONLY GLOBAL ARRAYS /(4 * NUM_X * NUM_Y * OUTER_LOOP_COUNT)
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(cl_mem), (void*)&cl_tmp);


        // 2. LOCAL: 8 * sizeof(REAL) * worksize
        const int PRIV_MULT = sizeof(REAL) * 8 * localWorkSize;
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, PRIV_MULT, NULL);


        // 1. GPU-ONLY GLOBAL ARRAYS //
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(uint), (void *) (&NUM_X) );
        oclCheckError(ciErr, CL_SUCCESS);
        shrLog(stdlog, "Before executing kernel!\n");


        { // Finally, the array size:
            ciErr |= clEnqueueNDRangeKernel(cqCommandQueue, ckTridagMatMultX, 1, NULL,
                                             &globalWorkSize, &localWorkSize, 0, NULL, NULL);

            //shrLog(stdlog, "ERROR: %d, %d!\n", ciErr, CL_INVALID_WORK_GROUP_SIZE);

            oclCheckError(ciErr, CL_SUCCESS);
            shrLog(stdlog, "After executing kernel!\n");
        }

        { // read back from GPU
            ciErr = clEnqueueReadBuffer(
                            cqCommandQueue, cl_tmp, CL_TRUE, 0,
                            4 * NUM_X * NUM_Y * sizeof(REAL), tmp, 0, NULL, NULL
                    );
            oclCheckError(ciErr, CL_SUCCESS);

            clReleaseKernel(ckTridagMatMultX);
            clReleaseMemObject(cl_tmp);
            clReleaseProgram(cpProgram);

            clReleaseCommandQueue(cqCommandQueue);
            clReleaseContext(cxGPUContext);

        }
        oclCheckError(ciErr, CL_SUCCESS);
    }

    // write:
    shrLog(stdlog, "\n\nResult GPGPU SCAN MAT MULT!!!!\n\n");
    for(int j=0; j<NUM_Y; j++) {
        shrLog(stdlog, "ROW NUM %d:\n\t", j);
        REAL* tmpj = tmp + j*NUM_X*4;
        for(int i=0; i<NUM_X; i++) {
            REAL* tmpji = tmpj + i*4;
            shrLog(stdlog, "[elem(%d)=(%f, %f, %f, %f)], ",
                        i, tmpji[0], tmpji[1], tmpji[2], tmpji[3]);
        }
    }


    { // on CPU
        for(int i=0; i<NUM_X*NUM_Y*4; i+=4) {
            REAL val = (REAL)(i/4 + 1);
            tmp[i]   = val;   //1.0 + val;
            tmp[i+1] = val+1; //0.0;
            tmp[i+2] = val+2; //1.0 + val;
            tmp[i+3] = val+3; //0.0;
        }

        for(int j=0; j<NUM_Y; j++) {
            REAL* tmpj = tmp + j*NUM_X*4;
            scan_matmult_2by2(tmpj, NUM_X);
        }
        // write:
        shrLog(stdlog, "\n\nResult CPU SCAN MAT MULT!!!!\n\n");
        for(int j=0; j<NUM_Y; j++) {
            shrLog(stdlog, "ROW NUM %d:\n\t", j);
            REAL* tmpj = tmp + j*NUM_X*4;
            for(int i=0; i<NUM_X; i++) {
                REAL* tmpji = tmpj + i*4;
                shrLog(stdlog, "[elem(%d)=(%f, %f, %f, %f)], ",
                            i, tmpji[0], tmpji[1], tmpji[2], tmpji[3]);
            }
        }
    }
}


void testFunCompScan (
        cl_command_queue&   cqCommandQueue,
        cl_program          cpProgram,
        cl_context&         cxGPUContext
) {
    const unsigned NUM_X = 256, NUM_Y = 1;

    unsigned long int ARR_SIZE = 2 * NUM_X * NUM_Y * sizeof(REAL);

    cl_mem cl_tmp; cl_int ciErr;
    REAL* tmp = new REAL[NUM_X*NUM_Y*4];

    for(int i=0; i<NUM_X*NUM_Y*4; i+=2) {
        REAL val = (REAL)1.0/(i/4 + 1);
        tmp[i]   = val;   //1.0 + val;
        tmp[i+1] = val+1; //0.0;
    }

    { // 1. enqueue buffer


        cl_tmp = clCreateBuffer( cxGPUContext,
                                  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                  ARR_SIZE, tmp, &ciErr
                               );
        //ciErr |= clEnqueueWriteBuffer(cqCommandQueue, cl_tmp, CL_TRUE, 0,
        //                              cur_size, tmp, 0, NULL, NULL);
        oclCheckError(ciErr, CL_SUCCESS);
    }
    { // 2. make kernel and execute
        size_t  globalWorkSize = NUM_X * NUM_Y;
        size_t  localWorkSize  = NUM_X;
        cl_kernel ckTridagMatMultX = NULL; // OpenCL kernel

        unsigned int  counter = 0; // ro_scals
        unsigned long int cur_size;

        ckTridagMatMultX = clCreateKernel(cpProgram, "scanLinFunCompIncl", &ciErr);
        oclCheckError(ciErr, CL_SUCCESS);

        // 1. GPU-ONLY GLOBAL ARRAYS /(4 * NUM_X * NUM_Y * OUTER_LOOP_COUNT)
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(cl_mem), (void*)&cl_tmp);


        // 2. LOCAL: 8 * sizeof(REAL) * worksize
        const int PRIV_MULT = sizeof(REAL) * 8 * localWorkSize;
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, PRIV_MULT, NULL);


        // 1. GPU-ONLY GLOBAL ARRAYS //
        ciErr |= clSetKernelArg(ckTridagMatMultX, counter++, sizeof(uint), (void *) (&NUM_X) );
        oclCheckError(ciErr, CL_SUCCESS);
        shrLog(stdlog, "Before executing kernel!\n");


        { // Finally, the array size:
            ciErr |= clEnqueueNDRangeKernel(cqCommandQueue, ckTridagMatMultX, 1, NULL,
                                             &globalWorkSize, &localWorkSize, 0, NULL, NULL);

            //shrLog(stdlog, "ERROR: %d, %d!\n", ciErr, CL_INVALID_WORK_GROUP_SIZE);

            oclCheckError(ciErr, CL_SUCCESS);
            shrLog(stdlog, "After executing kernel!\n");
        }

        { // read back from GPU
            ciErr = clEnqueueReadBuffer(
                            cqCommandQueue, cl_tmp, CL_TRUE, 0,
                            ARR_SIZE, tmp, 0, NULL, NULL
                    );
            oclCheckError(ciErr, CL_SUCCESS);

            clReleaseKernel(ckTridagMatMultX);
            clReleaseMemObject(cl_tmp);
            clReleaseProgram(cpProgram);

            clReleaseCommandQueue(cqCommandQueue);
            clReleaseContext(cxGPUContext);

        }
        oclCheckError(ciErr, CL_SUCCESS);
    }

    // write:
    shrLog(stdlog, "\n\nResult GPGPU SCAN MAT MULT!!!!\n\n");
    for(int j=0; j<NUM_Y; j++) {
        shrLog(stdlog, "ROW NUM %d:\n\t", j);
        REAL* tmpj = tmp + j*NUM_X*2;
        for(int i=0; i<NUM_X; i++) {
            REAL* tmpji = tmpj + i*2;
            shrLog(stdlog, "[elem(%d)=(%f, %f)], ", i, tmpji[0], tmpji[1]);
        }
    }


    { // on CPU
        for(int i=0; i<NUM_X*NUM_Y*4; i+=2) {
            REAL val = (REAL)1.0/(i/4 + 1);
            tmp[i]   = val;   //1.0 + val;
            tmp[i+1] = val+1; //0.0;
        }

        for(int j=0; j<NUM_Y; j++) {
            REAL* tmpj = tmp + j*NUM_X*4;
            scan_linfuncomp(tmpj, NUM_X);
        }

        // write:
        shrLog(stdlog, "\n\nResult CPU SCAN MAT MULT!!!!\n\n");
        for(int j=0; j<NUM_Y; j++) {
            shrLog(stdlog, "ROW NUM %d:\n\t", j);
            REAL* tmpj = tmp + j*NUM_X*2;
            for(int i=0; i<NUM_X; i++) {
                REAL* tmpji = tmpj + i*2;
                shrLog(stdlog, "[elem(%d)=(%f, %f)], ", i, tmpji[0], tmpji[1]);
            }
        }
    }
}

