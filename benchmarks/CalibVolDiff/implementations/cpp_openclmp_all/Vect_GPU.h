#ifndef NORDEA_GPU
#define NORDEA_GPU

#include "SDK_stub.h"

#include "PrepareKernels.h"

////////////////////////////////////////////////////////////
//// OpenCL internals!
////////////////////////////////////////////////////////////

void verifyEnoughResources(cl_device_id device) {
    cl_ulong glob_mem_size;
    clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(glob_mem_size),&glob_mem_size,NULL);
    { // check global mem matches the user-declared one
        const size_t USER_MEM = ((size_t)GPU_GLB_MEM) << 30;
        const size_t ub       = static_cast<size_t>( USER_MEM * 101.0 / 100.0 );
        const size_t lb       = static_cast<size_t>( USER_MEM *  99.0 / 100.0 );
        if ( !(glob_mem_size >= lb && glob_mem_size <= ub) ) {
            fprintf(stderr, "WARNING! Querried GPU global memory %lu DIFFERS from what user declared: [%ld,%ld]!\n", 
                             glob_mem_size, lb, ub);
            glob_mem_size = static_cast<cl_ulong>(USER_MEM);
            fprintf(stderr, "Using user-declared size for global memory: %lu\n", glob_mem_size);
        }
    }

    cl_ulong glob_mem_required = 15 * OUTER_LOOP_COUNT * NUM_X * NUM_Y * sizeof(REAL);
    assert( glob_mem_required <= glob_mem_size && "Not Enough Global Memory. A possible fix is to strip-mine all kernels.\n" );
    //fprintf(stderr, "GPU global memory: %lu, needed: %lu bytes\n", glob_mem_size, glob_mem_required);
}


void makeOclBuffers (
        cl_context          cxGPUContext,
        cl_command_queue&   cqCommandQueue,
        RWScalars&          ro_scal,
        NordeaArrays&       cpu_arrs,
        oclNordeaArrays&    ocl_arrs
) {
    cl_int ciErr;
    cl_int ciErr2;
    unsigned long int cur_size;

    /*********************/
    /*** 1. RO scalars ***/
    /*********************/
    {
        cur_size = sizeof(RWScalars);
        ocl_arrs.ro_scals = clCreateBuffer(
                    cxGPUContext, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2
                );

        ciErr = ciErr2;
        ciErr |= clEnqueueWriteBuffer(cqCommandQueue, ocl_arrs.ro_scals, CL_TRUE, 0,
                                      cur_size, &ro_scal, 0, NULL, NULL);
        //oclCheckError(ciErr, CL_SUCCESS);
    }

    /*********************/
    /*** 1. RO Arrays! ***/
    /*********************/
    {
        cur_size = NUM_T*sizeof(REAL);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(timeline);

        cur_size = NUM_X*sizeof(REAL);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myX);

        cur_size *= REAL3_CT;
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myDx);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myDxx);

        cur_size = NUM_Y*sizeof(REAL);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myY);

        cur_size *= REAL3_CT;
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myDy);
        CREATE_AND_ENQUEUE_RO_CL_BUFFER(myDyy);
    }

    //oclCheckError(ciErr, CL_SUCCESS);

    /*********************/
    /*** 1. RW Arrays! ***/
    /*********************/
    {
        cur_size = OUTER_LOOP_COUNT * NUM_X * NUM_Y * sizeof(REAL);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(a);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(b);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(c);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(y);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(u);
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(v);

        cur_size *= 2;
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(tmp2);

        cur_size *= 2;
        CREATE_AND_ENQUEUE_GPU_ONLY_CL_BUFFER(tmp4);
    }

    { // Finally the result array (i.e., Read-Write but needs to be brought back to CPU in the end)
        cur_size = OUTER_LOOP_COUNT * NUM_X * NUM_Y * sizeof(REAL);
        CREATE_AND_ENQUEUE_OUT_CL_BUFFER(res_arr);
    }

    oclCheckError(ciErr, CL_SUCCESS);

    //shrLog("Before EXITING ALLOC ARRAYS!!! ...\n");

}

#define USE_GPU_ITER 1

//extern "C"
unsigned long int
runOnGPU ( RWScalars& ro_scal, NordeaArrays& cpu_arrs, oclNordeaArrays& ocl_arrs ) {
    cl_context cxGPUContext;                        // OpenCL context
    cl_command_queue cqCommandQueue[16];            // OpenCL command que
    cl_uint nDevice;                                // OpenCL device count
    cl_device_id* cdDevices;                        // OpenCL device list
    cl_program cpProgram;                           // OpenCL program
    GPUkernels  kernels;                            // OpenCL kernel
    const unsigned int dev_id = GPU_DEV_ID;

    assert(GPU_DEV_ID >= 0 && "GPU DEVICE ID < 0 !\n");

    { // initialize the loop-variant scalars!
        ro_scal.t_ind      = NUM_T - 2;
        ro_scal.dtInv      = 1.0/(cpu_arrs.timeline[ro_scal.t_ind+1]-cpu_arrs.timeline[ro_scal.t_ind]);
        ro_scal.timeline_i = cpu_arrs.timeline[ro_scal.t_ind];
    }

    { // making command queue, building program, etc
        char  compile_opts[128];
        sprintf( compile_opts, "-D lgWARP=%d", lgWARP );

        build_for_GPU(
            cxGPUContext, cqCommandQueue,
            nDevice, cdDevices, cpProgram, dev_id, compile_opts, "", "CrankNicolson"
        );

        verifyEnoughResources(cdDevices[dev_id]);
    }

    // allocate space for the RO and RW arrays on GPU!
    makeOclBuffers (
            cxGPUContext, cqCommandQueue[dev_id], ro_scal, cpu_arrs, ocl_arrs
        );


    { // make all kernels!
        make_kernels( cqCommandQueue[dev_id], cpProgram, ro_scal, ocl_arrs, kernels );
    }

    unsigned long int elapsed;
    struct timeval t_start, t_end, t_diff;
    gettimeofday(&t_start, NULL);
    
    // now execute kernels and record the time!
    for(int t_ind = NUM_T-2; t_ind>=0; --t_ind) {
            run_GPUkernels_one_time_iteration ( cqCommandQueue[dev_id], kernels );
    } // END TIME LOOP!


    { // WRITE BACK THE RESULT ARRAY TO CPU !!! //
        cl_int  ciErr;
        const unsigned int ARR_SIZE = NUM_X * NUM_Y * OUTER_LOOP_COUNT * sizeof(REAL);
        ciErr = clEnqueueReadBuffer (
                            cqCommandQueue[dev_id], ocl_arrs.res_arr, CL_TRUE,
                            0, ARR_SIZE, cpu_arrs.res_arr, 0, NULL, NULL
                );
        oclCheckError(ciErr, CL_SUCCESS);

        // RELEASE ALL GPU RESOURCES
        release_all_GPU_resources (
                cqCommandQueue[dev_id], cxGPUContext, cpProgram, cdDevices, ocl_arrs, kernels
            );
    }

    gettimeofday(&t_end, NULL);
    timeval_subtract(&t_diff, &t_end, &t_start);
    elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;

    return elapsed;
}

#endif // end include NORDEA_GPU
