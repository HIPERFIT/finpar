#ifndef NORDEA_GPU
#define NORDEA_GPU

#include "../../include/SDK_stub.h"

#include "PrepareKernels.h"

////////////////////////////////////////////////////////////
//// OpenCL internals!
////////////////////////////////////////////////////////////

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
void runOnGPU ( RWScalars& ro_scal, NordeaArrays& cpu_arrs, oclNordeaArrays& ocl_arrs ) {
    cl_context cxGPUContext;                        // OpenCL context
    cl_command_queue cqCommandQueue[16];            // OpenCL command que
    cl_uint nDevice;                                // OpenCL device count
    cl_device_id* cdDevices;                        // OpenCL device list
    cl_program cpProgram;                           // OpenCL program
    GPUkernels  kernels;                            // OpenCL kernel
    unsigned int dev_id = 0;

    { // initialize the loop-variant scalars!
        ro_scal.t_ind      = NUM_T - 2;
        ro_scal.dtInv      = 1.0/(cpu_arrs.timeline[ro_scal.t_ind+1]-cpu_arrs.timeline[ro_scal.t_ind]);
        ro_scal.timeline_i = cpu_arrs.timeline[ro_scal.t_ind];
    }

    // making command queue, building program, etc
    build_for_GPU(
            cxGPUContext, cqCommandQueue,
            nDevice, cdDevices, cpProgram, dev_id, NULL, "CrankNicolson"
        );


    // allocate space for the RO and RW arrays on GPU!
    makeOclBuffers (
            cxGPUContext, cqCommandQueue[dev_id], ro_scal, cpu_arrs, ocl_arrs
        );


    { // make all kernels!
        make_kernels( cqCommandQueue[dev_id], cpProgram, ro_scal, ocl_arrs, kernels );
    }

    // now execute kernels!
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


}

#endif // end include NORDEA_GPU
