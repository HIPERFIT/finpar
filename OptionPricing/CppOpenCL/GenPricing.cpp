#include "ParseInput.h"
#include "StructGPU.h"
//#include "Contracts.h"
//#include <omp.h>

/*****************************************************/
/*****************************************************/
/*****************************************************/

double*
run_GPUkernel (
        LoopROScalars   &   ro_scal,
        SobolArrays     &   sob_arrs,
        ModelArrays     &   md_arrs,
        BrowBridgeArrays&   bb_arrs,
        int             &   num_cores,
        size_t          &   elapsed          
        

) {
    cl_context cxGPUContext;                        // OpenCL context
    cl_command_queue cqCommandQueue[16]; 			// OpenCL command que
    cl_uint nDevice;                                // OpenCL device count
    cl_device_id* cdDevices;                        // OpenCL device list
    cl_program cpProgram;                           // OpenCL program
    cl_int ciErr1;                                  // Error code var

    size_t globalWorkSize[1]; 	                    // 1D var for Total # of work items
    size_t localWorkSize [1]; 	                    // 1D var for # of work items in the work group, e.g., { 128 }

    REAL* glb_vhat = NULL;


    double* vhat_fin = (double*)malloc(ro_scal.num_models*sizeof(double));
    for( UINT i = 0; i < ro_scal.num_models; i++ ) vhat_fin[i] = 0.0;

    GPU_KERNEL  kernel_type = priv_or_vect_kernel( ro_scal );
#ifdef _OPTIMIZATION_MEM_COALES_ON
    const char* kernel_name = (kernel_type == PRIV) ?
        		            "GenericPricingPrivOpt" :   // .cl
                                    "GenericPricingVectOpt" ;   // .cl
#else
    const char* kernel_name = "GenericPricingVectUncoalesced";
#endif

    { // build the OpenCL program
        const char* preamble =  makeGPUprogPreamble( ro_scal, kernel_type );
        char  compile_opts[128];
        sprintf( compile_opts, "-D lgWARP=%d -D TILE=%d", lgWARP, TILE );
        build_for_GPU(  cxGPUContext, cqCommandQueue, nDevice, 
                        cdDevices,    cpProgram,      GPU_DEV_ID,
                        compile_opts, preamble,       kernel_name  );
    }
    // we do not measure the just-in-time compilation time!
    struct timeval t_start, t_end, t_diff;
    gettimeofday(&t_start, NULL);

    {
        num_cores = discriminate_cost_model(ro_scal, cdDevices[GPU_DEV_ID], kernel_type);
        computeSobolFixIndex(sob_arrs, ro_scal.chunk);

        globalWorkSize[0]    = getWorkSize( ro_scal.num_gpuits, ro_scal.chunk, ro_scal.BLOCK );
        localWorkSize [0]    = ro_scal.BLOCK;
    }

    oclLoopArrays ocl_arrs;
    if(kernel_type == PRIV) { // CALL THE PRIVATE KERNEL
        glb_vhat = oclAllocArrays_PrivKernel(
                        ocl_arrs, cxGPUContext, cqCommandQueue[GPU_DEV_ID],
                        ciErr1, ro_scal, sob_arrs, bb_arrs, md_arrs
                    );

        runGPU_PRIV (
                ro_scal, glb_vhat, ocl_arrs, cqCommandQueue[GPU_DEV_ID], cpProgram, 
                globalWorkSize, localWorkSize
            );

        reduceVHAT_CPU( ro_scal, glb_vhat, vhat_fin );

    } else { // VECT VERSION!

        glb_vhat = oclAllocArrays_VectKernel(
                    ocl_arrs, cxGPUContext, cqCommandQueue[GPU_DEV_ID],
                    ciErr1, ro_scal, sob_arrs, bb_arrs, md_arrs
                );

        UINT sob_ini_count = ro_scal.sobol_count_ini;

        for(UINT cur_iter = 0; cur_iter < ro_scal.num_mcits; cur_iter+=ro_scal.num_gpuits) {
            if(cur_iter + ro_scal.num_gpuits > ro_scal.num_mcits) {
                ro_scal.num_gpuits = ro_scal.num_mcits - cur_iter;
                globalWorkSize[0] = getWorkSize( ro_scal.num_gpuits,ro_scal.chunk,ro_scal.BLOCK);
            }
            if(cur_iter != 0) {
                size_t cur_size = sizeof(LoopROScalars);
                cl_int ciErr2;
                ro_scal.sobol_count_ini = cur_iter + sob_ini_count;
                clReleaseMemObject(ocl_arrs.ro_scals);
                ocl_arrs.ro_scals = clCreateBuffer(
                                        cxGPUContext, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2
                                    );
                ciErr2 = clEnqueueWriteBuffer(cqCommandQueue[GPU_DEV_ID], ocl_arrs.ro_scals, CL_TRUE, 0,
                                                cur_size, &ro_scal, 0, NULL, NULL);
                oclCheckError(ciErr2, CL_SUCCESS);
            }

            runGPU_VECT (
                    ro_scal, glb_vhat, ocl_arrs, cqCommandQueue[GPU_DEV_ID], cpProgram,
                    globalWorkSize, localWorkSize
                );

            reduceVHAT_CPU( ro_scal, glb_vhat, vhat_fin );
        }

        ro_scal.sobol_count_ini = sob_ini_count;
    }

    for(UINT ii = 0; ii<ro_scal.num_models; ii++) {
        vhat_fin[ii] = vhat_fin[ii] / ro_scal.num_mcits;
    }

    gettimeofday(&t_end, NULL);
    timeval_subtract(&t_diff, &t_end, &t_start);
    elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;

    { // free allocated space
        shrLog(stderr, "Release CPU buffers and OpenCL objects...\n");
        // clFinish(cqCommandQueue[0]);
        clReleaseProgram(cpProgram);
        if(kernel_type == PRIV) ocl_arrs.cleanupPRIV();
        else                    ocl_arrs.cleanupVECT();
        free(glb_vhat);
        clReleaseCommandQueue(cqCommandQueue[GPU_DEV_ID]);
        free(cdDevices);
        clReleaseContext(cxGPUContext);
    }

    return vhat_fin;
}




/*****************************************************/
/*****************************************************/
/*****************************************************/
int main() {
    LoopROScalars    scals;
    SobolArrays      sob_arrs;
    ModelArrays      md_arrs;
    BrowBridgeArrays bb_arrs;

    fprintf(stdout, "\n// Generic Pricing, Multi-Threaded Benchmark:\n");

    readDataSet(scals, sob_arrs, md_arrs, bb_arrs);

    fprintf(stdout, "// Contract: %d, MC Its#: %d, #Underlyings: %d, #Path Dates: %d, chunk: %d\n\n", 
            scals.contract, scals.num_mcits, scals.num_under, scals.num_dates, scals.chunk      );

    int     num_cores;
    double* prices;
    unsigned long int elapsed;
    { // run kernel
        prices = run_GPUkernel( scals, sob_arrs, md_arrs, bb_arrs, num_cores, elapsed );

        sob_arrs.cleanup();
         md_arrs.cleanup();
         bb_arrs.cleanup();
    }

    {   // validation and writeback of the result
        bool is_valid = validate   ( scals.num_models, prices );
        writeStatsAndResult( is_valid, scals.num_models, prices, true, num_cores, elapsed ); 
        free(prices);       
    }
}


// cat ../Data/Medium/input.data ../Data/Medium/output.data  | ./GenPricing 2> Debug.txt
