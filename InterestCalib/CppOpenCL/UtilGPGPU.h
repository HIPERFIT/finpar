#ifndef GPGPU_UTILITIES
#define GPGPU_UTILITIES

#include "SDK_stub.h"

/************************/
/*** Helper Functions ***/
/************************/
bool is_pow2(uint n) {
    uint c = 1;
    while(c < n) {
        c = c << 1;
    }
    return (c == n);
}


struct CpuArrays {
    // size of the shape array
    const uint SS;

    // shape of the irregular arrays for one genome
    short* shape;

    // [13 * POP_SIZE] = 
    // { a, a_p, b, b_p, rho, rho_p, nu, nu_p, 
    //   sigma, sigma_p, logLik, logLik_p, bf_rat }
    REAL* genomes; 

    // [ 4 * SS * POP_SIZE ]
    REAL* ci_t1cs_scale;  

    // [ 2 * NUM_SWAP_QUOTES * POP_SIZE ]
    REAL* new_quote_price;

    // [ NUM_SWAP_QUOTES*NUM_HERMITE ]
    REAL* accum0;

    // [ 8 * NUM_SWAP_QUOTES * POP_SIZE ]
    // { mux, muy, sqrt_sigmax = sqrt(2.0) * sigmax, 
    //   t2 = rhoxy / (sigmax*rhoxycs), sigmay_rhoxycs, zc_mat, f, df } 
    REAL* scalars; 

    REAL  gene_ranges[10];

    CpuArrays (const uint n, short* shp) : SS(n) {
        shape           = shp;
    
        genomes         = new REAL[ 13 * POP_SIZE ];
        ci_t1cs_scale   = new REAL[  4 * SS * POP_SIZE ];
        scalars         = new REAL[  8 * NUM_SWAP_QUOTES * POP_SIZE ];

        new_quote_price = new REAL [ 2 * NUM_SWAP_QUOTES * POP_SIZE ];
        accum0          = new REAL [ NUM_SWAP_QUOTES * NUM_HERMITE * POP_SIZE ];

        for(int i=0; i<5; i++) { gene_ranges[i]   = g_mins[i]; }
        for(int i=0; i<5; i++) { gene_ranges[i+5] = g_maxs[i]; }
    }

    void releaseResources() {
        delete[] shape;
        delete[] genomes;
        delete[] ci_t1cs_scale;
        delete[] scalars;
        delete[] new_quote_price;
        delete[] accum0;
    }

    // genome helpers
    REAL* get_a     () { return genomes; }
    REAL* get_b     () { return genomes + POP_SIZE*2; }
    REAL* get_rho   () { return genomes + POP_SIZE*4; }
    REAL* get_nu    () { return genomes + POP_SIZE*6; }
    REAL* get_sigma () { return genomes + POP_SIZE*8; }
    REAL* get_logLik() { return genomes + POP_SIZE*10;}
    REAL* get_bf_rat() { return genomes + POP_SIZE*12;}

    // get the shape of the irregular array (for one genome)
    short* get_shape (){ return shape; }
    
    // get the start iterator into arrays ci, t1cs, scale for genome i
    REAL* get_ci    (int i) { return ci_t1cs_scale + SS*i;                 }
    REAL* get_t1cs  (int i) { return ci_t1cs_scale + SS*i +   SS*POP_SIZE; }
    REAL* get_scale (int i) { return ci_t1cs_scale + SS*i + 2*SS*POP_SIZE; }
    REAL* get_bbi   (int i) { return ci_t1cs_scale + SS*i + 3*SS*POP_SIZE; }

    // get the start iterator into array-expanded scalars for genome i
    REAL* get_mux   (int i) { return scalars + i * NUM_SWAP_QUOTES ; }
    REAL* get_muy   (int i) { return scalars + i * NUM_SWAP_QUOTES +     NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_sqsigx(int i) { return scalars + i * NUM_SWAP_QUOTES + 2 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_t2    (int i) { return scalars + i * NUM_SWAP_QUOTES + 3 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_sigrho(int i) { return scalars + i * NUM_SWAP_QUOTES + 4 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_zcmat (int i) { return scalars + i * NUM_SWAP_QUOTES + 5 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_f     (int i) { return scalars + i * NUM_SWAP_QUOTES + 6 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_df    (int i) { return scalars + i * NUM_SWAP_QUOTES + 7 * NUM_SWAP_QUOTES * POP_SIZE; }

    // get the quotes prices for genome i
    REAL* get_quote (int i) { return new_quote_price + i * NUM_SWAP_QUOTES; }
    REAL* get_price (int i) { return new_quote_price + i * NUM_SWAP_QUOTES + NUM_SWAP_QUOTES * POP_SIZE; }

    // get the accumulator for the irregular array for genome i 
    REAL* get_accum (int i) { return accum0 + i * NUM_SWAP_QUOTES * NUM_HERMITE; }
};

/**************************************/
/*** OpenCL related Data Structures ***/
/**************************************/
struct OclBuffers {
    // Used for BestFit reduction (size LWG_FB)
    cl_mem best_ind;
    cl_mem best_val;

    cl_mem    shape;           // [ SS * 4 ] (short)
    cl_mem    swap_quotes;     // [ NUM_SWAP_QUOTES * 4 ] (real)
    cl_mem    genomes;         // [13 * POP_SIZE] (real)
    cl_mem    ci_t1cs_scale;   // [ 4 * SS * POP_SIZE ] (real)
    cl_mem    new_quote_price; // [ 2 * NUM_SWAP_QUOTES * POP_SIZE ] (real)
    cl_mem    accum0;          // [ NUM_SWAP_QUOTES * NUM_HERMITE ]  (real)
    cl_mem    scalars;         // [ 8 * NUM_SWAP_QUOTES * POP_SIZE ] (real)

    cl_mem    gauss_coefs;    // NUM_HERMITE * sizeof(REAL);
    cl_mem    gauss_weights;  // NUM_HERMITE * sizeof(REAL);

    cl_mem    gene_ranges;    // 10 * sizeof(REAL)
    cl_mem    sobol_dir_vct;  // 30 * sizeof(REAL)

    cl_mem    accept_cond;    // POP_SIZE*sizeof(char)


    void releaseResources (  ) {
        shrLog(stdlog, "Release oclBuffers ... ");

        clReleaseMemObject( best_ind        );
        clReleaseMemObject( best_val        );

        clReleaseMemObject( shape           );
        clReleaseMemObject( swap_quotes     );

        clReleaseMemObject( genomes         );
        clReleaseMemObject( ci_t1cs_scale   );
        clReleaseMemObject( new_quote_price );
        clReleaseMemObject( accum0          );
        clReleaseMemObject( scalars         );

        clReleaseMemObject( gauss_coefs     );
        clReleaseMemObject( gauss_weights   );

        clReleaseMemObject( gene_ranges     );
        clReleaseMemObject( sobol_dir_vct   );

        clReleaseMemObject( accept_cond     );
    }
};

struct OclObjects {
    UINT              dev_id;

    cl_command_queue  cmdQueue[16];  // OpenCL command queue    
    cl_context        context;       // OpenCL context
    cl_device_id*     devices;       // OpenCL device list
    cl_program        program;       // OpenCL program

    cl_command_queue& getCommandQueue() { 
        return cmdQueue[dev_id]; 
    }
    void releaseResources (  ) {
        shrLog(stdlog, "Release OpenCL objects ...");
        clReleaseProgram     ( program );
        clReleaseCommandQueue( getCommandQueue() );
        free(devices);
        clReleaseContext( context );
    }
};


/******************************************/
/*** Filling the OpenCL data structures ***/
/******************************************/
void compileGPUprog( OclObjects& objs ) {
    char compile_option[512];

#ifdef MAC
    // Build the program with 'mad' Optimization option
    char flags[] = "-cl-mad-enable -DMAC ";
#else
    char flags[] = "";
#endif

    //sprintf( compile_option, "-I. -I SrcCL -D TODAY=%10f -D LWG_FB=%d -D INFTY=%f", TODAY, LWG_FB, INFTY );
    sprintf( compile_option, "%s -I. -I SrcCL -I../includeC -D lgWARP=%d -D TODAY=%10f -D LWG_FB=%d -D INFTY=%f -D MAX_DATE=%f -D MIN_DATE=%f -D SOBOL_BITS_NUM=%d",
	     flags, lgWARP, TODAY, LWG_FB, INFTY, MAX_DATE, MIN_DATE, NUM_SOBOL_BITS );

    objs.dev_id = 0;
    cl_uint nDevice;
    build_for_GPU(
        objs.context, 
        objs.cmdQueue,
        nDevice, 
        objs.devices, 
        objs.program, 
        objs.dev_id, 
        compile_option, 
        "",
        "SrcCL/CalibKers"
    );
}

void makeOclBuffers ( CpuArrays& cpu_arrs, OclObjects& ocl_objs, OclBuffers& ocl_arrs ) {
    cl_int ciErr, ciErr2;
    size_t cur_size;
    cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();

    cur_size = LWG_FB * sizeof(REAL);
    ocl_arrs.best_val = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2); ciErr |= ciErr2;

    cur_size = LWG_FB * sizeof(uint);
    ocl_arrs.best_ind = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2); ciErr |= ciErr2;

    // for the main kernel(s)
    cur_size = 4 * cpu_arrs.SS * sizeof(short);
    ocl_arrs.shape = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.shape, CL_TRUE, 0, cur_size, cpu_arrs.shape, 0, NULL, NULL);

    // swaption quotes
    cur_size = 4 * NUM_SWAP_QUOTES * sizeof(REAL);
    ocl_arrs.swap_quotes = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.swap_quotes, CL_TRUE, 0, cur_size, SwaptionQuotes, 0, NULL, NULL);

    // gene ranges
    cur_size = 10 * sizeof(REAL);
    ocl_arrs.gene_ranges = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.gene_ranges, CL_TRUE, 0, cur_size, cpu_arrs.gene_ranges, 0, NULL, NULL);

    // sobol direction vector
    cur_size = 30 * sizeof(REAL);
    ocl_arrs.sobol_dir_vct = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.sobol_dir_vct, CL_TRUE, 0, cur_size, SobolDirVct, 0, NULL, NULL);

    // acceptance condition
    cur_size = POP_SIZE * sizeof(unsigned char);
    ocl_arrs.accept_cond = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
    
    // GaussHermite coefficients and weights!
    cur_size = NUM_HERMITE * sizeof(REAL);
    ocl_arrs.gauss_coefs = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.gauss_coefs, CL_TRUE, 0, cur_size, HermiteCoeffs, 0, NULL, NULL);
    ocl_arrs.gauss_weights = clCreateBuffer( ocl_objs.context, CL_MEM_READ_ONLY, cur_size, NULL, &ciErr2 );    ciErr |= ciErr2;
    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.gauss_weights, CL_TRUE, 0, cur_size, HermiteWeights, 0, NULL, NULL);

    cur_size = 13 * POP_SIZE * sizeof(REAL);
    ocl_arrs.genomes = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
//    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.genomes, CL_TRUE, 0, cur_size, cpu_arrs.genomes, 0, NULL, NULL);

    cur_size = 4 * cpu_arrs.SS * POP_SIZE * sizeof(REAL);
    ocl_arrs.ci_t1cs_scale = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
//    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.ci_t1cs_scale, CL_TRUE, 0, cur_size, cpu_arrs.ci_t1cs_scale, 0, NULL, NULL);

    cur_size = 2 * NUM_SWAP_QUOTES * POP_SIZE * sizeof(REAL);
    ocl_arrs.new_quote_price = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
//    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.new_quote_price, CL_TRUE, 0, cur_size, cpu_arrs.new_quote_price, 0, NULL, NULL);

    cur_size = 8 * NUM_SWAP_QUOTES * POP_SIZE * sizeof(REAL);
    ocl_arrs.scalars = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
//    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.scalars, CL_TRUE, 0, cur_size, cpu_arrs.scalars, 0, NULL, NULL);

    cur_size = NUM_SWAP_QUOTES * NUM_HERMITE * POP_SIZE * sizeof(REAL);
    ocl_arrs.accum0 = clCreateBuffer( ocl_objs.context, CL_MEM_READ_WRITE, cur_size, NULL, &ciErr2 ); ciErr |= ciErr2;
//    ciErr |= clEnqueueWriteBuffer(cmd_queue, ocl_arrs.accum0, CL_TRUE, 0, cur_size, cpu_arrs.accum0, 0, NULL, NULL);

    oclCheckError(ciErr, CL_SUCCESS);
}

void initGPUresources ( CpuArrays& cpu_arrs, OclObjects& ocl_objs, OclBuffers& ocl_arrs ) {
    compileGPUprog ( ocl_objs );
    makeOclBuffers ( cpu_arrs, ocl_objs, ocl_arrs );
}

void releaseGPUresources( OclObjects& ocl_objs, OclBuffers& ocl_arrs ) {
    fprintf(stderr, "Releasing GPU Resources ... ");
    ocl_arrs.releaseResources();
    ocl_objs.releaseResources();
    fprintf(stderr, "OK!\n");
}

#endif // GPGPU_UTILITIES
