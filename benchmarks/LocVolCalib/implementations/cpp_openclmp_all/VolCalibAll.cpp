#include "Util.h"
#include "Constants.h"
#include "DataStructConst.h"
#include "ParseInput.h"
#include "VolCalibInit.h"
#include "Vect_CPU.h"
#include "Vect_GPU.h"


unsigned long int 
whole_loop_nest (
        real_t*       res,
        const real_t* strikes,
        const real_t  s0,
        const real_t  t,
        const real_t  alpha,
        const real_t  nu,
        const real_t  beta
) {
    // loop index
    unsigned int i;
    unsigned long int  elapsed;
    //const unsigned int NUM_XY = NUM_X*NUM_Y;

    // arrays:
    real_t *a = NULL, *b = NULL, *c = NULL,
         *y = NULL, *u = NULL, *v = NULL,
         *scan_tmp = NULL;

    { // allocate arrays
        a = new real_t __attribute__ ((aligned (32))) [OUTER_LOOP_COUNT*NUM_XY];
        b = new real_t[OUTER_LOOP_COUNT*NUM_XY];
        c = new real_t[OUTER_LOOP_COUNT*NUM_XY];

        y = new real_t[OUTER_LOOP_COUNT*NUM_XY];
        u = new real_t[OUTER_LOOP_COUNT*NUM_XY];
        v = new real_t[OUTER_LOOP_COUNT*NUM_XY];

        scan_tmp = new real_t[4*OUTER_LOOP_COUNT*NUM_XY];
    }

    // parallel!
    whole_loop_nest_init ( s0, t, alpha, nu, beta );

    // parallel!
    for( i=0; i<OUTER_LOOP_COUNT; ++i ) {
        real_t* res_arr = myResArr+i*NUM_X*NUM_Y;
        setPayoff_expanded(strikes[i], res_arr);
    }

    if( IS_GPU > 0 ) {
        RWScalars         ro_scal;
        NordeaArrays      cpu_arrs;
        oclNordeaArrays   ocl_arrs;

        { // init arrays
            cpu_arrs.myX = myX; cpu_arrs.myDx = myDx; cpu_arrs.myDxx = myDxx;
            cpu_arrs.myY = myY; cpu_arrs.myDy = myDy; cpu_arrs.myDyy = myDyy;
            cpu_arrs.timeline = myTimeline;
            cpu_arrs.a = a; cpu_arrs.b = b; cpu_arrs.c = c;
            cpu_arrs.y = y; cpu_arrs.u = u; cpu_arrs.v = v;
            cpu_arrs.tmp = scan_tmp;
            cpu_arrs.res_arr  = &myResArr[0];
        }

        { // init scalars
            ro_scal.NUM_X = NUM_X; ro_scal.NUM_Y = NUM_Y; ro_scal.NUM_XY = NUM_X * NUM_Y;
            ro_scal.alpha = alpha; ro_scal.beta  = beta;  ro_scal.nu     = nu;
        }

        { // SAFETY CHECK!
            bool is_safe =  (NUM_X <= WORKGROUP_SIZE) && (NUM_Y <= WORKGROUP_SIZE) &&
                            is_pow2(NUM_X) && is_pow2(NUM_Y) &&
                            (WORKGROUP_SIZE % NUM_X == 0) && (WORKGROUP_SIZE % NUM_Y == 0);
            assert(is_safe && "NOT SAFE TO PARALLELISE ON GPU!");
        }

        elapsed = runOnGPU ( ro_scal, cpu_arrs, ocl_arrs );
    } else {
        // CPU instrumentation here: 

        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        for(int t_ind = NUM_T-2; t_ind>=0; --t_ind) {

            iteration_expanded_CPU (
                    t_ind, alpha, beta, nu,
                    a, b, c, y, u, v, scan_tmp
                );
        }

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    for( i=0; i<OUTER_LOOP_COUNT; ++i ) {
        real_t* res_arr = myResArr+i*NUM_X*NUM_Y;
        res[i] = res_arr[myYindex*NUM_X+myXindex];
    }

    { // de-allocate the arrays
        delete[] a;
        delete[] b;
        delete[] c;

        delete[] y;
        delete[] u;
        delete[] v;

        delete[] scan_tmp;
    }

    return elapsed;
}


int main() {
    const real_t s0 = 0.03, strike = 0.03, t = 5.0, alpha = 0.2, nu = 0.6, beta = 0.5;
    real_t *strikes, *res;

    if( IS_GPU ) { 
        fprintf(stdout, "\n// GPU All-Level Massively-Parallel ");
        fprintf(stdout, "Execution of Volatility Calibration Benchmark:\n"); 
    } else { 
        fprintf(stdout, "\n// CPU All-Level Multi-Threaded     ");
        fprintf(stdout, "Execution of Volatility Calibration Benchmark:\n"); 
    }


    readDataSet( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T );
    NUM_XY = NUM_X*NUM_Y;

    strikes = new real_t[OUTER_LOOP_COUNT];
    res     = new real_t[OUTER_LOOP_COUNT];
    allocGlobArrs();

    for(unsigned i=0; i<OUTER_LOOP_COUNT; ++i) {
        strikes[i] = 0.001*i;
    }

    unsigned long int elapsed_usec = 0;
    { // Instrumenting Runtime and Validation!
      // since OpenCL compilation or device querring hangs
      // for a long time many time we measure the runtime of
      // the time-series only in both CPU and GPU cases!

        elapsed_usec = whole_loop_nest( res, strikes, s0, t, alpha, nu, beta );
    }

    {
      FILE* runtime = fopen(getenv("HIPERMARK_RUNTIME"), "w");
      FILE* result = fopen(getenv("HIPERMARK_RESULT"), "w");
      fprintf(runtime, "%d\n", elapsed_usec);
      fclose(runtime);
      write_1Darr(result, res, OUTER_LOOP_COUNT);
      fclose(result);
    }

    delete[] strikes;
    delete[] res;
    deallocGlobArrs();

    return 0;
}
