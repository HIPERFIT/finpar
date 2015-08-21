/****************************************/
/*** Translated to C++ & parallelized ***/
/*** in OpenCL by Cosmin E. Oancea    ***/
/***   from an original OCaml code    ***/
/***    (sequential) by LexiFi        ***/
/***   and a sequential Python implem ***/
/***    using a different genetic alg.***/
/***    by Christian Andreetta        ***/  
/****************************************/ 


#include "Util.h"
#include "Constants.h"
#include "ParseInput.h"
#include "GenAlgFlat.h"
#include "WriteResult.h"

int main()
{
    real_t wg_a  = 0.0, wg_b   = 0.0, wg_sigma  = 0.0, 
           wg_nu = 0.0, wg_rho = 0.0, wg_logLik = 0.0;
    real_t* calib_arr = NULL;

    printf("\n// Running GPU-Parallel Swaption-Calibration Benchmark\n");

    readDataSet(    POP_SIZE, MCMC_LOOPS,
                    NUM_SWAP_QUOTES, SwaptionQuotes, 
                    NUM_HERMITE, HermiteCoeffs, HermiteWeights,
                    NUM_SOBOL_BITS, SobolDirVct
               ); 

    unsigned long int elapsed_usec = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        calib_arr = mainKernelGPU(wg_a, wg_b, wg_sigma, wg_nu, wg_rho, wg_logLik);

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed_usec = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {   // writeback of the result
        FILE* runtime = fopen(getenv("HIPERMARK_RUNTIME"), "w");
        FILE* result = fopen(getenv("HIPERMARK_RESULT"), "w");
        const int Ps = get_CPU_num_threads();
        fprintf(runtime, "%d\n", elapsed_usec / 1000);
        fclose(runtime);

        writeResult(result,
                    wg_a,
                    wg_b,
                    wg_sigma,
                    wg_nu,
                    wg_rho,
                    wg_logLik,
                    calib_arr);

        fclose(result);
    }

    return 0;
}
