/****************************************/
/*** Translated to C++ Cosmin Oancea  ***/
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
    Genome gene;
    REAL*  calib_arr = NULL;

    printf("\n// Running Original Sequential Swaption-Calibration Benchmark\n");

    readDataSet(    POP_SIZE, MCMC_LOOPS,
                    NUM_SWAP_QUOTES, SwaptionQuotes, 
                    NUM_HERMITE, HermiteCoeffs, HermiteWeights,
                    NUM_SOBOL_BITS, SobolDirVct
               ); 

    unsigned long int elapsed_usec = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        calib_arr = mainKernelSeqCPU(gene);

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed_usec = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {   // validation and writeback of the result
        FILE* runtime = fopen("runtime.txt", "w");
        FILE* result = fopen("result.json", "w");
        const int Ps = get_CPU_num_threads();
        fprintf(runtime, "%d\n", elapsed_usec / 1000);
        fclose(runtime);

        writeResult(result,
                    gene.a,
                    gene.b,
                    gene.sigma,
                    gene.nu,
                    gene.rho,
                    gene.logLik,
                    calib_arr);

        fclose(result);
    }

	return 1;
}
