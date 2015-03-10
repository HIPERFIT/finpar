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

    unsigned long int elapsed = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        calib_arr = mainKernelSeqCPU(gene);

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {   // validation and writeback of the result
        bool is_valid = validate( gene.logLik, calib_arr, NUM_SWAP_QUOTES );
        writeStatsAndResult(    is_valid, gene.a, gene.b, gene.sigma, gene.nu, gene.rho, 
                                gene.logLik, calib_arr, NUM_SWAP_QUOTES, 
                                false, 1, elapsed );        
    }

	return 1;
}

