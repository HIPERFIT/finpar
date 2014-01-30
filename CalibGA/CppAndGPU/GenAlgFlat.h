#ifndef GEN_ALG_FLAT
#define GEN_ALG_FLAT

#include "Constants.h"
#include "GenAlgUtil.h"
#include "Candidate.h"

#include "IrregShape.h"
#include "UtilGPGPU.h"
#include "EvalGenomeInl.h"

#include "FindBestKer.h"
#include "MainKer.h"
#include "GenAlgKers.h"

/**
 * Printing Swaption / Calibrated Price / Black Price / RMS
 */
void printSummary(uint winner, CpuArrays& arrs) {
    REAL rms = 0.0;

    for( int i = 0; i < NUM_SWAP_QUOTES; i ++ ) {
        REAL black_price = arrs.get_quote(winner)[i];
        REAL calib_price = arrs.get_price(winner)[i];
        REAL err_ratio   =  (calib_price - black_price) / black_price;

        rms += err_ratio * err_ratio;

        printf( "Swaption %d: { { %f, %f, %f}, %f }, Calibrated Price: %f, Black Price: %f, PERCENT: %f\n",
                i, SwaptionQuotes[i][0], SwaptionQuotes[i][1], SwaptionQuotes[i][2], SwaptionQuotes[i][3], 
                10000.0*calib_price, 10000*black_price, 100.0*fabs(err_ratio) );
    }

    rms = 100.0 * sqrt ( rms / NUM_SWAP_QUOTES );

    printf( "Best Genome (index %d) is: { a = %f, b = %f, sigma = %f, nu = %f, rho = %f  }, Likelihood: %f, RMS: %f!\n",
            winner, arrs.get_a()[winner], arrs.get_b()[winner], arrs.get_sigma()[winner], 
            arrs.get_nu()[winner], arrs.get_rho()[winner], arrs.get_logLik()[winner], rms );
}


/**
 * Utility function: find the genome with the best likelihood:
 *    scans the logLik array and fill in the index and likelihood
 *    of the best genome.
 */
void find_best(const REAL* logLik, int& best_ind, REAL& best_lik) {
    bool sanity = true;

    best_lik = -INFINITY;
    best_ind = 0;

    // this is in fact a reduction, but POP_SIZE is
    //     not big enough to warrant a parallel execution.
    for ( UINT i = 0; i < POP_SIZE; i++ ) {  // parallel reduction with MAX
        REAL val = logLik[i];

        sanity = !( isnan(val) || isinf(val) );
        assert( sanity && "val is NaN in find_best" );

        if( val > best_lik ) { best_ind = i; best_lik = val; }
    }
}

void find_best_GPU(
            FindBestKer& bestker, OclObjects& ocl_objs, 
            OclBuffers& ocl_arrs, const REAL* logLik, 
            int& best_ind,        REAL& best_lik
) {
    bool sanity = true;

    bestker.run(best_ind, best_lik);

    sanity = ! ( isnan(best_lik) || isinf(best_lik) );
    assert( sanity && "val is NaN in find_best" );
}

void print_best(    REAL* g_a,     // genome begins
                    REAL* g_b, 
                    REAL* g_rho, 
                    REAL* g_nu, 
                    REAL* g_sigma,
                    REAL* logLik,
              const int   best_ind,
              const REAL  best_lik,
              const int   j
) {
    // print best candidate for the current iteration:
    printf("Finished iteration: %d, best cand (ind, likelihood): (%d, %.9f)\n", j, best_ind, best_lik);
    printf("Best genome: (a,b,rho,nu,sigma) : (%f,%f,%f,%f,%f)\n", 
            g_a[best_ind], g_b[best_ind], g_rho[best_ind], g_nu[best_ind], g_sigma[best_ind] );
}

Move_Type selectMoveType(REAL move_selected) {
        Move_Type move_type     = NONE;

        REAL      prob;
        Move_Type type;
        UINT      k = 0;
        do {
            prob = mcmc_moves_selection_cumdensfct[k].fst;
            type = mcmc_moves_selection_cumdensfct[k].snd;

            if( move_selected <= prob ) {
                move_type = type;
            }

            k ++;
        } while ( move_selected > prob && k < CUMDENSFCT_CARD );

        assert(move_type != NONE && "MOVE_TYPE == NONE is ILLEGAL!");

        return move_type;
}

void mainKernel (  ) {
    uint   FLAT_SZ;
    short* shape      = getIregShapeAdjusted( LWG_EG, FLAT_SZ );
    int  * start_inds = getStartInd(FLAT_SZ, shape);

    CpuArrays cpu_arrs(FLAT_SZ, shape);
    REAL *g_a, *g_b, *g_rho, *g_nu, *g_sigma, *logLik, *bf_rat;
    { // getting the cpu arrays
        g_a     = cpu_arrs.get_a     ();
        g_b     = cpu_arrs.get_b     ();
        g_rho   = cpu_arrs.get_rho   ();
        g_nu    = cpu_arrs.get_nu    ();
        g_sigma = cpu_arrs.get_sigma ();
        logLik  = cpu_arrs.get_logLik();
        bf_rat  = cpu_arrs.get_bf_rat();
    }

    OclObjects ocl_objs; OclBuffers ocl_arrs; 
    { // Initialization of the GPU resources:
        initGPUresources ( cpu_arrs, ocl_objs, ocl_arrs );
    }

    FindBestKer bestker  (POP_SIZE, 10*POP_SIZE, ocl_arrs.genomes, ocl_arrs.best_ind, ocl_arrs.best_val, ocl_objs, 128);
    MainKer     mainker  ( cpu_arrs, ocl_arrs, ocl_objs );
    GenAlgKers  genalgker( cpu_arrs, ocl_arrs, ocl_objs );

    srand   ( SEED );
    srand48 ( SEED );

    mlfi_timeb  t_start, t_end;
    unsigned long int elapsed;
    mlfi_ftime(&t_start);

    // initialized the genomes with random numbers inside 
    //   their acceptable bounds => requires 5*POP_SIZE randoms
#if (GPU_VERSION == 2)
    genalgker.run_init( );
#else

    for( int i = 0; i < POP_SIZE; i++ ) {
        REAL r01, tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[0] - g_mins[0]) + g_mins[0];
        g_a     [i           ] = tmp;
        g_a     [i + POP_SIZE] = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[1] - g_mins[1]) + g_mins[1];
        g_b     [i           ] = tmp;
        g_b     [i + POP_SIZE] = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[2] - g_mins[2]) + g_mins[2];
        g_rho   [i           ] = tmp;
        g_rho   [i + POP_SIZE] = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[3] - g_mins[3]) + g_mins[3];
        g_nu    [i           ] = tmp;
        g_nu    [i + POP_SIZE] = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[4] - g_mins[4]) + g_mins[4];
        g_sigma [i           ] = tmp;
        g_sigma [i + POP_SIZE] = tmp;
        
        bf_rat  [i           ] = 1.0;
    }
#endif

    // Initial evaluation of the genomes!
#if (GPU_VERSION == 2)
          mainker.run      (); 
        genalgker.run_cpLik();

#elif (GPU_VERSION == 1)
        mainker.run(); //eval_genome_GPU( mainker, cpu_arrs );        
        for( int i = 0; i < POP_SIZE; i++ ) {
            logLik[i] = logLik[i+POP_SIZE];
        }
#else
        {
            int i;
#pragma omp parallel for default(shared) schedule(dynamic) private(i) if(POP_SIZE>3)
            for( i = 0; i < POP_SIZE; i++ ) {
                logLik[i] = eval_genome( g_a[i], g_b[i], g_rho[i], g_nu[i], g_sigma[i], 
                                         cpu_arrs.get_quote(i), cpu_arrs.get_price(i), 
                                         FLAT_SZ, shape, start_inds );
            }
        }
#endif

    // sequential loop that runs the genetic algorithms
    // (takes a population and returns a population)
    for( int j = 0; j < MCMC_LOOPS; j++ ) {

        // select which move to perform: BEGIN
        // Note: this block can also be a loop (in fixed order) 
        //       over the various move types
        REAL      move_selected = getRandUnifNorm();
        Move_Type move_type     = selectMoveType(move_selected);

#if ( GPU_VERSION == 2 && WITH_SOBOL == 0 )
        sobol_offset ++;
#endif

        if ( move_type == DIMS_ALL ) {
#if (GPU_VERSION == 2)
            genalgker.run_mutate( 33 );
#else
            for( int i = 0; i < POP_SIZE; i++ ) {
                mutate_dims_all( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i );  // default ampl_ratio
            }
#endif

       } else if ( move_type == DIMS_ONE ) {

            UINT dim_j = getRandIntNorm(GENOME_DIM); 

#if (GPU_VERSION == 2)
    #if (WITH_SOBOL == 0)
            sobol_offset++;
    #endif
            genalgker.run_mutate( dim_j );
#else
            for( int i = 0; i < POP_SIZE; i++ ) {
                mutate_dims_one( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i, dim_j ); // default ampl_ratio
            }
#endif
        } else /* if ( move_type == DEMCMC   ) */ {
#if (GPU_VERSION == 2)
            genalgker.run_McMcDc( );
#else
            for ( int i = 0; i < POP_SIZE; i++ ) { // parallel modulo random nums
                mcmc_DE( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i );
            }
#endif
        }

        // evaluate the proposal
#if (GPU_VERSION == 2 || GPU_VERSION == 1)
        mainker.run(); 
#else
        {
            int i;
#pragma omp parallel for default(shared) schedule(dynamic) private(i) if(POP_SIZE>3)
            for( i = 0; i < POP_SIZE; i++ ) {
                logLik[i+POP_SIZE] = eval_genome( g_a [i+POP_SIZE], g_b    [i+POP_SIZE], g_rho[i+POP_SIZE], 
                                                  g_nu[i+POP_SIZE], g_sigma[i+POP_SIZE], 
                                                  cpu_arrs.get_quote(i), cpu_arrs.get_price(i), 
                                                  FLAT_SZ, shape, start_inds );
            }
        }
#endif

        // mcmc_acceptance_rejection();
        // Deciding whether to accept or reject the proposal,
        // obtained by mutating/crossover of the individual.
#if (GPU_VERSION == 2)
        genalgker.run_accept_cond( );
        genalgker.run_accept_prop( );
#else
        for ( int i = 0; i < POP_SIZE; i++ ) {  // parallel
            // Metropolis: get a random U[0,1) for each candidate
            // rand=N.random.uniform()
            REAL rand = getRandUnifNorm();
            
            // selection: dimensions considered independent
            // acceptance = min( 1, N.exp(c.logLik_proposal-c.logLik) * c.backward_forward_ratio )
            REAL acceptance = std::min( 1.0, exp( logLik[i+POP_SIZE] - logLik[i] ) * bf_rat[i] );
          
            // if acceptance criterion is met then p->p' else does nothing 
            if ( rand < acceptance ) accept( i, g_a, g_b, g_rho, g_nu, g_sigma, logLik );
        }
#endif

        // print best candidate for the current iteration:
#if (GPU_VERSION == 2)
        if ( (j % 16) == 0 ) { 
            int best_ind; REAL best_lik;
            find_best_GPU( bestker, ocl_objs, ocl_arrs, logLik, best_ind, best_lik );
            printf( "\n Iteration: %d: Best Likelihood: %f, genome index: %d!\n", 
                    j, best_lik, best_ind );
        }
#else
        if ( (j % 16) == 0 ){ 
            int best_ind; REAL best_lik;
            find_best(logLik, best_ind, best_lik);
            printf( "\n Iteration: %d: Best Likelihood: %f, genome index: %d!\n", 
                    j, best_lik, best_ind );
        }
#endif
    }

    mlfi_ftime(&t_end);
    elapsed = mlfi_diff_time(t_end,t_start);

#if (GPU_VERSION == 2)
    { // copy back genomes
        uint cur_size = 12 * POP_SIZE * sizeof(REAL);
        uint ciErr;
        cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();
        ciErr = clEnqueueReadBuffer(cmd_queue, ocl_arrs.genomes, CL_TRUE, 0, cur_size, cpu_arrs.genomes, 0, NULL, NULL);
        oclCheckError( ciErr, CL_SUCCESS );
    }    
#endif

    { // print best candidate for the current iteration:
        int best_ind; REAL best_lik;
        printf("\n\n");
        find_best(logLik, best_ind, best_lik);
//        print_best( g_a, g_b, g_rho, g_nu, g_sigma, logLik, best_ind, best_lik, MCMC_LOOPS-1 );

        // recompute the calibrated price and the black price!
        logLik[best_ind] = eval_genome( g_a[best_ind], g_b[best_ind], g_rho[best_ind], g_nu[best_ind], g_sigma[best_ind], 
                                         cpu_arrs.get_quote(best_ind), cpu_arrs.get_price(best_ind), FLAT_SZ, shape, start_inds );

        // Finally, print result!
        printSummary(best_ind, cpu_arrs);
    }

#if   (GPU_VERSION == 2)
    printf("\nALL-ON-GPU Calibration Time: %lu !\n\n", elapsed);
#elif (GPU_VERSION == 1)
    printf("\nMOSTLY-GPU Calibration Time: %lu !\n\n", elapsed);
#else
    printf("\nPAR-ON-CPU Calibration Time: %lu !\n\n", elapsed);
#endif


    delete[] start_inds;

    // Releasing the GPU resources:
    releaseGPUresources ( ocl_objs, ocl_arrs );
    
    // Releasing the CPU resources:
    cpu_arrs.releaseResources();
}

#endif // end ifndef

