#ifndef GEN_ALG_FLAT
#define GEN_ALG_FLAT

#include "Constants.h"
#include "GenAlgUtil.h"
#include "Candidate.h"

#include "IrregShape.h"
#include "EvalGenomeInl.h"
#include "UtilCPU.h"

/**
 * Printing Swaption / Calibrated Price / Black Price / RMS
 * The result is an array REAL[NUM_SWAP_QUOTES, 3] recording
 *   for each swaption the calibrated price, the black price 
 *   and the percentagewise difference between the two.
 */
REAL* makeSummary(uint winner, CpuArrays& arrs) {
    REAL* res = (REAL*) malloc( 3*NUM_SWAP_QUOTES*sizeof(REAL) );

    REAL rms = 0.0;

    fprintf(stderr, "\n\nCALIBRATION RESULT: best genome is at index %d: ", winner);
    fprintf(stderr, "{ a = %f, b = %f, sigma = %f, nu = %f, rho = %f  }, Likelihood: %f!\n",
                    arrs.get_a() [winner], arrs.get_b()  [winner], arrs.get_sigma ()[winner], 
                    arrs.get_nu()[winner], arrs.get_rho()[winner], arrs.get_logLik()[winner] );
    fprintf(stderr, "\nPer-Swaption Approximation w.r.t. Black Price:\n\n");

    for( int i = 0; i < NUM_SWAP_QUOTES; i ++ ) {
        REAL black_price = arrs.get_quote(winner)[i];
        REAL calib_price = arrs.get_price(winner)[i];
        REAL err_ratio   =  (calib_price - black_price) / black_price;

        res[3*i + 0] = 10000.0*calib_price;
        res[3*i + 1] = 10000.0*black_price;
        res[3*i + 2] = 100.0*fabs(err_ratio);

        rms += err_ratio * err_ratio;

        fprintf(stderr,"Swaption %d: {{%f, %f, %f},%f}, CalibratedPrice: %f, BlackPrice: %f, DiffPerc: %f\n",
                i, SwaptionQuotes[4*i+0], SwaptionQuotes[4*i+1], SwaptionQuotes[4*i+2], SwaptionQuotes[4*i+3], 
                res[3*i + 0], res[3*i + 1], res[3*i + 2] );
    }

    rms = 100.0 * sqrt ( rms / NUM_SWAP_QUOTES );
    fprintf(stderr, "\n\n Best Genome RMS: %f\n\n", rms);

    return res;
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

/**
 * Main Entry Point for Swaption Calibration
 * The (out) arguments are in fact the resulted 
 *      winning genome result & its likelihood.
 * The array result has size 3*NUM_SWAP_QUOTES
 *      and records for each swaption the calibrated 
 *      price, the black price and the percentage 
 *      difference between the two.
 */
REAL* mainKernelCPU(REAL& wg_a, 
                    REAL& wg_b, 
                    REAL& wg_sigma, 
                    REAL& wg_nu, 
                    REAL& wg_rho, 
                    REAL& wg_logLik
) {
    uint   FLAT_SZ;
    short* shape      = getIregShapeAdjusted( LWG_EG, FLAT_SZ );
    int  * start_inds = getStartInd( FLAT_SZ, shape, NUM_SWAP_QUOTES );

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

    srand   ( SEED );
    srand48 ( SEED );

    // initialized the genomes with random numbers inside 
    //   their acceptable bounds => requires 5*POP_SIZE randoms
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

    // Initial evaluation of the genomes!           
#pragma omp parallel for default(shared) schedule(dynamic) if(POP_SIZE>3)
    for( int i = 0; i < POP_SIZE; i++ ) {
        logLik[i] = eval_genome( g_a[i], g_b[i], g_rho[i], g_nu[i], g_sigma[i], 
                                 cpu_arrs.get_quote(i), cpu_arrs.get_price(i), 
                                 FLAT_SZ, shape, start_inds 
                               );
    }

    // convergence loop that runs the genetic algorithms
    // (takes a population and returns a population)
    for( int j = 0; j < MCMC_LOOPS; j++ ) {

        // select which move to perform.
        // Note: this block can also be a loop (in fixed order) 
        //       over the various move types
        REAL      move_selected = getRandUnifNorm();
        Move_Type move_type     = selectMoveType(move_selected);

        if ( move_type == DIMS_ALL ) {
            for( int i = 0; i < POP_SIZE; i++ ) {
                mutate_dims_all( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i );
            }

        } else if ( move_type == DIMS_ONE ) {

            UINT dim_j = getRandIntNorm(GENOME_DIM); 

            for( int i = 0; i < POP_SIZE; i++ ) {
                mutate_dims_one( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i, dim_j );
            }
        } else /* if ( move_type == DEMCMC   ) */ {

            for ( int i = 0; i < POP_SIZE; i++ ) { // parallel modulo random nums
                mcmc_DE( g_a, g_b, g_rho, g_nu, g_sigma, bf_rat, i );
            }
        }

        // evaluate the proposal
#pragma omp parallel for default(shared) schedule(dynamic) if(POP_SIZE>3)
        for( int i = 0; i < POP_SIZE; i++ ) {
            logLik[i+POP_SIZE] = 
                eval_genome( g_a [i+POP_SIZE], g_b    [i+POP_SIZE], g_rho[i+POP_SIZE], 
                             g_nu[i+POP_SIZE], g_sigma[i+POP_SIZE], 
                             cpu_arrs.get_quote(i), cpu_arrs.get_price(i), 
                             FLAT_SZ, shape, start_inds );
        }

        // mcmc_acceptance_rejection();
        // Deciding whether to accept or reject the proposal,
        // obtained by mutating/crossover of the individual.
        for ( int i = 0; i < POP_SIZE; i++ ) {  // parallel
            // Metropolis: get a random U[0,1) for each candidate
            REAL rand = getRandUnifNorm();
            
            // selection: dimensions considered independent
            // acceptance = min( 1, N.exp(c.logLik_proposal-c.logLik) * c.backward_forward_ratio )
            REAL acceptance = std::min( 1.0, exp( logLik[i+POP_SIZE] - logLik[i] ) * bf_rat[i] );
          
            // if acceptance criterion is met then p->p' else does nothing 
            if ( rand < acceptance ) accept( i, g_a, g_b, g_rho, g_nu, g_sigma, logLik );
        }

        // print best candidate for the current iteration:
        if ( (j % 16) == 0 ){ 
            int best_ind; REAL best_lik;
            find_best(logLik, best_ind, best_lik);
            fprintf(stderr, "\n Iteration: %d: Best Likelihood: %f, genome index: %d!\n", 
                            j, best_lik, best_ind );
        }
    }

    REAL* result;
    { // print best candidate for the current iteration:
        int best_ind; REAL best_lik;
        find_best(logLik, best_ind, best_lik);

        // recompute the calibrated price and the black price!
        logLik[best_ind] = 
            eval_genome( g_a[best_ind], g_b[best_ind], g_rho[best_ind], g_nu[best_ind], 
                         g_sigma[best_ind], cpu_arrs.get_quote(best_ind), 
                         cpu_arrs.get_price(best_ind), FLAT_SZ, shape, start_inds );

        // Finally, make/print the result!
        wg_a      = cpu_arrs.get_a()     [best_ind];
        wg_b      = cpu_arrs.get_b()     [best_ind]; 
        wg_sigma  = cpu_arrs.get_sigma ()[best_ind]; 
        wg_nu     = cpu_arrs.get_nu()    [best_ind];
        wg_rho    = cpu_arrs.get_rho()   [best_ind];
        wg_logLik = cpu_arrs.get_logLik()[best_ind];

        result = makeSummary( best_ind, cpu_arrs );
    }

    delete[] start_inds;    
    // Releasing the CPU resources:
    cpu_arrs.releaseResources();

    return result;
}

#endif // end ifndef

