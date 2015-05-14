#ifndef GEN_ALG_FLAT
#define GEN_ALG_FLAT

using namespace std;

#include "Constants.h"
#include "GenAlgUtil.h"
#include "Genome.h"
#include "UtilCPU.h"
#include "EvalGenomeOrig.h"


/**
 * Printing Swaption / Calibrated Price / Black Price / RMS
 * The result is an array REAL[NUM_SWAP_QUOTES, 3] recording
 *   for each swaption the calibrated price, the black price 
 *   and the percentagewise difference between the two.
 */
void makeSummary(   const int     winner, 
                    const Genome& wgene, 
                    const REAL*   new_quote, 
                    const REAL*   new_price,
                          REAL*   result // output
) {
    REAL rms = 0.0;

    fprintf(stderr, "\n\nCALIBRATION RESULT: best genome is at index %d: ", winner);
    fprintf(stderr, "{ a = %f, b = %f, sigma = %f, nu = %f, rho = %f  }, Likelihood: %f!\n",
                    wgene.a, wgene.b, wgene.sigma, wgene.nu, wgene.rho, wgene.logLik );
    fprintf(stderr, "\nPer-Swaption Approximation w.r.t. Black Price:\n\n");

    for( int i = 0; i < NUM_SWAP_QUOTES; i ++ ) {
        REAL black_price = new_quote[i];
        REAL calib_price = new_price[i];
        REAL err_ratio   =  (calib_price - black_price) / black_price;

        result[3*i + 0] = 10000.0*calib_price;
        result[3*i + 1] = 10000.0*black_price;
        result[3*i + 2] = 100.0*fabs(err_ratio);

        rms += err_ratio * err_ratio;

        fprintf(stderr,"Swaption %d: {{%f, %f, %f},%f}, CalibratedPrice: %f, BlackPrice: %f, DiffPerc: %f\n",
                i, SwaptionQuotes[4*i+0], SwaptionQuotes[4*i+1], SwaptionQuotes[4*i+2], SwaptionQuotes[4*i+3], 
                result[3*i + 0], result[3*i + 1], result[3*i + 2] );
    }

    rms = 100.0 * sqrt ( rms / NUM_SWAP_QUOTES );
    fprintf(stderr, "\n\n Best Genome RMS: %f\n\n", rms);
}


/**
 * Utility function: find the genome with the best likelihood:
 *    scans the logLik array and fill in the index and likelihood
 *    of the best genome.
 */
void find_best(const Genome* genomes, int& best_ind, REAL& best_lik) {
    bool sanity = true;

    best_lik = -INFINITY;
    best_ind = 0;

    // this is in fact a reduction, but POP_SIZE is
    //     not big enough to warrant a parallel execution.
    for ( UINT i = 0; i < POP_SIZE; i++ ) {  // parallel reduction with MAX
        REAL val = genomes[i].logLik;

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
REAL* mainKernelSeqCPU( Genome& winner ) {

    UINT n_schedi_max = 0;
    for( UINT ttt = 0; ttt < NUM_SWAP_QUOTES; ttt++ ) {
        const UINT n_schedi = static_cast<UINT>(12.0 * SwaptionQuotes[4*ttt+2] / SwaptionQuotes[4*ttt+1]);
        n_schedi_max = max(n_schedi_max, n_schedi);
    }

    SeqArrays cpu_arrs(n_schedi_max);
    Genome* genomes = cpu_arrs.genomes;

    srand   ( SEED );
    srand48 ( SEED );

    // initialized the genomes with random numbers inside 
    //   their acceptable bounds => requires 5*POP_SIZE randoms
    for( int i = 0; i < POP_SIZE; i++ ) {
        REAL r01, tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[0] - g_mins[0]) + g_mins[0];
        genomes[i].a = genomes[i+POP_SIZE].a = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[1] - g_mins[1]) + g_mins[1];
        genomes[i].b = genomes[i+POP_SIZE].b = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[2] - g_mins[2]) + g_mins[2];
        genomes[i].rho = genomes[i+POP_SIZE].rho = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[3] - g_mins[3]) + g_mins[3];
        genomes[i].nu = genomes[i+POP_SIZE].nu = tmp;

        r01 = getRandRandNorm();
        tmp = r01 * ( g_maxs[4] - g_mins[4]) + g_mins[4];
        genomes[i].sigma = genomes[i+POP_SIZE].sigma = tmp;
        
//        genomes[i].fbRat = genomes[i+POP_SIZE].fbRat = 1.0;
    }

    // Initial evaluation of the genomes!
    REAL quote, price;         
    for( int i = 0; i < POP_SIZE; i++ ) {
        Genome& gene = genomes[i];
        REAL rms = 0.0;
        for( UINT ttt = 0; ttt < NUM_SWAP_QUOTES; ttt++ ) {
            eval_genome_new ( 
                        gene.a, gene.b, gene.rho, gene.nu, gene.sigma,
                        SwaptionQuotes+4*ttt, cpu_arrs.tmp_arrs, quote, price
                    );
            rms += logLikelihood( quote, price );
        }
        gene.logLik = rms;
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
                genomes[i+POP_SIZE].mutate_dims_all( genomes[i] );
            }

        } else if ( move_type == DIMS_ONE ) {

            UINT dim_j = getRandIntNorm(GENOME_DIM); 

            for( int i = 0; i < POP_SIZE; i++ ) {
                genomes[i+POP_SIZE].mutate_dims_one( genomes[i], dim_j );
            }
        } else /* if ( move_type == DEMCMC   ) */ {

            for ( int i = 0; i < POP_SIZE; i++ ) { // parallel modulo random nums
                // compute the k^th and the l^th genes
                UINT cand_UB = POP_SIZE - 1;
                UINT k = getRandIntNorm(cand_UB); // random in [0,pop_size-1)
                if ( k == i ) {
                    k = cand_UB;
                    cand_UB -= 1;
                }
                UINT l = getRandIntNorm(cand_UB); 
                if ( l == i || l == k ) {
                    l = cand_UB;
                }
 
                // do DEMCMC base on k^th and l^th genomes 
                genomes[i+POP_SIZE].mcmc_DE( genomes[i], genomes[k], genomes[l] );
            }
        }

        // evaluate the proposals
        for( int i = 0; i < POP_SIZE; i++ ) {
            Genome& gene = genomes[i+POP_SIZE];

            REAL rms = 0.0;
            for( UINT ttt = 0; ttt < NUM_SWAP_QUOTES; ttt++ ) {
                eval_genome_new ( 
                            gene.a, gene.b, gene.rho, gene.nu, gene.sigma,
                            SwaptionQuotes+4*ttt, cpu_arrs.tmp_arrs, quote, price
                        );
                rms += logLikelihood( quote, price );
            }
            gene.logLik = rms;
        }

        // mcmc_acceptance_rejection();
        // Deciding whether to accept or reject the proposal,
        // obtained by mutating/crossover of the individual.
        for ( int i = 0; i < POP_SIZE; i++ ) {  // parallel
            // Metropolis: get a random U[0,1) for each candidate
            REAL rand = getRandUnifNorm();
            
            Genome& orig = genomes[i];
            Genome& muta = genomes[i+POP_SIZE];
            // selection: dimensions considered independent
            // acceptance = min( 1, N.exp(c.logLik_proposal-c.logLik) * c.backward_forward_ratio )
            REAL acceptance = std::min( 1.0, exp( muta.logLik - orig.logLik ) * muta.fbRat );
          
            // if acceptance criterion is met then p->p' else does nothing 
            if ( rand < acceptance ) 
                orig = muta;
        }

        // print best candidate for the current iteration:
        if ( (j % 16) == 0 ){ 
            int best_ind; REAL best_lik;
            find_best(genomes, best_ind, best_lik);
            fprintf(stderr, "\n Iteration: %d: Best Likelihood: %f, genome index: %d!\n", 
                            j, best_lik, best_ind );
        }
    }

    REAL* result = new REAL[3*NUM_SWAP_QUOTES];
    { // print best candidate for the current iteration:
        int best_ind; REAL best_lik;
        find_best(genomes, best_ind, best_lik);
        //winner = genomes[best_ind];
        winner.a   = genomes[best_ind].a;   winner.b  = genomes[best_ind].b; 
        winner.rho = genomes[best_ind].rho; winner.nu = genomes[best_ind].nu;
        winner.sigma = genomes[best_ind].sigma;  winner.logLik = genomes[best_ind].logLik; 

        // recompute the calibrated price and the black price!
        REAL* quotes = cpu_arrs.get_quote();
        REAL* prices = cpu_arrs.get_price();
        for( UINT ttt = 0; ttt < NUM_SWAP_QUOTES; ttt++ ) {
            eval_genome_new ( 
                        winner.a, winner.b, winner.rho, winner.nu, winner.sigma,
                        SwaptionQuotes+4*ttt, cpu_arrs.tmp_arrs, quote, price
                    );
            quotes[ttt] = quote; 
            prices[ttt] = price;
        }

        

        // write summary
        makeSummary( best_ind, winner, quotes, prices, result );
    }

    // Releasing the CPU resources:
    cpu_arrs.releaseResources();

    return result;
}

#endif // end ifndef GEN_ALG_FLAT

