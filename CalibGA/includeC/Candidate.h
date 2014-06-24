#ifndef CANDIDATE_CLASS
#define CANDIDATE_CLASS

#include "Constants.h"
#include "GenAlgUtil.h"

//=========================================================================
/**
 * Candidate:
 *   - identified by parameters (target, and nuisance)
 *   - holds one proposal, to be accepted or rejected
 *   - holds fitness results
 *   - holds MCMC selection regions
 * Population:
 *   - container of candidates
 *   - container of proposals
 *   - caller of mutation functions
 *   - select if to keep/update candidates and proposals
 * MCMC: moves:
 *   - TODO: Gibbs sampling: single dimension update
 *   - all dimensions update
 *   - crossover
 */
//=========================================================================

/**
 * Genome min/max ranges and proposed initial values (not used)
 * { g_a, g_b, g_rho, g_nu, g_sigma }
 */
const REAL g_mins [GENOME_DIM] = { EPS0,      EPS0,      -1.0+EPS0, EPS0, EPS0  }; 
const REAL g_maxs [GENOME_DIM] = { 1.0-EPS0,  1.0-EPS0,   1.0-EPS0, 0.2,  0.2   };
const REAL g_inis [GENOME_DIM] = { 0.02,      0.02,       0.0,      0.01, 0.04  };
//const REAL Candidate::p_ref[GENOME_DIM] = { 1.0, -2.0, 0.5, -0.3, -0.5, 0.1 };


/**
 * perturbing a genome; requires one random uniform number in [0,1)
 */
inline
REAL perturbation(  
                const REAL  gene,
                const REAL  gene_k,
                const REAL  gene_l,
                const UINT  i, 
                const REAL  gamma1,
                const REAL  amplitude_ratio
) {
    REAL amplitude     = fabs( (g_maxs[i] - g_mins[i]) * amplitude_ratio );
    REAL semiamplitude = amplitude / 2.0;
    REAL r01           = getRandRandNorm();
    REAL perturb       = ( amplitude * r01 - semiamplitude );

    return ( gene + perturb + gamma1 * ( gene_k - gene_l ) );
}

/*
 * contraining gene i to be within accepted bounds!
 */
inline 
REAL constrain_dim1( const int i, const REAL gene ) {
    return std::max( g_mins[i], std::min( g_maxs[i], gene ) );
}

/**
 * accepting the proposal
 */
void accept(const int i,  
            REAL* g_a,     // genome begins
            REAL* g_b, 
            REAL* g_rho, 
            REAL* g_nu, 
            REAL* g_sigma,
            REAL* logLik
) {
    g_a     [i] = g_a    [i+POP_SIZE];
    g_b     [i] = g_b    [i+POP_SIZE];
    g_rho   [i] = g_rho  [i+POP_SIZE];
    g_nu    [i] = g_nu   [i+POP_SIZE];
    g_sigma [i] = g_sigma[i+POP_SIZE];
    logLik  [i] = logLik [i+POP_SIZE];
}

/**
 * Helper mutate function: requires one uniform random number
 */
Tuple<REAL,REAL> mutate_helper( const REAL gene, const REAL gene_prop, const int i, const REAL amplitude_ratio ) {
    REAL forward_range, backward_range;
    REAL tmp_min_max,   tmp_max_min;

    const REAL amplitude     = fabs( (g_maxs[i] - g_mins[i]) * amplitude_ratio );
    const REAL semiamplitude = amplitude / 2.0;

    tmp_min_max = std::min( g_maxs[i], gene + semiamplitude );
    tmp_max_min = std::max( g_mins[i], gene - semiamplitude );
    forward_range = tmp_min_max - tmp_max_min;

    tmp_min_max = std::min( g_maxs[i], gene_prop + semiamplitude );
    tmp_max_min = std::max( g_mins[i], gene_prop - semiamplitude );
    backward_range = tmp_min_max - tmp_max_min;

    const REAL bf_fact = ( semiamplitude > 0.0 ) ? (backward_range / forward_range) : 1.0;

    // assign p'
    REAL r01 = getRandRandNorm();
    REAL diff = amplitude * r01 - semiamplitude;
    return Tuple<REAL,REAL>(gene + diff, bf_fact);
}

/**
 * mutate all dimensions. martingale with uniform prior.
 *   requires 5 * POP_SIZE uniform random numbers!
 * I took out the evaluation of the proposal, since it should be
 *   done for all types of crossover/mutation/etc
 */
void mutate_dims_all(   REAL* g_a,     // genome begins
                        REAL* g_b, 
                        REAL* g_rho, 
                        REAL* g_nu, 
                        REAL* g_sigma,
                        REAL* g_fb_rat,
                  const int   i,
                  const REAL  amplitude_ratio = MOVES_UNIF_AMPL_RATIO ) {
    REAL fb_rat = 1.0;

    Tuple<REAL, REAL> tmp(EPS0, 1.0);
    tmp = mutate_helper( g_a    [i], g_a    [i+POP_SIZE], 0, amplitude_ratio );
    g_a     [i + POP_SIZE] = constrain_dim1(0, tmp.fst); fb_rat *= tmp.snd;
    
    tmp = mutate_helper( g_b    [i], g_b    [i+POP_SIZE], 1, amplitude_ratio );
    g_b     [i + POP_SIZE] = constrain_dim1(1, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_rho  [i], g_rho  [i+POP_SIZE], 2, amplitude_ratio );
    g_rho   [i + POP_SIZE] = constrain_dim1(2, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_nu   [i], g_nu   [i+POP_SIZE], 3, amplitude_ratio );
    g_nu    [i + POP_SIZE] = constrain_dim1(3, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_sigma[i], g_sigma[i+POP_SIZE], 4, amplitude_ratio );
    g_sigma [i + POP_SIZE] = constrain_dim1(4, tmp.fst); fb_rat *= tmp.snd;

    g_fb_rat[i] = fb_rat;
}


/**
 * mutate ONE dimension (Gibbs sampling). martingale with uniform prior
 *   dimension to be changed is given as input: has to be the same for 
 *   all population. Memory coalescence?
 * Requires POP_SIZE random uniform numbers!
 */
void mutate_dims_one(   REAL* g_a,     // genome begins
                        REAL* g_b, 
                        REAL* g_rho, 
                        REAL* g_nu, 
                        REAL* g_sigma,
                        REAL* g_fb_rat,
                  const int   i,
                  const int   dim_j,
                  const REAL  amplitude_ratio = MOVES_UNIF_AMPL_RATIO 
) {    
    REAL fb_rat = 1.0;

    Tuple<REAL, REAL> tmp(EPS0, 1.0);
    tmp = mutate_helper( g_a    [i], g_a    [i+POP_SIZE], 0, (dim_j == 0) ? amplitude_ratio : 0.0 );
    g_a     [i + POP_SIZE] = constrain_dim1(0, tmp.fst); fb_rat *= tmp.snd;
    
    tmp = mutate_helper( g_b    [i], g_b    [i+POP_SIZE], 1, (dim_j == 1) ? amplitude_ratio : 0.0 );
    g_b     [i + POP_SIZE] = constrain_dim1(1, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_rho  [i], g_rho  [i+POP_SIZE], 2, (dim_j == 2) ? amplitude_ratio : 0.0 );
    g_rho   [i + POP_SIZE] = constrain_dim1(2, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_nu   [i], g_nu   [i+POP_SIZE], 3, (dim_j == 3) ? amplitude_ratio : 0.0 );
    g_nu    [i + POP_SIZE] = constrain_dim1(3, tmp.fst); fb_rat *= tmp.snd;

    tmp = mutate_helper( g_sigma[i], g_sigma[i+POP_SIZE], 4, (dim_j == 4) ? amplitude_ratio : 0.0 );
    g_sigma [i + POP_SIZE] = constrain_dim1(4, tmp.fst); fb_rat *= tmp.snd;

    g_fb_rat[i] = fb_rat;
}

/**
 * Crossover
 * Requires a vector of size [ POP_SIZE*8 ] of random uniform numbers.
 */
void mcmc_DE(               REAL* g_a,     // genome begins
                            REAL* g_b, 
                            REAL* g_rho, 
                            REAL* g_nu, 
                            REAL* g_sigma, // genome ends
                            REAL* bf_rat,  // fwd_bwd_ratio
                      const int   j,
                      const REAL  gamma_avg = 2.38 / sqrt(2.0*GENOME_DIM),
                      const REAL  ampl_ratio = 0.1 * MOVES_UNIF_AMPL_RATIO
) {
        UINT cand_UB = POP_SIZE - 1;
        UINT k = getRandIntNorm(cand_UB); // random in [0,pop_size-1)
        if ( k == j ) {
            k = cand_UB;
            cand_UB -= 1;
        }
        UINT l = getRandIntNorm(cand_UB); 
        if ( l == j || l == k ) {
            l = cand_UB;
        }            

        // proposal
        //   gamma: integrated out from the adviced 
        //          Multivariate Gaussian with Gaussian target (Braak, 2006)
        REAL gamma1 = gamma_avg - 0.5 + getRandUnifNorm();

        g_a     [j + POP_SIZE] = constrain_dim1( 0, perturbation( g_a    [j], g_a    [k], g_a    [l], 0, gamma1, ampl_ratio ) );
        g_b     [j + POP_SIZE] = constrain_dim1( 1, perturbation( g_b    [j], g_b    [k], g_b    [l], 1, gamma1, ampl_ratio ) );
        g_rho   [j + POP_SIZE] = constrain_dim1( 2, perturbation( g_rho  [j], g_rho  [k], g_rho  [l], 2, gamma1, ampl_ratio ) );
        g_nu    [j + POP_SIZE] = constrain_dim1( 3, perturbation( g_nu   [j], g_nu   [k], g_nu   [l], 3, gamma1, ampl_ratio ) );
        g_sigma [j + POP_SIZE] = constrain_dim1( 4, perturbation( g_sigma[j], g_sigma[k], g_sigma[l], 4, gamma1, ampl_ratio ) );
        
        bf_rat[j] = 1.0; // TODO: fix
}

#endif // #ifndef CANDIDATE_CLASS




