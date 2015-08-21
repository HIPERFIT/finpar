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
const real_t g_mins [GENOME_DIM] = { EPS0,      EPS0,      -1.0+EPS0, EPS0, EPS0  }; 
const real_t g_maxs [GENOME_DIM] = { 1.0-EPS0,  1.0-EPS0,   1.0-EPS0, 0.2,  0.2   };
const real_t g_inis [GENOME_DIM] = { 0.02,      0.02,       0.0,      0.01, 0.04  };
//const real_t Candidate::p_ref[GENOME_DIM] = { 1.0, -2.0, 0.5, -0.3, -0.5, 0.1 };


/**
 * perturbing a genome; requires one random uniform number in [0,1)
 */
inline
real_t perturbation(  
                const real_t  gene,
                const real_t  gene_k,
                const real_t  gene_l,
                const UINT  i, 
                const real_t  gamma1,
                const real_t  amplitude_ratio
) {
    real_t amplitude     = fabs( (g_maxs[i] - g_mins[i]) * amplitude_ratio );
    real_t semiamplitude = amplitude / 2.0;
    real_t r01           = getRandRandNorm();
    real_t perturb       = ( amplitude * r01 - semiamplitude );

    return ( gene + perturb + gamma1 * ( gene_k - gene_l ) );
}

/*
 * contraining gene i to be within accepted bounds!
 */
inline 
real_t constrain_dim1( const int i, const real_t gene ) {
    return std::max( g_mins[i], std::min( g_maxs[i], gene ) );
}


/**
 * Helper mutate function: requires one uniform random number
 */
Tuple<real_t,real_t> mutate_helper( const real_t gene, const real_t gene_prop, const int i, const real_t amplitude_ratio ) {
    real_t forward_range, backward_range;
    real_t tmp_min_max,   tmp_max_min;

    const real_t amplitude     = fabs( (g_maxs[i] - g_mins[i]) * amplitude_ratio );
    const real_t semiamplitude = amplitude / 2.0;

    tmp_min_max = std::min( g_maxs[i], gene + semiamplitude );
    tmp_max_min = std::max( g_mins[i], gene - semiamplitude );
    forward_range = tmp_min_max - tmp_max_min;

    tmp_min_max = std::min( g_maxs[i], gene_prop + semiamplitude );
    tmp_max_min = std::max( g_mins[i], gene_prop - semiamplitude );
    backward_range = tmp_min_max - tmp_max_min;

    const real_t bf_fact = ( semiamplitude > 0.0 ) ? (backward_range / forward_range) : 1.0;

    // assign p'
    real_t r01 = getRandRandNorm();
    real_t diff = amplitude * r01 - semiamplitude;
    return Tuple<real_t,real_t>(gene + diff, bf_fact);
}


struct Genome {
    real_t a; 
    real_t b;
    real_t rho; 
    real_t nu;
    real_t sigma;
    real_t logLik;
    real_t fbRat;
    real_t padding;

    Genome() {
        this->a  = 0.0; this->b     = 0.0; this->rho    = 0.0; 
        this->nu = 0.0; this->sigma = 0.0; this->logLik = 0.0;
        this->fbRat = 0.0; 
    }
    Genome( const real_t a, const real_t b, const real_t rho, const real_t nu, 
            const real_t sigma, const real_t logLik, const real_t fbRat) {
        this->a  = a;  this->b     = b;     this->rho    = rho; 
        this->nu = nu; this->sigma = sigma; this->logLik = logLik;
        this->fbRat = fbRat; 
    } 
    Genome(const Genome& gene) {
        this->a  = gene.a;  this->b     = gene.b;     this->rho    = gene.rho; 
        this->nu = gene.nu; this->sigma = gene.sigma; this->logLik = gene.logLik; 
        this->fbRat = gene.fbRat;
    }
    Genome& operator=(const Genome& gene) {
        this->a  = gene.a;  this->b     = gene.b;     this->rho    = gene.rho; 
        this->nu = gene.nu; this->sigma = gene.sigma; this->logLik = gene.logLik;
        this->fbRat = gene.fbRat; 
        return (*this);
    }


    /**
     * mutate all dimensions. martingale with uniform prior.
     *   requires 5 * POP_SIZE uniform random numbers!
     * I took out the evaluation of the proposal, since it should be
     *   done for all types of crossover/mutation/etc
     */
    void mutate_dims_all( const Genome&  orig,
                          const real_t     amplitude_ratio = MOVES_UNIF_AMPL_RATIO ) {
        real_t fb_rat = 1.0;
        Tuple<real_t, real_t> tmp(EPS0, 1.0);

        tmp         = mutate_helper( orig.a,     this->a,     0, amplitude_ratio );
        this->a     = constrain_dim1(0, tmp.fst);           
        fb_rat     *= tmp.snd;
    
        tmp         = mutate_helper( orig.b,     this->b,     1, amplitude_ratio );
        this->b     = constrain_dim1(1, tmp.fst);
        fb_rat     *= tmp.snd;

        tmp         = mutate_helper( orig.rho,   this->rho,   2, amplitude_ratio );
        this->rho   = constrain_dim1(2, tmp.fst); 
        fb_rat     *= tmp.snd;

        tmp         = mutate_helper( orig.nu,    this->nu,    3, amplitude_ratio );
        this->nu    = constrain_dim1(3, tmp.fst); 
        fb_rat     *= tmp.snd;

        tmp         = mutate_helper( orig.sigma, this->sigma, 4, amplitude_ratio );
        this->sigma = constrain_dim1(4, tmp.fst);  
        fb_rat     *= tmp.snd;

        this->fbRat = fb_rat;
    }

    /**
     * mutate ONE dimension (Gibbs sampling). martingale with uniform prior
     *   dimension to be changed is given as input: has to be the same for 
     *   all population. Memory coalescence?
     * Requires POP_SIZE random uniform numbers!
     */
    void mutate_dims_one(   const Genome& orig,
                            const int     dim_j,
                            const real_t    amplitude_ratio = MOVES_UNIF_AMPL_RATIO 
    ) {    
        real_t fb_rat = 1.0;

        Tuple<real_t, real_t> tmp(EPS0, 1.0);
        tmp   = mutate_helper( orig.a, this->a, 0, (dim_j == 0) ? amplitude_ratio : 0.0 );
        this->a = constrain_dim1(0, tmp.fst); 
        fb_rat *= tmp.snd;
        
        tmp = mutate_helper( orig.b, this->b, 1, (dim_j == 1) ? amplitude_ratio : 0.0 );
        this->b = constrain_dim1(1, tmp.fst); 
        fb_rat *= tmp.snd;
    
        tmp = mutate_helper( orig.rho, this->rho, 2, (dim_j == 2) ? amplitude_ratio : 0.0 );
        this->rho = constrain_dim1(2, tmp.fst); 
        fb_rat *= tmp.snd;
    
        tmp = mutate_helper( orig.nu, this->nu, 3, (dim_j == 3) ? amplitude_ratio : 0.0 );
        this->nu = constrain_dim1(3, tmp.fst); 
        fb_rat *= tmp.snd;
    
        tmp         = mutate_helper( orig.sigma, this->sigma, 4, (dim_j == 4) ? amplitude_ratio : 0.0 );
        this->sigma = constrain_dim1(4, tmp.fst); 
        fb_rat     *= tmp.snd;

        this->fbRat = fb_rat;
    }

    /**
     * Crossover
     */
    void mcmc_DE(   const Genome& orig,
                    const Genome& gene_k,
                    const Genome& gene_l,
                    const real_t  gamma_avg = 2.38 / sqrt(2.0*GENOME_DIM),
                    const real_t  ampl_ratio = 0.1 * MOVES_UNIF_AMPL_RATIO
    ) {
        // proposal
        //   gamma: integrated out from the adviced 
        //          Multivariate Gaussian with Gaussian target (Braak, 2006)
        real_t gamma1 = gamma_avg - 0.5 + getRandUnifNorm();

        this->a     = constrain_dim1( 0, perturbation( orig.a,     gene_k.a,     gene_l.a,     0, gamma1, ampl_ratio ) );
        this->b     = constrain_dim1( 1, perturbation( orig.b,     gene_k.b,     gene_l.b,     1, gamma1, ampl_ratio ) );
        this->rho   = constrain_dim1( 2, perturbation( orig.rho,   gene_k.rho,   gene_l.rho,   2, gamma1, ampl_ratio ) );
        this->nu    = constrain_dim1( 3, perturbation( orig.nu,    gene_k.nu,    gene_l.nu,    3, gamma1, ampl_ratio ) );
        this->sigma = constrain_dim1( 4, perturbation( orig.sigma, gene_k.sigma, gene_l.sigma, 4, gamma1, ampl_ratio ) );
        
        this->fbRat = 1.0; // TODO: fix
    }

    /**
     * accepting the proposal
     */
    void accept( const Genome& mutated ) {
        this->a      = mutated.a;
        this->b      = mutated.b;
        this->rho    = mutated.rho;
        this->nu     = mutated.nu;
        this->sigma  = mutated.sigma;
        this->logLik = mutated.logLik;
        this->fbRat  = mutated.fbRat;
    }
};

#endif // #ifndef CANDIDATE_CLASS




