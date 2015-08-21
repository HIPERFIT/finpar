#ifndef GEN_ALG_UTIL   
#define GEN_ALG_UTIL   

#include "Constants.h"

//=========================================================================

/*****************************************/
/************** RANDOM NUMBERS ***********/
/*****************************************/

real_t getSobolNum( uint n ) {
    int n_gray = (n >> 1) ^ n;
    int res = 0;
    for( int i=0; i < NUM_SOBOL_BITS; i++ ) {
        int  t    = (1 << i);
        bool cond = ( (n_gray & t) == t ); 
        if ( cond ) {
            res = res ^ SobolDirVct[i];
        }
    }
    real_t rres = static_cast<real_t>(res) / ( (1<<NUM_SOBOL_BITS) + 1.0);
    if ( rres >= 1.0 || rres < 0 ) { printf("sobol(%d) = %f. Exiting!\n", n, rres); exit(0);}
    return rres;
}

//unsigned int rand_count = 1;
uint sobol_offset = 1; //11;

real_t getRandRandNorm() {
#if WITH_SOBOL
    return getSobolNum( sobol_offset++ );
#endif
    double r   = static_cast<double>(std::rand()); // drand48()
    double d   = static_cast<double>(RAND_MAX)+0.1;
    return static_cast<real_t>(r / d);
}

real_t getRandUnifNorm() {
    return getRandRandNorm();
    //return static_cast<real_t>(drand48());
}

//Returns a (pseudo) random integer in [0, ub)
UINT getRandIntNorm(long int ub) {  
    real_t r01 = getRandRandNorm();
    return static_cast<UINT>(r01 * ub);
}



/********************************************/
/************ Other Things!     *************/
/********************************************/

//=========================================================================
const real_t sqrtTwoPi = sqrt(2*PI);

real_t normal_pdf( const real_t& z, const real_t& mu, const real_t& sigma ) {
    real_t sigmap = fabs(sigma);
    real_t res    = 1.0 / (sigmap * sqrtTwoPi);
    real_t ecf    = (z-mu) * (z-mu) / (2.0 * sigmap * sigmap);
    return res * exp( -ecf );
}


real_t cauchy_pdf( const real_t& z, const real_t mu = 0.0, const real_t gamma = 4.0) {
    real_t x = (z-mu) / gamma;
    return 1.0 / ( PI * gamma * (1+x*x) );
}

real_t logLikelihood_normal(const real_t& y_ref, const real_t& y) {
    //real_t mu    = y_ref;
    real_t sigma = (y_ref / 50.0) * LLHOOD_NORMAL_OFFS;
    real_t pdfs  = normal_pdf( y, y_ref, sigma );
    pdfs += 1.0e-20; // avoid NaNs

    return log(pdfs);
}

real_t logLikelihood_cauchy(const real_t& y_ref, const real_t& y) {
    real_t gamma = ( fabs(y_ref) / 50.0 ) * LLHOOD_CAUCHY_OFFS + 0.01;
    real_t pdfs  = cauchy_pdf( y, y_ref, gamma );
    pdfs += 1.0e-20; // avoid NaNs
    return log(pdfs);
}

real_t logLikelihood(const real_t& y_ref, const real_t& y) {
    if   ( LLHOOD_TYPE == NORMAL ) return logLikelihood_normal(y_ref, y);
    else/* LLHOOD_TYPE == CAUCHY */return logLikelihood_cauchy(y_ref, y);
}
//=========================================================================

#endif  // #ifndef GEN_ALG_UTIL

