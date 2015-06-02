#ifndef GEN_ALG_UTIL   
#define GEN_ALG_UTIL   

#include "Constants.h"

//=========================================================================

/*****************************************/
/************** RANDOM NUMBERS ***********/
/*****************************************/

REAL getSobolNum( uint n ) {
    int n_gray = (n >> 1) ^ n;
    int res = 0;
    for( int i=0; i < NUM_SOBOL_BITS; i++ ) {
        int  t    = (1 << i);
        bool cond = ( (n_gray & t) == t ); 
        if ( cond ) {
            res = res ^ SobolDirVct[i];
        }
    }
    REAL rres = static_cast<REAL>(res) / ( (1<<NUM_SOBOL_BITS) + 1.0);
    if ( rres >= 1.0 || rres < 0 ) { printf("sobol(%d) = %f. Exiting!\n", n, rres); exit(0);}
    return rres;
}

//unsigned int rand_count = 1;
uint sobol_offset = 1; //11;

REAL getRandRandNorm() {
#if WITH_SOBOL
    return getSobolNum( sobol_offset++ );
#endif
    double r   = static_cast<double>(std::rand()); // drand48()
    double d   = static_cast<double>(RAND_MAX)+0.1;
    return static_cast<REAL>(r / d);
}

REAL getRandUnifNorm() {
    return getRandRandNorm();
    //return static_cast<REAL>(drand48());
}

//Returns a (pseudo) random integer in [0, ub)
UINT getRandIntNorm(long int ub) {  
    REAL r01 = getRandRandNorm();
    return static_cast<UINT>(r01 * ub);
}



/********************************************/
/************ Other Things!     *************/
/********************************************/

//=========================================================================
const REAL sqrtTwoPi = sqrt(2*PI);

REAL normal_pdf( const REAL& z, const REAL& mu, const REAL& sigma ) {
    REAL sigmap = fabs(sigma);
    REAL res    = 1.0 / (sigmap * sqrtTwoPi);
    REAL ecf    = (z-mu) * (z-mu) / (2.0 * sigmap * sigmap);
    return res * exp( -ecf );
}


REAL cauchy_pdf( const REAL& z, const REAL mu = 0.0, const REAL gamma = 4.0) {
    REAL x = (z-mu) / gamma;
    return 1.0 / ( PI * gamma * (1+x*x) );
}

REAL logLikelihood_normal(const REAL& y_ref, const REAL& y) {
    //REAL mu    = y_ref;
    REAL sigma = (y_ref / 50.0) * LLHOOD_NORMAL_OFFS;
    REAL pdfs  = normal_pdf( y, y_ref, sigma );
    pdfs += 1.0e-20; // avoid NaNs

    return log(pdfs);
}

REAL logLikelihood_cauchy(const REAL& y_ref, const REAL& y) {
    REAL gamma = ( fabs(y_ref) / 50.0 ) * LLHOOD_CAUCHY_OFFS + 0.01;
    REAL pdfs  = cauchy_pdf( y, y_ref, gamma );
    pdfs += 1.0e-20; // avoid NaNs
    return log(pdfs);
}

REAL logLikelihood(const REAL& y_ref, const REAL& y) {
    if   ( LLHOOD_TYPE == NORMAL ) return logLikelihood_normal(y_ref, y);
    else/* LLHOOD_TYPE == CAUCHY */return logLikelihood_cauchy(y_ref, y);
}
//=========================================================================

#endif  // #ifndef GEN_ALG_UTIL

