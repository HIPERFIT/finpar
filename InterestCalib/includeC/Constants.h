#ifndef CONSTS
#define CONSTS

#include <stdio.h>
#include <math.h>
#include "KerConsts.h"

#ifdef __APPLE__
typedef unsigned uint;
#endif

typedef uint UINT;

/*************************************/
/******** DATE RELATED STUFF *********/
/*************************************/

#if WITH_FLOAT
    const REAL INFTY = 1.0e38;
#else
    const REAL INFTY = 1.0e49;
#endif

bool equalEps(REAL x1, REAL x2) {
    return ( fabs(x1-x2) <= 1.0e-8 );
}

struct Date {
    int year;
    int month;
    int day;
    int hour;
    int min;

    Date(int y, int mo, int d, int h, int mn) : 
        year(y), month(mo), day(d), hour(h), min(mn) { } 
};

//"2299-12-31T23:59:59"
REAL MAX_DATE = 168307199.0;

//"1980-01-01T12:00:00"
REAL MIN_DATE = 3600.0;


#include "Date.h"

// Date.of_string("2012-01-01")
const REAL TODAY = static_cast<REAL>( date_of_gregorian( Date(2012, 1, 1, 12, 0) ) );


/*************************************/
/*******      GRID PARAMS      *******/
/*************************************/
const uint WPTH     = 8;
const uint lgLWG_FB = 8;
const uint   LWG_FB = 256;


/*************************************/
/******** MATH RELATED STUFF *********/
/*************************************/
//const REAL R     = 0.03;
const uint IT_MAX= 10000;  // a.k.a. itMax

template<class T>
struct Triple {
    T fst; T snd; T thrd;
    Triple(T a, T b, T c) : fst(a), snd(b), thrd(c) { }
};

struct SwapOfSwap {
    const UINT  n; // size of swap_sched
    REAL* swap_sched1;
    REAL* swap_sched2;
    REAL  maturity;
    REAL  strike;

    SwapOfSwap(const UINT len) : n(len) {
        swap_sched1 = new REAL[n];
        swap_sched2 = new REAL[n];
    } 
    void cleanUp() { 
        delete[] swap_sched1; 
        delete[] swap_sched2; 
    }
};

//-------------------------------------------------------------------------
// Gaussian Quadrature with Hermite linear expansion: cmulative distribution function of Normal distribution

UINT  NUM_HERMITE;    // 11
REAL* HermiteCoeffs; // [NUM_HERMITE]
REAL* HermiteWeights;// [NUM_HERMITE]

//=========================================================================
UINT NUM_SWAP_QUOTES; //196
// 196 x { maturity_in_year * swap_frequency * swap_term_in_year * volatility }
REAL* SwaptionQuotes; //[NUM_SWAP_QUOTES,4]

UINT NUM_SOBOL_BITS; // 30
int* SobolDirVct;  // [SOBOL_BITS_NUM]

int PROCS;

/******************************************/
/******* Genetic Algorithm Stuff **********/
/******************************************/

const UINT SEED         = 12345;
const UINT GENOME_DIM   = 5;
const REAL GENOME_SCALE = 1.0;


// Likelihood Parameters
enum Likelihood_Type { CAUCHY, NORMAL };
const Likelihood_Type LLHOOD_TYPE = CAUCHY; //NORMAL; 
const REAL LLHOOD_CAUCHY_OFFS     = 5.0;
const REAL LLHOOD_NORMAL_OFFS     = 1.0;


// population
//const UINT POP_SIZE = 128;
UINT POP_SIZE;

// MCMC
//const UINT MCMC_LOOPS = (UINT) ( 100.0 * (1.0 + log(GENOME_SCALE)) * 
//                                         (2.0 + log(GENOME_DIM  )) + 64);//was 2.0+...
UINT MCMC_LOOPS;
const REAL MOVES_UNIF_AMPL_RATIO = 0.005 / GENOME_SCALE;

// Generic Tuple
template<class T1, class T2>
struct Tuple {
    T1 fst;
    T2 snd;
    Tuple(const T1& t1, const T2& t2) : fst(t1), snd(t2) {}
};

enum Move_Type { DIMS_ALL, DIMS_ONE, DEMCMC, NONE };
typedef Tuple<REAL,Move_Type> MV_EL_TYPE;
const UINT CUMDENSFCT_CARD = 3;
MV_EL_TYPE mcmc_moves_selection_cumdensfct[3] = 
//{ MV_EL_TYPE(0.0,DIMS_ALL), MV_EL_TYPE(0.0,DIMS_ONE), MV_EL_TYPE(1.0,DEMCMC) };
    { MV_EL_TYPE(0.2,DIMS_ALL), MV_EL_TYPE(0.5,DIMS_ONE), MV_EL_TYPE(1.0,DEMCMC) };
// was: 0.2                       0.5                       1.0


const int LWG = 128;

#endif // ifndef CONSTS

