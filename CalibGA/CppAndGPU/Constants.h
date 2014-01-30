#ifndef CONSTS
#define CONSTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <time.h>
#include <sys/timeb.h>
#include <omp.h>


#include "KerConsts.h"

#ifdef __APPLE__
typedef unsigned uint;
#endif

typedef uint UINT;

typedef struct timeb mlfi_timeb;
#define mlfi_ftime ftime

#define mlfi_diff_time(t1,t2) \
  (t1.time - t2.time) * 1000 + (t1.millitm - t2.millitm)


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

//assert max_date==168307199, a.k.a. "2299-12-31T23:59:59"
REAL MAX_DATE = 168307199.0;

//assert min_date==3600, a.k.a., "1980-01-01T12:00:00"
REAL MIN_DATE = 3600.0;


#include "Date.h"

// Date.of_string("2012-01-01")
const REAL TODAY = static_cast<REAL>( date_of_gregorian( Date(2012, 1, 1, 12, 0) ) );


/*************************************/
/*******      GRID PARAMS      *******/
/*************************************/
const uint WPTH     = 8;
const uint lgLWG_FB = 7;
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

const UINT Gauss_DIM = 11;
class GaussHermiteCoefs {
  public:
    static const REAL coefs  [Gauss_DIM];
    static const REAL weights[Gauss_DIM];
};

const REAL GaussHermiteCoefs::coefs[Gauss_DIM] = {  
            0.0, 0.6568095668820999044613, -0.6568095668820997934390, -1.3265570844949334805563,
                 1.3265570844949330364670,  2.0259480158257567872226, -2.0259480158257558990442,
                -2.7832900997816496513337,  2.7832900997816474308877,  3.6684708465595856630159, 
                -3.6684708465595838866591 
        };
const REAL GaussHermiteCoefs::weights[Gauss_DIM] = {  
            0.6547592869145917315876, 0.6609604194409607336169, 0.6609604194409606225946, 
            0.6812118810666693002887, 0.6812118810666689672217, 0.7219536247283847574252, 
            0.7219536247283852015144, 0.8025168688510405656800, 0.8025168688510396775015, 
            1.0065267861723647957461, 1.0065267861723774522886 
        };
//=========================================================================
const UINT NUM_SWAP_QUOTES = 196;
// 196 x { maturity_in_year * swap_frequency * swap_term_in_year * volatility }
const REAL SwaptionQuotes[NUM_SWAP_QUOTES][4] = {
        {1.0, 6.0, 1.0, 1.052},
        {2.0, 6.0, 1.0, 0.81485},
        {3.0, 6.0, 1.0, 0.6165},
        {4.0, 6.0, 1.0, 0.46995},
        {5.0, 6.0, 1.0, 0.38295},
        {6.0, 6.0, 1.0, 0.3325},
        {7.0, 6.0, 1.0, 0.3016},
        {8.0, 6.0, 1.0, 0.2815},
        {9.0, 6.0, 1.0, 0.26435},
        {10.0,6.0, 1.0, 0.2496},
        {15.0,6.0, 1.0, 0.2516},
        {20.0,6.0, 1.0, 0.28835},
        {25.0,6.0, 1.0, 0.27155},
        {30.0,6.0, 1.0, 0.23465},

        {1.0, 6.0, 2.0, 0.61445},
        {2.0, 6.0, 2.0, 0.54805},
        {3.0, 6.0, 2.0, 0.46795},
        {4.0, 6.0, 2.0, 0.3919},
        {5.0, 6.0, 2.0, 0.3434},
        {6.0, 6.0, 2.0, 0.3083},
        {7.0, 6.0, 2.0, 0.28655},
        {8.0, 6.0, 2.0, 0.2697},
        {9.0, 6.0, 2.0, 0.25775},
        {10.0,6.0, 2.0, 0.2443},
        {15.0,6.0, 2.0, 0.26495},
        {20.0,6.0, 2.0, 0.28195},
        {25.0,6.0, 2.0, 0.26845},
        {30.0,6.0, 2.0, 0.20995},

        {1.0, 6.0, 3.0, 0.5835},
        {2.0, 6.0, 3.0, 0.49255},
        {3.0, 6.0, 3.0, 0.42825},
        {4.0, 6.0, 3.0, 0.3695},
        {5.0, 6.0, 3.0, 0.329},
        {6.0, 6.0, 3.0, 0.3022},
        {7.0, 6.0, 3.0, 0.28165},
        {8.0, 6.0, 3.0, 0.26615},
        {9.0, 6.0, 3.0, 0.25485},
        {10.0,6.0, 3.0, 0.24375},
        {15.0,6.0, 3.0, 0.2718},
        {20.0,6.0, 3.0, 0.28135},
        {25.0,6.0, 3.0, 0.26865},
        {30.0,6.0, 3.0, 0.2131},

        {1.0, 6.0, 4.0, 0.5415},
        {2.0, 6.0, 4.0, 0.46235},
        {3.0, 6.0, 4.0, 0.403},
        {4.0, 6.0, 4.0, 0.3559},
        {5.0, 6.0, 4.0, 0.3232},
        {6.0, 6.0, 4.0, 0.29675},
        {7.0, 6.0, 4.0, 0.27715},
        {8.0, 6.0, 4.0, 0.26385},
        {9.0, 6.0, 4.0, 0.254},
        {10.0,6.0, 4.0, 0.2454},
        {15.0,6.0, 4.0, 0.27845},
        {20.0,6.0, 4.0, 0.2821},
        {25.0,6.0, 4.0, 0.2678},
        {30.0,6.0, 4.0, 0.2131},

        {1.0, 6.0, 5.0, 0.517},
        {2.0, 6.0, 5.0, 0.446},
        {3.0, 6.0, 5.0, 0.3903},
        {4.0, 6.0, 5.0, 0.34755},
        {5.0, 6.0, 5.0, 0.3166},
        {6.0, 6.0, 5.0, 0.29305},
        {7.0, 6.0, 5.0, 0.2745},
        {8.0, 6.0, 5.0, 0.2639},
        {9.0, 6.0, 5.0, 0.2534},
        {10.0, 6.0, 5.0, 0.2499},
        {15.0, 6.0, 5.0, 0.28315},
        {20.0, 6.0, 5.0, 0.2825},
        {25.0, 6.0, 5.0, 0.277},
        {30.0, 6.0, 5.0, 0.21175},

        {1.0, 6.0, 6.0, 0.478},
        {2.0, 6.0, 6.0, 0.42105},
        {3.0, 6.0, 6.0, 0.37715},
        {4.0, 6.0, 6.0, 0.3378},
        {5.0, 6.0, 6.0, 0.311},
        {6.0, 6.0, 6.0, 0.2895},
        {7.0, 6.0, 6.0, 0.2745},
        {8.0, 6.0, 6.0, 0.264},
        {9.0, 6.0, 6.0, 0.2573},
        {10.0, 6.0, 6.0, 0.25475},
        {15.0, 6.0, 6.0, 0.28815},
        {20.0, 6.0, 6.0, 0.28195},
        {25.0, 6.0, 6.0, 0.26015},
        {30.0, 6.0, 6.0, 0.2097},

        {1.0, 6.0, 7.0, 0.452},
        {2.0, 6.0, 7.0, 0.4074},
        {3.0, 6.0, 7.0, 0.368},
        {4.0, 6.0, 7.0, 0.3307},
        {5.0, 6.0, 7.0, 0.30645},
        {6.0, 6.0, 7.0, 0.2877},
        {7.0, 6.0, 7.0, 0.27475},
        {8.0, 6.0, 7.0, 0.2664},
        {9.0, 6.0, 7.0, 0.26155},
        {10.0, 6.0, 7.0, 0.26035},
        {15.0, 6.0, 7.0, 0.292},
        {20.0, 6.0, 7.0, 0.2825},
        {25.0, 6.0, 7.0, 0.25685},
        {30.0, 6.0, 7.0, 0.2081},

        {1.0, 6.0, 8.0, 0.43395},
        {2.0, 6.0, 8.0, 0.39445},
        {3.0, 6.0, 8.0, 0.35885},
        {4.0, 6.0, 8.0, 0.3281},
        {5.0, 6.0, 8.0, 0.30395},
        {6.0, 6.0, 8.0, 0.28745},
        {7.0, 6.0, 8.0, 0.2767},
        {8.0, 6.0, 8.0, 0.27065},
        {9.0, 6.0, 8.0, 0.26625},
        {10.0, 6.0, 8.0, 0.26625},
        {15.0, 6.0, 8.0, 0.2921},
        {20.0, 6.0, 8.0, 0.2814},
        {25.0, 6.0, 8.0, 0.25265},
        {30.0, 6.0, 8.0, 0.2083},

        {1.0, 6.0, 9.0, 0.42285},
        {2.0, 6.0, 9.0, 0.3857},
        {3.0, 6.0, 9.0, 0.3521},
        {4.0, 6.0, 9.0, 0.3239},
        {5.0, 6.0, 9.0, 0.30285},
        {6.0, 6.0, 9.0, 0.2895},
        {7.0, 6.0, 9.0, 0.2799},
        {8.0, 6.0, 9.0, 0.27485},
        {9.0, 6.0, 9.0, 0.2712},
        {10.0, 6.0, 9.0, 0.27205},
        {15.0, 6.0, 9.0, 0.29205},
        {20.0, 6.0, 9.0, 0.27855},
        {25.0, 6.0, 9.0, 0.24945},
        {30.0, 6.0, 9.0, 0.219},

        {1.0, 6.0, 10.0, 0.41765},
        {2.0, 6.0, 10.0, 0.38095},
        {3.0, 6.0, 10.0, 0.34795},
        {4.0, 6.0, 10.0, 0.3217},
        {5.0, 6.0, 10.0, 0.30365},
        {6.0, 6.0, 10.0, 0.2916},
        {7.0, 6.0, 10.0, 0.2842},
        {8.0, 6.0, 10.0, 0.27985},
        {9.0, 6.0, 10.0, 0.2769},
        {10.0, 6.0, 10.0, 0.2775},
        {15.0, 6.0, 10.0, 0.306},
        {20.0, 6.0, 10.0, 0.2763},
        {25.0, 6.0, 10.0, 0.2458},
        {30.0, 6.0, 10.0, 0.22},

        {1.0, 6.0, 15.0, 0.37905},
        {2.0, 6.0, 15.0, 0.35465},
        {3.0, 6.0, 15.0, 0.33505},
        {4.0, 6.0, 15.0, 0.31725},
        {5.0, 6.0, 15.0, 0.3008},
        {6.0, 6.0, 15.0, 0.29075},
        {7.0, 6.0, 15.0, 0.28365},
        {8.0, 6.0, 15.0, 0.2787},
        {9.0, 6.0, 15.0, 0.27385},
        {10.0, 6.0, 15.0, 0.2709},
        {15.0, 6.0, 15.0, 0.2689},
        {20.0, 6.0, 15.0, 0.24225},
        {25.0, 6.0, 15.0, 0.2096},
        {30.0, 6.0, 15.0, 0.18285},

        {1.0, 6.0, 20.0, 0.37975},
        {2.0, 6.0, 20.0, 0.3605},
        {3.0, 6.0, 20.0, 0.3407},
        {4.0, 6.0, 20.0, 0.321},
        {5.0, 6.0, 20.0, 0.3063},
        {6.0, 6.0, 20.0, 0.29315},
        {7.0, 6.0, 20.0, 0.28395},
        {8.0, 6.0, 20.0, 0.2777},
        {9.0, 6.0, 20.0, 0.27205},
        {10.0, 6.0, 20.0, 0.26675},
        {15.0, 6.0, 20.0, 0.24875},
        {20.0, 6.0, 20.0, 0.21735},
        {25.0, 6.0, 20.0, 0.1939},
        {30.0, 6.0, 20.0, 0.17205},

        {1.0, 6.0, 25.0, 0.38115},
        {2.0, 6.0, 25.0, 0.3627},
        {3.0, 6.0, 25.0, 0.34425},
        {4.0, 6.0, 25.0, 0.3222},
        {5.0, 6.0, 25.0, 0.3084},
        {6.0, 6.0, 25.0, 0.2941},
        {7.0, 6.0, 25.0, 0.28285},
        {8.0, 6.0, 25.0, 0.2751},
        {9.0, 6.0, 25.0, 0.2663},
        {10.0, 6.0, 25.0, 0.26055},
        {15.0, 6.0, 25.0, 0.2338},
        {20.0, 6.0, 25.0, 0.20735},
        {25.0, 6.0, 25.0, 0.1823},
        {30.0, 6.0, 25.0, 0.1686},

        {1.0, 6.0, 30.0, 0.38285},
        {2.0, 6.0, 30.0, 0.3633},
        {3.0, 6.0, 30.0, 0.34125},
        {4.0, 6.0, 30.0, 0.3188},
        {5.0, 6.0, 30.0, 0.30305},
        {6.0, 6.0, 30.0, 0.2888},
        {7.0, 6.0, 30.0, 0.2748},
        {8.0, 6.0, 30.0, 0.26725},
        {9.0, 6.0, 30.0, 0.25985},
        {10.0, 6.0, 30.0, 0.25165},
        {15.0, 6.0, 30.0, 0.2267},
        {20.0, 6.0, 30.0, 0.1989},
        {25.0, 6.0, 30.0, 0.18115},
        {30.0, 6.0, 30.0, 0.16355}
};

const int SOBOL_BITS_NUM = 30;
int  SOBOL_DIR_VCT[30] =  {   536870912,268435456,134217728,67108864,33554432,16777216,8388608,4194304,
                        2097152,1048576,524288,262144,131072,65536,32768,16384,8192,
                        4096,2048,1024,512,256,128,64,32,16,8,4,2,1
                    }; 

int PROCS;

/******************************************/
/******* Genetic Algorithm Stuff **********/
/******************************************/

const UINT SEED         = 12345;
const UINT GENOME_DIM   = 5;
const REAL GENOME_SCALE = 1.0;


// Likelihood Parameters
enum Likelihood_Type { CAUCHY, NORMAL };
const Likelihood_Type LLHOOD_TYPE = CAUCHY; //NORMAL; //CAUCHY;
const REAL LLHOOD_CAUCHY_OFFS     = 5.0;
const REAL LLHOOD_NORMAL_OFFS     = 1.0;


// population
const UINT POP_SIZE = 128; //128; //128; //1024;//2048; //100;

// MCMC
const UINT MCMC_LOOPS = (UINT) ( 100.0 * (1.0 + log(GENOME_SCALE)) * (2.0 + log(GENOME_DIM)) );   //128;  // 160; //
const REAL MOVES_UNIF_AMPL_RATIO = 0.005 / GENOME_SCALE;

// Generic Tuple
template<class T1, class T2>
struct Tuple {
    T1 fst;
    T2 snd;

    Tuple(const T1& t1, const T2& t2) : fst(t1), snd(t2) {}
};

// 
enum Move_Type { DIMS_ALL, DIMS_ONE, DEMCMC, NONE };
typedef Tuple<REAL,Move_Type> MV_EL_TYPE;
const UINT CUMDENSFCT_CARD = 3;
MV_EL_TYPE mcmc_moves_selection_cumdensfct[3] = 
//                        { MV_EL_TYPE(0.0,DIMS_ALL), MV_EL_TYPE(0.0,DIMS_ONE), MV_EL_TYPE(1.0,DEMCMC) };
                        { MV_EL_TYPE(0.2,DIMS_ALL), MV_EL_TYPE(0.5,DIMS_ONE), MV_EL_TYPE(1.0,DEMCMC) };
// was:                                  0.2                       0.5                       1.0


const int LWG = 128;

#endif // ifndef CONSTS

