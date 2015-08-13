#ifndef PARSE_INPUT
#define PARSE_INPUT

#include "ParserC.h"
#include <assert.h>
#include <math.h>
#include <algorithm>
#include "GenAlgUtil.h"

//using namespace std;

#if WITH_FLOAT
    #define read_real read_float
#else
    #define read_real read_double
#endif

/***********************************/
/********** READ DATA SET **********/
/***********************************/

void readDataSet(   UINT&  pop_size,
                    UINT&  num_conv_iters,   
                    UINT&  num_swap_quotes, 
                    REAL*& swaption_quotes,  // [num_swap_quotes,4]
                    
                    UINT&  num_hermitians,
                    REAL*& hermite_coefs,   // [num_hermite]
                    REAL*& hermite_weights,  // [num_hermite]
                
                    UINT&  num_sobol_bits,
                    int*&  sobol_dirs_vct   // [num_sobol_bits]
) {
    int64_t shape[3];
    bool atr_ok = true;

    if ( read_int( static_cast<UINT*>( &pop_size       ) ) ||
         read_int( static_cast<UINT*>( &num_conv_iters ) )  ) {
            fprintf(stderr, "Syntax error when reading population size or convergence-loop count!\n");
            exit(1);
    }
    atr_ok  = (pop_size >= 64) && (num_conv_iters >= 100);
    assert(atr_ok && "Population size < 64 OR convergence-loop count < 100!\n");

    { // swaption quotes
        if ( read_int( static_cast<UINT*>( &num_swap_quotes ) ) ) {
            fprintf(stderr, "Syntax error when reading NUM_SWAP_QUOTES.\n");
            exit(1);
        }
        atr_ok  = num_swap_quotes > 0;
        assert(atr_ok && "Number of swaptions LESS or equal to ZERO!\n");

        if ( read_array(sizeof(REAL), read_real, (void**)&swaption_quotes, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading the swaption quotes.\n");
            exit(1);
        }
        atr_ok = ( shape[1] == 4 && shape[0] == num_swap_quotes );
        assert(atr_ok && "Incorrect shape of the swaption quotes!");
    }

    { // hermitean coefficients
        if ( read_int( static_cast<UINT*>( &num_hermitians ) ) ) {
            fprintf(stderr, "Syntax error when reading number of hermitians.\n");
            exit(1);
        }
        atr_ok  = num_hermitians > 0;
        assert(atr_ok && "Number of hermitians LESS or equal to ZERO!\n");

        if ( read_array(sizeof(REAL), read_real, (void**)&hermite_coefs, shape, 1) ) {
            fprintf(stderr, "Syntax error when reading the hermitian coeficients.\n");
            exit(1);
        }
        atr_ok = ( shape[0] == num_hermitians );
        assert(atr_ok && "Incorrect shape of the hermitian coefficients!");

        if ( read_array(sizeof(REAL), read_real, (void**)&hermite_weights, shape, 1) ) {
            fprintf(stderr, "Syntax error when reading the hermitian weights.\n");
            exit(1);
        }
        atr_ok = ( shape[0] == num_hermitians );
        assert(atr_ok && "Incorrect shape of the hermitian weights!");
    }

    { // sobol direction vectors!
        if ( read_int( static_cast<UINT*>( &num_sobol_bits ) ) ) {
            fprintf(stderr, "Syntax error when reading NUM_SOBOL_BITS.\n");
            exit(1);
        }
        atr_ok  = num_sobol_bits > 0;
        assert(atr_ok && "Number of sobol bits LESS or equal to ZERO!\n");

        if ( read_array(sizeof(int), read_int, (void**)&sobol_dirs_vct, shape, 1) ) {
            fprintf(stderr, "Syntax error when reading sobol direction vectors.\n");
            exit(1);
        }
        atr_ok = ( shape[0] == num_sobol_bits );
        assert( atr_ok && "Incorrect shape of sobol direction vectors!" );
    }
}

REAL* readOutput( const int& N, REAL& ref_logLik ) {
    REAL*   result;
    int64_t shape[3];
    
    REAL ref_a, ref_b, ref_sigma, ref_nu, ref_rho;

    // reading the reference genome!
    if(     read_real( &ref_a )     || read_real( &ref_b  )     ||
            read_real( &ref_sigma ) || read_real( &ref_nu )     ||
            read_real( &ref_rho )   || read_real( &ref_logLik )  ) {
        fprintf(stderr, "Syntax error when reading the reference genome or its likelihood.\n");
        exit(1);
    }

    // reading the standard output
    if ( read_array(sizeof(REAL), read_real, (void**)&result, shape, 2) ) {
        fprintf(stderr, "Syntax error when reading the output.\n");
        exit(1);
    }
    bool ok = ( shape[0] == N ) && ( shape[1] == 3 );
    assert(ok && "Incorrect shape of the reference array result!");

    return result;
}

bool validate( const REAL& logLik, const REAL* res, const int& N ) {
    bool  is_valid = true;
    
    REAL  ref_logLik;
    REAL* ref_res = readOutput( N, ref_logLik );

    if( fabs(ref_logLik - logLik) > 1.0 ) {
        is_valid = false;
        fprintf(stderr, "Difference in logLikelihood > 1.0: result logLik %f, reference logLik = %f!\n", 
                        logLik, ref_logLik );
        return is_valid;
    } else if ( logLik < ref_logLik && fabs(ref_logLik - logLik) > 0.5 ) {
        is_valid = false;
        fprintf(stderr, "Result likelihood is worse than reference with more than 0.5: result logLik %f, reference logLik = %f!\n", 
                        logLik, ref_logLik );
        return is_valid;
    }

    for ( int i = 0; i < N; i ++ ) {
        REAL lgLik_cur = logLikelihood(    res[3*i+1],    res[3*i]);
        REAL lgLik_ref = logLikelihood(ref_res[3*i+1],ref_res[3*i]);
        REAL diff_curr_lgLik = (lgLik_cur > lgLik_ref) ? 
                                0.0 : lgLik_ref - lgLik_cur ;

        if (diff_curr_lgLik > 0.7) {
            is_valid = false;
            fprintf(stderr, "Error[%d] = %f, Acceptable = 1.1!\n", 
                            i, diff_curr_lgLik );
            break;
        }
    }

    return is_valid;
}

void writeStatsAndResult(   const bool& valid,
                            const REAL& wg_a, 
                            const REAL& wg_b, 
                            const REAL& wg_sigma, 
                            const REAL& wg_nu, 
                            const REAL& wg_rho, 
                            const REAL& wg_logLik,
                            const REAL* calib_arr,
                            const int & NUM_SWAP_QUOTES,
                            const bool& is_gpu,
                            const int & P, 
                            const unsigned long int& elapsed  
) {
    // print stats to stdout
    fprintf(stdout, "// Dataset with Number of Swaption Quotes = %d.\n", 
                    NUM_SWAP_QUOTES );

    if(valid) { fprintf(stdout, "1\t\t// VALID   Result,\n"); } 
    else      { fprintf(stdout, "0\t\t// INVALID Result,\n"); }

    fprintf(stdout, "%ld\t\t// Runtime in microseconds,\n", elapsed);
    if(is_gpu) fprintf(stdout, "%d\t\t// GPU Threads,\n\n", P);
    else       fprintf(stdout, "%d\t\t// CPU Threads,\n\n", P);

    // write the genome
    write_scal(&wg_a,     "a-field      of the winning genome");
    write_scal(&wg_b,     "b-field      of the winning genome");
    write_scal(&wg_sigma, "sigma-field  of the winning genome");
    write_scal(&wg_nu,    "nu-field     of the winning genome");
    write_scal(&wg_rho,   "rho-field    of the winning genome");
    write_scal(&wg_logLik,"LgLikelihood of the winning genome");

    // write the result
    write_2Darr( calib_arr, static_cast<int>(NUM_SWAP_QUOTES), 3, 
                 "Swaption Calibration Result: foreach swaption [CalibPrice,BlackPrice,PercentDiff]" );
}

#endif // PARSE_INPUT

