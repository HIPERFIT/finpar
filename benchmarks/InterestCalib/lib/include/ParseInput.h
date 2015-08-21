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
                    real_t*& swaption_quotes,  // [num_swap_quotes,4]
                    
                    UINT&  num_hermitians,
                    real_t*& hermite_coefs,   // [num_hermite]
                    real_t*& hermite_weights,  // [num_hermite]
                
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

        if ( read_array(sizeof(real_t), read_real, (void**)&swaption_quotes, shape, 2) ) {
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

        if ( read_array(sizeof(real_t), read_real, (void**)&hermite_coefs, shape, 1) ) {
            fprintf(stderr, "Syntax error when reading the hermitian coeficients.\n");
            exit(1);
        }
        atr_ok = ( shape[0] == num_hermitians );
        assert(atr_ok && "Incorrect shape of the hermitian coefficients!");

        if ( read_array(sizeof(real_t), read_real, (void**)&hermite_weights, shape, 1) ) {
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

#endif // PARSE_INPUT

